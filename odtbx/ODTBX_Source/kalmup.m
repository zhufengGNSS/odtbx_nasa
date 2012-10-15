function [x,P,eflag,dy,Pdy,K] = kalmup(datfun,t,x,P,y,options,eflag,...
    eratio,datarg,isel,S,C)
%KALMUP Kalman update with sigma-edit.
% KALMUP  Kalman update with sigma-edit.
%   [X,P] = KALMUP(DATFUN,T,X,P,Y) updates the state and covariance X and P
%   respectively with the observation Y using Kalman's gain matrix.  The
%   function OMINUSC computes the innovation and its covariance using
%   DATFUN. If Y (and the output of DATFUN) is a vector, KALMUP performs a
%   vector update.
%
%   [X,P] = KALMUP(DATFUN,T,X,P,Y,OPTIONS) performs as above with default
%   properties replaced by values in OPTIONS, an argument created with the
%   SETODTBXOPTIONS function.  See ODTBXOPTIONS for details.
%
%   [X,P,EFLAG] = KALMUP(DATFUN,T,X,P,Y,OPTIONS,EFLAG,ERATIO) performs as
%   above, except that the measurements are processed according to EFLAG as
%   follows:
%      EFLAG  PROCEDURE
%        2    Force the measurement to be processed, regardless of ERATIO.
%        1    Accept the measurement, if the innovations ratio < ERATIO.
%   In the latter case, a rejection of the measurement due to a failure of
%   the ERATIO test is reported by setting EFLAG to 0.  Note that ERATIO
%   is compared to the quadratic form dy'*inv(Pdy)*dy, so for example a
%   "3-sigma edit" requires ERATIO = 9.
%
%   If eflag and eratio are both empty ([]), then measurements are edited
%   with the default eratio value of 9.
%
%   [X,P,EFLAG] = KALMUP(DATFUN,T,X,P,Y,OPTIONS,EFLAG,ERATIO,DATARG)
%   passes DATARG to DATFUN as DATFUN(T,X,DATARG).
%
%   [X,P,EFLAG] = KALMUP(DATFUN,T,X,P,Y,OPTIONS,EFLAG,ERATIO,DATARG,ISEL)
%   processes only those measurements selected by the index array ISEL. Use
%   [] as a placeholder for OPTIONS, EFLAG, ERATIO, DATARG, and/or ISEL.
%
%   [X,P,EFLAG] = KALMUP(DATFUN,T,X,P,Y,OPTIONS,EFLAG,ERATIO,DATARG,ISEL,S,C)
%   passes the solve-for and consider selction matrices such that xs=S*x,
%   xc=C*x, [S;C] is full rank, and Schmidt Kalman filter is implemented
%   whereby the full covariances are propagated and used in computing the gain,
%   but only the solve-for states are updated.
%
%   [X,P,EFLAG,DY] = KALMUP(DATFUN,T,X,P,Y,...) also returns the
%   innovations DY, i.e. the observed minus computed measurement prior to
%   the state update.  To get the residual, i.e. the observed minus
%   computed measurement after the state update, issue the command
%   OMINUSC(DATFUN,T,X,Y,OPTIONS) after issuing KALMUP.
%
%   [X,P,EFLAG,DY,PDY] = KALMUP(DATFUN,T,X,P,Y,...) also returns the
%   covariance matrix associated with DY.
%
%   [X,P,EFLAG,DY,PDY,K] = KALMUP(DATFUN,T,X,P,Y,...) also returns the
%   Kalman gain matrix.
%
%   keyword: Estimation,
%   See also
%      ODEAS observations:    OMINUSC
%
%   Examples
%      Given [y,H,R]=range2D(t,x,sig), and an observation, y:
%         [x,P] = kalmup(@range2D,t,x,P,y,[],[],[],sig)
%      With eflag = 1, eratio = 9:
%         [x,P,eflag] = kalmup(@range2D,t,x,P,y,eflag,eratio,[],sig)
%      If y is a vector, force a scalar update:
%         for k=1:length(y),
%            [x,P,eflag] = kalmup(@range2D,t,x,P,y(k),eflag,eratio,[],sig)
%         end
%      With eflag=[1;2;1;...]
%         [x,P,eflag] = kalmup(@range2D,t,x,P,y,eflag,eratio,[],sig)
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2011 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

% Russell Carpenter
% NASA Goddard Space Flight Center
% Modification History
% ---------------------
% 2009/01/13 Sun Hur-Diaz   Added Schmidt-Kalman filter capability
% 2009/09/26 Sun Hur-Diaz   Added a check for NaN in measurements and set
%                           corresponding elements of the Kalman gain to 
%                           zero so that the covariance matrix is not
%                           updated
% 2009/10/29 Kevin Berry    Corrected thresh calculation in the measurement
%                           editting section

[m,n] = size(y);
if n > 1,
    error('KALMUP:matrixMeas','Matrix measurements not supported.');
end
if n < 1,
    error('KALMUP:emptyMeas','Empty measurements not supported.');
end
%TODO: could do more consistency checking, e.g. dim(x) = dim(P), etc.
if nargin < 6,
    options = [];
end

process = true(m,1); %TODO: use true and false in options also
if(exist('isel', 'var'))
    if ~isempty(isel),
        process = process(isel);
        if length(eflag) > 1
            eflag = eflag(isel);
            eratio = eratio(isel);
        end
    end
else
    isel = (1:length(y));
end
nx=length(x);

if nargin < 11
    S=eye(nx);
    ns=nx;
else
    [ns,m2,~]=size(S);
    if m2 ~= nx
        error('Dimension 2 of S is inconsistent with length of x')
    end
end

if ns < nx
    if nargin < 12
        error('Consider selection matrix C must be specified.')
    else
        [nc,m2,~]=size(C);
        if m2 ~= nx
            error('Dimension 2 of C is inconsistent with length of x')
        elseif nc+ns ~= nx
            error('[S;C] must be full rank.')
        end
    end
elseif nargin >= 12 && ~isempty(C)
    error('C must be empty since S is full rank')
elseif ns == nx
    C=[];
end

Minv = inv([S;C]);
Stilde = Minv(:,1:ns);
Ctilde = Minv(:,ns+1:nx);

P = (P+P')/2;
if(exist('datarg', 'var'))
    [dy,H,R,Pdy] = ominusc(datfun,t,x,y,options,P,datarg,isel);
else
    [dy,H,R,Pdy] = ominusc(datfun,t,x,y,options,P,[],isel);
end

% TODO: Need to add something like the following for consider:
% Rhat = R{1};
% R = R{2};


% Check for NaNs in the measurement that are selected for processing
process(isnan(dy)) = false;

% Perform measurement editing
[process,eflag] = editmeas(process,dy,Pdy,eflag,eratio);

ikeep = find(process == true);

if ~isempty(ikeep)
    
    H = H(ikeep,:);
    R = R(ikeep,ikeep);
    
    PHt = P*H';
    
    K = S*PHt/Pdy(ikeep,ikeep);
    
    % Solved-for states
    xs = S*x;
    % Apply the gain on the solved-for states only
    xs = xs + K*dy(ikeep);
    % Combine with the consider states that are unchanged
    if ns < nx
        x = Stilde*xs + Ctilde*C*x;
    else
        x=xs;
    end
    ImKH = eye(length(x)) - Stilde*K*H;
    P = ImKH*P*ImKH' + Stilde*K*R*K'*Stilde';
    P = (P+P')/2;
    % TODO: Need to replace above with the following for consider (possibly
    % do this in a new update function):
    % Hs = H*Stilde;
    % ImKH = eye(length(s)) - K*Hs;
    % Phatssa = ImKH*Phatssa*ImKH';
    % Phatssv = ImKH*Phatssv*ImKH' + K*Rhat*K';
    % Phatssw = ImKH*Phatssw*ImKH';
    % ImSKH = eye(length(x)) - Stilde*K*H;
    % Pa = ImSKH*Pa*ImSKH';
    % Pv = ImSKH*Pv*ImSKH' + Stilde*K*R*K'*Stilde';
    % Pw = ImSKH*Pw*ImSKH';
    % Phatss = Phatssa + Phatssv + Phatssw;
    % dPssa = S*Pa*S' - Phatssa;
    % dPssv = S*Pv*S' - Phatssv;
end

