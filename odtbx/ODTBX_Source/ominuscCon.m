function [dy,H,R,Pdy] = ominuscCon(datfun,t,x,y,options,P,datarg,isel)
% OMINUSCCON  Innovations, meas. partials, and possibly innovation covariance.
%   [DY,H,R] = OMINUSC(DATFUN,T,X,Y,OPTIONS) uses the supplied measurement
%   data function DATFUN to compute the difference between observed and
%   computed measurements, dy = y - h(t,x), known as the "innovation" if
%   the computed measurement is based on an "a priori" value of the state,
%   x, and known as the "residual" if the computed measurment is based on
%   the "a posteriori" value of the state (after it is updated with the
%   measurement(s)).  Function DATFUN must return a column vector
%   corresponding to h(t,x), and two additional outputs corresponding to
%   the measurement partials, H(t) = dh(t,x)/dx, and the measurement noise
%   covariance, R = E[vv'], where y = h(t,x) + v.  As an alternate to
%   supplying H(t), DATFUN may return an empty matrix as its second output,
%   in which case OMINUSC will numerically compute H(t) using NUMJAC.  If
%   DATFUN is "vectorized," then h(t,x) must return as its first output a
%   2-D array with each column corresponding to h(t(i),x(t(i)); its next
%   two outputs must be 3-D arrays with each "slice" corresponding to
%   H(t(i)) and R(t(i)), respectively.
%
%   [DY,H,R,PDY] = OMINUSC(DATFUN,T,X,Y,OPTIONS,P) will also compute the
%   covariance matrix associated with DY, if P, the error covariance of the
%   state, is also supplied.  If DATFUN is vectorized, then the additional
%   input P must be a 3-D array with each slice corresponding to P(t(i)),
%   and the additional output Pdy will be of the same form.  Use OPTIONS =
%   [] as a place holder if no options are set.
%
%   [DY,H,R,PDY] = OMINUSC(DATFUN,T,X,Y,OPTIONS,P,ISEL) returns only
%   those elements of the outputs selected by the index array ISEL.
%
%   [DY,H,R,PDY] = OMINUSC(DATFUN,T,X,Y,OPTIONS,P,ISEL,DATARG) passes
%   DATARG to DATFUN as DATFUN(T,X,DATARG).  
%
%   Examples (" means same as above)
%      Given [y,H,R] = range2d(t,x,sig), range2d not vectorized:
%         [dy,H,R] = ominusc(@range2d,t,x,y,[],[],[],sig) % range2d returns H
%         ["] = ominusc(") % range2d returns H=[]
%         [dy,H,R,Pdy] = ominusc(@range2d,t,x,y,[],P,[],sig) 
%      With opts = setOdtbxOptions('DatVectorized',1) and range2D vectorized:
%         [dy,H,R] = ominusc(@range2d,t,x,y,opts,[],[],sig)
%         ["] = ominusc(") % range2d returns H
%         ["] = ominusc(") % range2d returns H=[]
%         [dy,H,R,Pdy] = ominusc(@range2d,t,x,y,opts,P,[],sig)
%
%   keyword: Estimation,
%
%   See also
%      Numerical Jacobian:    NUMJAC
%      options handling:      ODTBXOPTIONS, SETODTBXOPTIONS,
%                             GETODTBXOPTIONS
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

% ODEAS Toolbox Spiral 0

% Russell Carpenter
% NASA Goddard Space Flight Center

% datfun
% size(t)
% size(x)
% size(y)
% options
% size(P)
% datarg

if nargin < 7,
    datarg = [];
end
mt = length(t);

[h,H,R] = feval(datfun,t,x,datarg);
if isempty(H),
    for k = mt:-1:1,
        H(:,:,k) = estjac(datfun,t(k),x(:,k),h(:,k),options,datarg);
    end
    H(isnan(H)) = 0;
end
dy = y - h;

if nargout == 4,
    if nargin < 6,
        error('OMINUSC:noStateCov',['Residual covariance output ', ...
            'requires state covariance input.'])
    end
    [m,n] = size(dy);
    Pdy = NaN(m,m,mt);
    ig = [];
    for k = 1:m
        ig=union(ig,find(~isnan(dy(k,:))));
    end
    for k = ig
        if nargin == 8
            ikeep = ~isnan(dy(:,k)) & isel(:,k);
        else
            ikeep = ~isnan(dy(:,k));
        end
%         size(ikeep)
%         size(k)
%         size(H)
%         disp('-')
        if numel(ikeep) < 1
            keyboard;
        end
        Hnn = H(ikeep(1),:,k(1));
        Rnn = R(ikeep(1),ikeep(1),k(1));
        Pdynn = Hnn*P(:,:,k(1))*Hnn' + Rnn;
        Pdy(ikeep(1),ikeep(1),k(1)) = Pdynn;
    end
end
if nargin == 8,
    dy = dy(isel,:);
    H = H(isel,:,:);
    R = R(isel,isel,:);
    Pdy = Pdy(isel,isel,:);
end

%-------------------------------------------------------------------------
function H = estjac(g,t,x,xdot,options,datarg)
% ESTJAC  Helper function for numerical Jacobian function evaluation.
%   Note this is slightly different from the ESTJAC subfunction that
%   INTEG uses, due to difference in the fields of the OPTIONS structure.
%   Also, keeping them separate keeps their persistent variables distinct.
%
%   See also NUMJAC, ODTBXOPTIONS

persistent FAC THRESH VECTRIZD JPATTRN JG
if length(THRESH)~=length(x)
    THRESH = [];
    FAC = [];
end
if isempty(THRESH),
    atol = getOdtbxOptions(options,'DatJTolerance',1e-6);
    THRESH = zeros(size(x))+ atol(:);
end
if isempty(VECTRIZD),
    VECTRIZD = getOdtbxOptions(options,'DatVectorized','off');
end
if isempty(JPATTRN),
    JPATTRN = getOdtbxOptions(options,'DatJPattern',[]);
end
[H,FAC,JG] = numjac(g,t,x,xdot,THRESH,FAC,VECTRIZD,JPATTRN,JG,datarg);