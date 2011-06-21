function [M,N] = observ(dynfun,datfun,tspan,x0,options,P0,y,dynarg,datarg)
% OBSERV  Observability grammian, and possibly mapped innovations.
%   OBSERV(DYNFUN,DATFUN,TSPAN,X0,OPTIONS) returns the observability
%   grammian corresponding to the dynamics specified by DYNFUN with initial
%   condition X0, and the observations specified by DATFUN, over the time
%   interval TSPAN.  OBSERV calls ESTINT to integrate the dynamics, and
%   calls ESTINV to gets the measurement partials and any data weighting.
%   Refer to help for ESTINT and ESTINV for information about the
%   requirements on DYNFUN and DATFUN.  See ODTBXOPTIONS for information
%   on OPTIONS.
%
%   OBSERV(DYNFUN,DATFUN,TSPAN,X0,OPTIONS,P0) adds the inverse of the
%   matrix P0, which represents any "a priori" information about the
%   initial state, to the grammian.  To avoid including OPTIONS, use
%   OPTIONS = [].
%
%   [M,N] = OBSERV(DYNFUN,DATFUN,TSPAN,X0,OPTIONS,P0,Y) returns the
%   grammian matrix in M, and the mapping of the observation residuals to
%   the initial time in N, so that a least squares correction to the
%   x0 may be formed from inv(M)*N.  Use [] as a place holder for
%   OPTIONS and/or P0.
%
%   [M,N] = OBSERV(DYNFUN,DATFUN,TSPAN,X0,OPTIONS,P0,Y,DYNARG,DATARG)
%   passes DYNARG to DYNFUN and DATARG to DATFUN as DYNFUN(T,X,DYNARG) and
%   DATFUN(T,X,DATARG), respectively.  Use [] as a place holder for
%   OPTIONS, P0, and/or Y.
%
%   Examples
%      Given xdot = pr2bp(t,x,mu), [y,H,R] = range2d(t,x,sig):
%         observ(@pr2bp,@range2d,tspan,x0,[],[],[],mu,sig)
%         observ(@pr2bp,@range2d,tspan,x0,[],P0,[],mu,sig)
%         [M,N] = observ(@pr2bp,@range2d,tspan,x0,[],[],y,mu,sig)
%         [M,N] = observ(@pr2bp,@range2d,tspan,x0,[],P0,y,mu,sig)
%      With opts = setOdtbxOptions('OdeSolver',@ode45,...):
%        observ(@pr2bp,@range2d,tspan,x0,opts,[],[],mu,sig)
%
%   keyword: Estimation,
%
%   See also
%      ODEAS dynamics:        INTEG
%      ODEAS observations:    OMINUSC
%      Numerical Jacobian:    NUMJAC
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
% Russell Carpenter: Original
% Kevin Berry: Added isnan(y) check so only the visible measurements are summed
% Derek Surka: Added loop over all measurements, corrected initialization
%              of N, and added isemtpy(R) check.
% Keith Speckman: Changed isnan(y) check to isnan(dy)
%   (for other changes, see the svn repository)

m = length(tspan)-1;
[t,x,Phi] = integ(dynfun,tspan,x0,options,dynarg);
if m == 1,
    t = t(end);
    x = x(:,end);
    Phi = Phi(:,:,end);
end
if isempty(y) | nargin < 7,
    y = 0;
end
[dy,H,R] = ominusc(datfun,t,x,y,options,[],datarg);
if isempty(P0) | nargin < 6,
    M = 0;
else,
    M = inv(P0);
end
N = zeros(size(M,1),1);

nMeas = size(y,1);
if( isempty(R) )              
    R = repmat(eye(nMeas),[1 1 m+1]);
end

for j = 1:nMeas
    for k = find(~isnan(dy(j,:))),   
        HPhi = H(j,:,k)*Phi(:,:,k);
        HPhiTW = HPhi'/R(j,j,k);
        M = M + HPhiTW*HPhi;
        if nargout == 2,
            dz=dy(j,k);
            N = N + HPhiTW*dz;
        end
    end
end
end
