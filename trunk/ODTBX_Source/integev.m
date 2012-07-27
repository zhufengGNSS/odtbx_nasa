function [t,x,te,xe,Phi,Phie,S,Se] = integev(dynfun,tspan,x0,options,dynarg,evtfun)
% INTEGEV  Integration of state, and possibly STM, etc. w/event handling.
%   [T,X] = INTEGEV(DYNFUN,TSPAN,X0,OPTIONS) uses the supplied dynamics
%   function DYNFUN to integrate the ordinary differential equation,
%   dx(t)/dt = f(t,x), with x(t0) = x0.  Function DYNFUN(T,X) must return a
%   column vector corresponding to f(t,x).  If DYNFUN is "vectorized," then
%   f(t,x) must be a 2-D array with each column corresponding to
%   f(t(i),x(t(i)).  
%
%   [T,X,...] = INTEGEV(DYNFUN,TSPAN,X0,OPTIONS,DYNARG) passes DYNARG to
%   DYNFUN as DYNFUN(T,X,DYNARG).  Use OPTIONS = [] as a place holder if no
%   options are set.
%
%   [T,X,TE,XE] = INTEGEV(DYNFUN,TSPAN,X0,OPTIONS,DYNARG,EVTFUN) will detect
%   events specifed by EVTFUN and return the event times and states in TE
%   and XE.  Event functions have the form 
%      [VALUE,ISTERMINAL,DIRECTION] = EVTFUN(T,X)
%   where VALUE, ISTERMINAL, and DIRECTION are vectors for which the ith
%   element corresponds to the ith event function:
%      VALUE(i) is the value of the ith event function.
%      ISTERMINAL(i) = 1, if the integration stops at a zero of the function, 
%                    = 0, otherwise.
%      DIRECTION(i) = 0 if all zeros are to be located (the default), 
%                   = +1 if only zeros where the event function is increasing, 
%                   = -1 if only zeros where the event function is decreasing. 
%
%   [T,X,TE,XE,PHI,PHIE] = INTEGEV(DYNFUN,TSPAN,X0,OPTIONS) will also
%   integrate the differential equation of the state transition matrix,
%   d/dt[Phi'(t,t0)] = A(t) Phi(t,t0), where A(t) = df(t,x)/dx.  Function
%   DYNFUN(T,X) must return an additional output if called with two output
%   arguments, which may either be a matrix corresponding to A(t), or else
%   an empty matrix, in which case INTEGEV will numerically compute A(t)
%   using NUMJAC. If A(t) is supplied, it must be a 3-D array for the
%   vectorized case, with each "slice" corresponding to A(t(i)).
%
%   [T,X,TE,XE,PHI,PHIE,S,SE] = INTEGEV(DYNFUN,TSPAN,X0,OPTIONS) will also
%   integrate the differential equation of the discrete process noise
%   matrix, dS(t)/dt = A(t)S(t) + S(t)A'(t) + Q(t), if DYNFUN returns the
%   process noise spectral density matrix, Q = E[ww'], where x' = f(t,x) +
%   w, as an additional output.  If DYNFUN is vectorized, then Q must be a
%   3-D array, with each "slice" corresponding to Q(t(i)).
%
%   See also
%      INTEG
%      Numerical Jacobian:    NUMJAC
%      ODE solvers:           ODE113, ODE23, ODE45
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
% Date: 2005-02-22 12:07:47 -0500 (Tue, 22 Feb 2005)

nx = length(x0);
nx2 = nx^2;
if nargin < 4 || isempty(options)
    options = odtbxOptions;
end
if nargin == 6
    options = setOdtbxOptions(options,'OdeSolvOpts',...
        odeset(options.OdeSolvOpts,'events',evtfun));
end
if nargout > 6 % Output process noise
    z0(nx+nx2+1:nx+2*nx2,:) = zeros(nx2,1);
end
if nargout > 4 % Output STM
    z0(nx+1:nx+nx2,:) = reshape(eye(nx),nx2,1);
end
z0(1:nx,:) = x0;
odesolv = getOdtbxOptions(options,'OdeSolver',@ode113);
odeopts = getOdtbxOptions(options,'OdeSolvOpts',...
    odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10)); 
if nargin == 6
    [t,z,te,ze] = feval(odesolv,@odefun,tspan,z0,odeopts,options,dynfun,nx,dynarg);
else
    [t,z] = feval(odesolv,@odefun,tspan,z0,odeopts,options,dynfun,nx,dynarg);
end
x = z(:,1:nx)';
if nargout > 4 % STM
    nt = length(t);
    Phi = reshape(z(:,nx+1:nx+nx2)',nx,nx,nt);
end
if nargout > 6 % Process noise
    S = reshape(z(:,nx+nx2+1:end)',nx,nx,nt);
end
if nargin == 6 && ~isempty(ze)
    xe = ze(:,1:nx)';
    if nargout > 5
        nte = length(te);
        Phie = reshape(ze(:,nx+1:nx+nx2)',nx,nx,nte);
    end
    if nargout > 7
        Se = reshape(ze(:,nx+nx2+1:end)',nx,nx,nte);
    end
else
    te = [];
    xe = [];
    Phie = [];
    Se = [];
end

%-------------------------------------------------------------------------
function zdot = odefun(t,z,options,dynfun,nx,dynarg)
% ODEFUN  Derivative function for ODE Solvers

% FYI: nx = (realsqrt(1+8*nz)-1)/4;

[nz,mz] = size(z);
nx2 = nx^2;
x = z(1:nx,:);
if nz > nx,
    Phi = reshape(z(nx+(1:nx2),:),nx,nx,mz);
    [xdot,A,Q] = feval(dynfun,t,x,dynarg);
    if isempty(A),
        for k = mz:-1:1,
            A(:,:,k) = estjac(dynfun,t(k),x(:,k),xdot(:,k),options,dynarg);
        end
    end
    if nz > nx+nx2 && ~isempty(Q),
        S = reshape(z(nx2+nx+1:end,:),nx,nx,mz);
        for k = mz:-1:1,
            Sdot = A(:,:,k)*S(:,:,k) + S(:,:,k)*A(:,:,k)' + Q(:,:,k);
            zdot(nx2+nx+1:2*nx2+nx,k) = Sdot(:);
        end
    end
    for k = mz:-1:1,
        Phidot = A(:,:,k)*Phi;
        zdot(nx+(1:nx2),k) = Phidot(:);
    end
else
    xdot = feval(dynfun,t,x,dynarg);
end
zdot(1:nx,:) = xdot;

%-------------------------------------------------------------------------
function A = estjac(f,t,x,xdot,options,dynarg)
% ESTJAC  Helper function for numerical Jacobian function evaluation.
%   Note this is slightly different from the ESTJAC subfunction that
%   OMINUSC uses, due to difference in the fields of the OPTIONS structure.
%   Also, keeping them separate keeps their persistent variables distinct.
%
%   See also NUMJAC

persistent FAC THRESH VECTRIZD JPATTRN JG
if length(THRESH)~=length(x)
    THRESH = [];
    FAC = [];
end
if isempty(options) || ~isfield('OdeSolvOpts',options)
    options.OdeSolvOpts = [];
end
if isempty(THRESH),
    atol = odeget(options.OdeSolvOpts,'AbsTol',1e-6,'fast');
    THRESH = zeros(size(x))+ atol(:);
end
if isempty(VECTRIZD),
    VECTRIZD = odeget(options.OdeSolvOpts,'Vectorized','off','fast');
end
if isempty(JPATTRN),
    JPATTRN = odeget(options.OdeSolvOpts,'Jpattern',[],'fast');
end
[A,FAC,JG] = numjac(f,t,x,xdot,THRESH,FAC,VECTRIZD,JPATTRN,JG,dynarg);