function [t,x,te,xe,Phi,S] = integ(dynfun,tspan,x0,options,dynarg)
% INTEG  Integration of state, and possibly STM and process noise.
%   [T,X] = INTEG(DYNFUN,TSPAN,X0,OPTIONS) uses the supplied dynamics
%   function DYNFUN to integrate the ordinary differential equation,
%   dx(t)/dt = f(t,x), with x(t0) = x0.  Function DYNFUN(T,X) must return a
%   column vector corresponding to f(t,x).  If DYNFUN is "vectorized," then
%   f(t,x) must be a 2-D array with each column corresponding to
%   f(t(i),x(t(i)).  
%
%   [T,X,PHI] = INTEG(DYNFUN,TSPAN,X0,OPTIONS) will also integrate the
%   differential equation of the state transition matrix, d/dt[Phi'(t,t0)]
%   = A(t) Phi(t,t0), where A(t) = df(t,x)/dx.  Function DYNFUN(T,X) must
%   return an additional output if called with two output arguments, which
%   may either be a matrix corresponding to A(t), or else an empty matrix,
%   in which case INTEG will numerically compute A(t) using NUMJAC. If
%   A(t) is supplied, it must be a 3-D array for the vectorized case, with
%   each "slice" corresponding to A(t(i)).
%
%   [T,X,PHI,S] = INTEG(DYNFUN,TSPAN,X0,OPTIONS) will also integrate the
%   differential equation of the discrete process noise matrix, dS(t)/dt =
%   A(t)S(t) + S(t)A'(t) + Q(t), if DYNFUN returns the process noise
%   spectral density matrix, Q = E[ww'], where x' = f(t,x) + w, as an
%   additional output.  If DYNFUN is vectorized, then Q must be a 3-D
%   array, with each "slice" corresponding to Q(t(i)).
%
%   [T,X,...] = INTEG(DYNFUN,TSPAN,X0,OPTIONS,DYNARG) passes DYNARG to
%   DYNFUN as DYNFUN(T,X,DYNARG).  Use OPTIONS = [] as a place holder if no
%   options are set.
%
%   Examples (" means same as above)
%      Given xdot = pr2bp(t,x,mu), pr2bp not vectorized
%         [t,x] = integ(@pr2bp,tspan,x0,[],mu)
%      Given [xdot,A] = pr2bp("):
%         [t,x,Phi] = integ(") % pr2bp returns A
%         [t,x,Phi] = integ(") % pr2bp returns A=[]
%      Given [xdot,A,Q] = pr2bp("):
%         [t,x,Phi,S] = integ(") % pr2bp returns Q
%      With opts = setOdtbxOptions('OdeSolvOpts',odeset('Vectorized','on')) and
%      pr2bp vectorized:
%         [t,x] = integ(@pr2bp,tspan,x0,opts,mu) 
%         [t,x,Phi] = integ(") % pr2bp returns A
%         [t,x,Phi] = integ(") % pr2bp returns A=[]
%         [t,x,Phi,S] = integ(") % pr2bp returns Q
% 
%   keyword: Integrators, 
%
%   See also
%      INTEGEV
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

switch nargout
    case 2
        [t,x] = integev(dynfun,tspan,x0,options,dynarg);
    case 3
        [t,x,~,~,Phi] = integev(dynfun,tspan,x0,options,dynarg);
    case 4
        [t,x,~,~,Phi,~,S] = integev(dynfun,tspan,x0,options,dynarg);
    case 5
        [t,x,te,xe,Phi] = integev(dynfun,tspan,x0,options,dynarg);
    case 6
        [t,x,te,xe,Phi,~,S] = integev(dynfun,tspan,x0,options,dynarg);
end