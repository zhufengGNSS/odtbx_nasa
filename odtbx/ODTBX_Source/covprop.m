function [P,S,Phi,Gam] = covprop(dynfun,tspan,P0,X0,options,dynarg)
% COVPROP  Discrete propagation of covariance matrix.
%   COVPROP(DYNFUN,TSPAN,P) propagates a covariance matrix, P, over the
%   intervals specified in TSPAN, using the dynamics model specified by
%   DYNFUN, via the method given by van Loan ('78).
%
%   [P,S,PHI,GAM] = COVPROP(DYNFUN,TSPAN,P) also returns the process noise
%   covariance matrix, S, if any, that was used in the covariance
%   propagation, the state transition matrix that was used, and the
%   corresponding discrete noise input matrix.
%
%   [P,S,PHI] = COVPROP(DYNFUN,TSPAN,X,OPTIONS,DYNARG) allows for the input
%   of additional data that may be required, including a state vector, X,
%   in case DYNFUN is nonlinear function, an options struction, and any
%   additional arguments, DYNARG, that might be required for DYNFUN.
%
%   DYNFUN must support the calling syntax [XDOT,A,Q] = DYNFUN(T,X,DYNARG),
%   where xdot = f(t,x), A(t) = df(t,x)/dx.  If there is a process noise
%   input to the differential equation, w, then Q = E[ww'] should be the
%   power spectral density of the noise input.
%
%   Note: if f(t,x) is a nonlinear function, then the method of van Loan is
%   at best approximate, and then its accuracy will decline if the time
%   steps in TSPAN are not "small" compared with the modes of the dynamics.
%
%   Example
%      mu = 3.986e5;
%      a = 6878;
%      v = 8.339;
%      tspan = [0 100];
%      X0 = [a;0;0;v];
%      P0 = diag([1 1 1e-3 1e-3].^2);
%      tic
%      [t,X,Phi_i,S_i] = integ(@pr2bp,tspan,X0,[],mu); % pr2bp is on path
%      P_i = Phi_i(:,:,end)*P0*Phi_i(:,:,end)' + S_i(:,:,end);
%      toc
%      tic
%      [P,S,Phi,X] = covprop(@pr2bp,tspan,P0,X(:,[1 end]),[],mu);
%      toc
%      corrmat(P_i)
%      corrmat(P(:,:,end))
%
%   See also
%      integ
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

% OD Toolbox

% Russell Carpenter
% NASA Goddard Space Flight Center

[Xdot,A,Q] = feval(dynfun,tspan,X0,dynarg);
n = size(Xdot,1);
O = zeros(n,n);
m = length(tspan);
P = NaN([size(P0),m]);
P(:,:,1) = P0;
S = NaN([size(P0),m-1]);
Phi = NaN([size(P0),m-1]);
if nargout > 3,
    Gam = S;
end
for k = 1:m-1,
    dt = tspan(k+1) - tspan(k);
    F = [-A(:,:,k) Q(:,:,k); O A(:,:,k)'];
    Psi = expm(F*dt);
    Phi(:,:,k) = Psi(n+1:2*n,n+1:2*n)';
    S(:,:,k) = Phi(:,:,k)*Psi(1:n,n+1:2*n);
    S(:,:,k) = (S(:,:,k) + S(:,:,k)')/2; % Added 6/26/2008
    P(:,:,k+1) = Phi(:,:,k)*P(:,:,k)*Phi(:,:,k)' + S(:,:,k);
    if exist('Gam','var'),
        C = chol(S(:,:,k));
        W = diag(diag(C));
        Gam(:,:,k) = C'/W;
    end
end
P(:,:,1) = [];
