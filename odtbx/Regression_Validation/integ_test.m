function fail = integ_test(odesolv)
% Test function for integ, using any ODTBX-supported integrator.
%
% Test case is a random run, which has an analytical solution:
%
% Phi = [1 dt; 0 1];
% Qd = q*[dt^3/3 dt^2/2; dt^2/2 dt];
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% ODTBX: Orbit Determination Toolbox
%
% Copyright (c) 2003-2012 United States Government as represented by the
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

% Version History:
% Original - Russell Carpenter
% (for additional changes, see the git repository)

q = 1;
tspan = [0 30];
x0 = [1; 1];

% Set the integrator if specified
options = [];
if (nargin >= 1)
    options = setOdtbxOptions('OdeSolver', odesolv);
end

[~,x,Phi,Qd] = integ(@dynfun,tspan,x0,options,q);
dt = tspan(end) - tspan(1);
Phichk = [1 dt; 0 1];
Qdchk = q*[dt^3/3 dt^2/2; dt^2/2 dt];
xchk = Phichk*x0;
dx = x(:,end) - xchk;
dPhi = Phi(:,:,end) - Phichk;
dQd = Qd(:,:,end) - Qdchk;

fail = (norm(dx) > epstol(xchk)) || (norm(dPhi) > epstol(Phichk)) || ...
    (norm(dQd) > epstol(Qdchk));

function [xdot,A,Q] = dynfun(~,x,q)
% Random run, i.e. a bias whose drift is a random walk.
A = [0, 1; 0, 0];
Q = diag([0,q]);
xdot = A*x;

function tol = epstol(A)
% This tolerance is used in Mathworks' rank check.
tol = max(size(A)) * eps(norm(A));
