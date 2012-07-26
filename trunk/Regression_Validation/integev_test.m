function fail = integev_test
% Test function for integev, i.e. integ with events.
%
% Test case is a random run, which has an analytical solution:
%
% Phi = [1 dt; 0 1];
% Qd = q*[dt^3/3 dt^2/2; dt^2/2 dt];
%
% Hence, an initial condition of x0 = [dt; -1] will reach E[x(1)] = 0 at
% time dt.  The test event function checks for this.
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
% (for additional changes, see the svn repository)

q = 1;
tspan = [0 30];
x0 = [30; -1];

[~,x,Phi,Qd,te,xe,Phie,Qde] = integev(@dynfun,tspan,x0,[],q,@evtfun);
dt = tspan(end) - tspan(1);
Phichk = [1 dt; 0 1];
Qdchk = q*[dt^3/3 dt^2/2; dt^2/2 dt];
xchk = Phichk*x0;
dx = x(:,end) - xchk;
dPhi = Phi(:,:,end) - Phichk;
dQd = Qd(:,:,end) - Qdchk;
dxe = xe - xchk;
dPhie = Phie - Phichk;
dQde = Qde - Qdchk;
dte = te - x0(1);

fail = (norm(dx) > epstol(xchk)) || (norm(dPhi) > epstol(Phichk)) || ...
    (norm(dQd) > epstol(Qdchk)) || (norm(dxe) > epstol(xchk)) || ...
    (norm(dPhie) > epstol(Phichk)) || (norm(dQde) > epstol(Qdchk)) || ...
    (dte > epstol(te));

function [xdot,A,Q] = dynfun(~,x,q)
% Random run, i.e. a bias whose drift is a random walk.
A = [0, 1; 0, 0];
Q = diag([0,q]);
xdot = A*x;

function [value,isterminal,direction] = evtfun(~,x,varargin)
% Locate the time when x(1) passes through zero in a 
% decreasing direction and stop integration.
value = x(1);     % Detect x(1) = 0
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only

function tol = epstol(A)
% This tolerance is used in Mathworks' rank check.
tol = max(size(A)) * eps(norm(A));
