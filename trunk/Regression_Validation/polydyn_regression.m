function f = polydyn_regression
% POLYDYN_REGRESSION Regression test for polydyn.
%
% F = POLYDYN_REGRESSION() runs the regression test
%
% keyword: measurement
% See also OPNAVMEAS, CAMERA, ATTITUDE, BODY
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Kenneth Getzandanner    05/24/2012      Original polydyn_example.m


% Set up polyhedron
X = [1 -.5 -.5 0 0;
    0 0.866 -0.866 0 0;
    0 0 0 1.732 -1.732];

x = X(1,:)';
y = X(2,:)';
z = X(3,:)';

ic = [mean(x);mean(y);mean(z)];

x = x-ic(1);
y = y-ic(2);
z = z-ic(3);
dt = DelaunayTri(x,y,z);

% Create TriRep object
[tri Xb] = freeBoundary(dt);
tr = TriRep(tri, Xb);

% Define body mass, volume, and density
volume = 1.4999;
rho = 5.1977e+11;

mass = rho*volume;
GM = mass*6.67300e-20;

% Initial spacecraft state
x_kep.sma  = 2;
x_kep.ecc  = 0.01;
x_kep.incl = 35*(pi/180);
x_kep.raan = 250*(pi/180);
x_kep.argp = 270*(pi/180);
x_kep.tran = 180*(pi/180);

x0 = kep2cart(x_kep,GM);

% Set up polydyn options
dynOpts.tr = tr;
dynOpts.rho = rho;

dynOpts.w = 2*pi/650;
dynOpts.PRA = 0;
dynOpts.DEC = pi/2;
dynOpts.RA = 0;

% Define simulation time space
tspan = [0 86400/2];

% Set integrator options
odeOpts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Propagate spacecraft
tic
[~,x] = ode113(@polydyn,tspan,x0,odeOpts,dynOpts);
toc

% Uncomment to generate regression test data
% x_test = x;
% save polydyn_data x_test

% Load test data
load polydyn_data

% Compare generated and test data
tol = 1e-9;

dx = abs(x-x_test);

if dx<tol
    f = 0;
else
    f = 1;
    fprintf('POLYDYN Regression test failed! dx = %g\n\n',max(max(dx)))
end

end
