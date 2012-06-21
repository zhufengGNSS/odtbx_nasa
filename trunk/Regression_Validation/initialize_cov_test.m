function f = initialize_cov_test
% INITIALIZE_COV_TEST Regression test for intialize_cov.
%
% F = INITIALIZE_COV_TEST runs the regression test
%
% keyword: covariance
% See also INITIALIZE_COV
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
%   Benjamin Asher          06/22/2012      Original


% Initial conditions
r=[4387.08;4806.18;-852.852];
v=[-5.41772;4.31046;-3.57805];
sig_pos = 1e-1;
sig_vel = 1e-4;
sig_rad = 1e-3;
sig_spd = 1e-6;
mu = 3.986004418e5;

% Generate test data
Po = initialize_cov(r,v,sig_pos,sig_vel,sig_rad,sig_spd,mu);

% Uncomment to generate regression test data
% Po_test = Po;
% save initialize_cov_data.mat Po_test

% Load test data
load initialize_cov_data.mat

% Compare generated and test data
tol = 1e-9;
dP = Po_test-Po;

if abs(dP)<tol
    f = 0;
else
    f = 1;
    fprintf('INITIALIZE_COV regression test failed! dP = %g\n\n',max(max(dP)))
end
end