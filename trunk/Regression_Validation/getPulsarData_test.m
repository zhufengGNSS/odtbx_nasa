function failed = getPulsarData_test

% Test utilizing the pulsar database to read and assign xnav data for use
% with xnavmeas.
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

% Sun Hur-Diaz
% Emergent Space Technologies, Inc.
% 9/23/2010
%
% 
%% Set the measurement options
measOptions = odtbxOptions('measurement');

% Create xnav options structure and add xnav specific fields
xnavOpts.obs_time = 1e5;
xnavOpts.epoch = datenum('11 Aug 2010 12:00:00');
xnavOpts.num_meas = 2;
xnavOpts.sigma_r = [];
xnavOpts.sigma_rr = [];
xnavOpts.useRangeRate = true;

pulsar_id = {'B1937+21','B0531+21'};
pulsar_database = 'pulsars.txt';
[rso_eci,C_1] = getPulsarData(pulsar_database,pulsar_id);

xnavOpts.rso_eci = rso_eci;
xnavOpts.C_1 = C_1;

measOptions = setOdtbxOptions(measOptions,'xnav',xnavOpts);

t = 0:1e5:2e5;
x0 = [5263998.56;-2252740.29;3541580.75;-336.973421;6257.216826;4465.898658]*1e-3; % km & km/sec
[t,x] = integ(@r2bp,t,x0,[],3e5);

[Y,H,R] = xnavmeas(t,x,measOptions);

% Compare with test data
testdata = load('getPulsarData_test_data');

failed = 0;
tol = 1e-12;

if any(any(testdata.rso_eci-rso_eci)) > tol
    failed = 1;
    disp('Failed getPulsarData_test in rso_eci comparison')
end
if any(any(testdata.C_1-C_1)) > tol
    failed = 1;
    disp('Failed getPulsarData_test in rso_eci comparison')
end
if any(any(testdata.Y-Y)) > tol
    failed = 1;
    disp('Failed getPulsarData_test in rso_eci comparison')
end
if any(any(testdata.H-H)) > tol
    failed = 1;
    disp('Failed getPulsarData_test in rso_eci comparison')
end
if any(any(testdata.R-R)) > tol
    failed = 1;
    disp('Failed getPulsarData_test in rso_eci comparison')
end





