function failed = lnrmeas_test(type)
%
% lnrmeas_test Regression test for lnrmeas
% See also: lnrmeas.m
%
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%
%   Ravi Mathur         08/28/2012      Extracted from lnrmeas.m


disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
% Input Relay States as STKEphem type ephemeris files
relay.type            = 'STKEphem';
if strcmpi(type,'ValidationTest')
    relay.sat(1).filename = 'Relay1.e'; %Lunar Relay Satellite #1
    relay.sat(2).filename = 'Relay2.e'; %Lunar Relay Satellite #2
    if ~exist(relay.sat(1).filename)||~exist(relay.sat(2).filename)
        error('Validation test files could not be found.');
    end
elseif strcmpi(type,'RegressionTest')
    relay.sat(1).filename = 'Relay1_regress.e'; %Lunar Relay Satellite #1
    relay.sat(2).filename = 'Relay2_regress.e'; %Lunar Relay Satellite #2
    if ~exist(relay.sat(1).filename)||~exist(relay.sat(2).filename)
        msg = ['Regression test files could not be found. Do you have the ',...
            'regression test folder in your path? Running the validation ',...
            'test instead...'];
        warning('LNRMEAS:Regression1',msg)
        lnrmeas_test('ValidationTest');
        msg = ['Regression test could not be ran because the regression ',...
            'ephems could not be found. The validation case was ran ',...
            'instead. The regression test is being considered a failure ',...
            'regardless of validation case results!'];
        warning('LNRMEAS:Regression2',msg)
        failed = 1;
        return
    end
end
    
% Set the radius of Moon Obscuration
relay.MoonMaskRadius = 20+JATConstant('meanRadius','Moon')/1000;
%Using the Mean Radius of the Moon plus 20 km to account for lunar mountains

% Input relay antenna data
relay.sat(1).antenna.maxrange  = 100000; %km
relay.sat(2).antenna.maxrange  = 100000; %km

% Set propagator information
epoch  = datenum('1 Dec 2016 00:00:00.000'); %UTC
tspan  = 0:600:7200;
x0     = [-3.10570021720167e+004
    -3.80943874664242e+005
    -1.28680660459686e+005
    1.05814791201683e+000
    1.44219161241440e+000
    6.24467648736737e-001]; %in km & km/sec

jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %UTC in datenum format
jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM2');
jOptions = setOdtbxOptions(jOptions, 'gravDeg', 2, 'gravOrder', 0); %max 20x20
jOptions = setOdtbxOptions(jOptions, 'useSolarGravity', true);
jOptions = setOdtbxOptions(jOptions, 'useLunarGravity', true);
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
    % default information for propagating Earth-centric orbits using JAT 
    % force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

% Measurement Options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %UTC
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useLightTime', false);
measOptions = setOdtbxOptions(measOptions,'relay', relay);

% Load Expected Results
% NOTE: These results are based on the most recent July 2012 leap second. If a
% more recent leap second is added, then these results must be updated.
load lnrmeas_expected.mat;

% Display Expected Results
fprintf('%s\n',char(ones(1,66)*'-'));
disp('Expected Relay Measurements by time:')
fprintf('%s\n',char(ones(1,66)*'-'));
fprintf('%-6s %14s %14s %14s %14s\n','Time','#1 Range','#1 Doppler','#2 Range','#2 Doppler');
fprintf('%4s %13s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,66)*'-'));
for n=1:13
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f\n',t(n),y_expected(:,n))
end         
disp(' ')
disp(' ')

% Run lnrmeas
[y,H,R] = lnrmeas(t,x,measOptions);

% Display Results
fprintf('%s\n',char(ones(1,66)*'-'));
disp('Calculated Relay Measurements by time:')
fprintf('%s\n',char(ones(1,66)*'-'));
fprintf('%-6s %14s %14s %14s %14s\n','Time','#1 Range','#1 Doppler','#2 Range','#2 Doppler');
fprintf('%4s %13s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,66)*'-'));
for n=1:13
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f\n',t(n),y(:,n))
    %fprintf('%-6i %20.11f %20.11f %20.11f %20.11f\n',t(n),y(:,n))
end     

passed = tol > max( max( abs( y - y_expected ) ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed! Please verify that the most recent leap second is July 2012.')
else
    disp(' ')
    disp('Test Passed.')
end

end
