function failed = tdrssmeas_test()
%
% tdrssmeas_test Regression test for tdrssmeas
% See also: tdrssmeas.m
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
%   Ravi Mathur         08/28/2012      Extracted from tdrssmeas.m

%gsList = createGroundStationList();
disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
tdrss.type  = 'Keplerian'; % Input TDRS States as Keplarian Elements
% TDRS East
tdrss.sat(1).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
tdrss.sat(1).sma  = 42165.3431; %semi-major axis
tdrss.sat(1).ecc  = 0.00026248; %eccentricity
tdrss.sat(1).incl = 1.7350*pi/180; %inclination
tdrss.sat(1).raan = 278.2107*pi/180; %right ascension of the ascending node
tdrss.sat(1).argp = 270.1285*pi/180; %argument of periapse
tdrss.sat(1).mean = 135.9127*pi/180; %mean anomaly
% TDRS West
tdrss.sat(2).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
tdrss.sat(2).sma  = 42166.4487; %semi-major axis
tdrss.sat(2).ecc  = 0.00030059; %eccentricity
tdrss.sat(2).incl = 8.5680*pi/180; %inclination
tdrss.sat(2).raan = 61.9693*pi/180; %right ascension of the ascending node
tdrss.sat(2).argp = 92.0025*pi/180; %argument of periapse
tdrss.sat(2).mean = 37.3454*pi/180; %mean anomaly
% TDRS Spare
tdrss.sat(3).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
tdrss.sat(3).sma  = 42167.6194; %semi-major axis
tdrss.sat(3).ecc  = 0.00025831; %eccentricity
tdrss.sat(3).incl = 9.8973*pi/180; %inclination
tdrss.sat(3).raan = 54.8622*pi/180; %right ascension of the ascending node
tdrss.sat(3).argp = 138.8066*pi/180; %argument of periapse
tdrss.sat(3).mean = 125.6643*pi/180; %mean anomaly
% Convert mean anomaly to true anomaly
for n=1:length(tdrss.sat)
    S.M = tdrss.sat(n).mean;
    S.ecc = tdrss.sat(n).ecc;
    tdrss.sat(n).tran = kepanom(S, 'nu'); %true anomaly
end

% Input TDRS Antenna Data
tdrss.sat(1).antenna.maxrange  = 70000; %km
tdrss.sat(1).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(2).antenna.maxrange  = 70000; %km
tdrss.sat(2).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(3).antenna.maxrange  = 70000; %km
tdrss.sat(3).antenna.halfangle = 13*pi/180; %radians

% Set propagator information
epoch  = datenum([2008  9 26  0  0   .000]); %UTC epoch
tspan  = 0:600:3600;
x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
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
measOptions = setOdtbxOptions(measOptions,'EarthAtmMaskRadius',6478.12); %km
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);

% Set Expected Results
y_expected = [
    77143.8676829857  2835.46917821229               NaN               NaN  77941.4862100122 -6249.96951871453;
    78228.5123284356 -20780.2992395048               NaN               NaN  79843.3003057668 -25420.8497100346;
    81412.4520066554 -32498.1098036057   81376.442501345   34456.437489037  83240.3532719398 -31819.6509120546;
    85144.6557209073 -31141.3680452223  77380.2041620424  33700.8409138789  86685.8437713931 -27168.8751418441;
    NaN               NaN  74065.7701163095  22895.1417907037               NaN               NaN;
    NaN               NaN  72406.4859646846  5497.55474322249               NaN               NaN;
    NaN               NaN  72840.7066213739 -12693.5108798617               NaN               NaN]';

% Display Expected Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Expected TDRSS Measurements (East, West, and Spare) by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','East Range','East Doppler','West Range','West Doppler','Spare Range','Spare Doppler');
fprintf('%4s %13s %14s %14s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:7
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),y_expected(:,n))
end
disp(' ')
disp(' ')

% Run tdrssmeas
[y,H,R] = tdrssmeas(t,x,measOptions); %#ok<NASGU>
% Even though I don't use H and R here, I still need to generate them to
% test the portion of tdrssmeas that calculates them. If they aren't
% requested as outputs, then the "nargout > 2" check won't be true.

% Display Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Calculated TDRSS Measurements (East, West, and Spare) by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','East Range','East Doppler','West Range','West Doppler','Spare Range','Spare Doppler');
fprintf('%4s %13s %14s %14s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:7
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),y(:,n))
end

passed = tol > max( max( abs( y - y_expected ) ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end

end
