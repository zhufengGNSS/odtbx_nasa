%% This file is an example for the tdrss measurement model, tdrssmeas.m
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

tic

%% Input TDRS States as Keplarian Elements

% Units assumed are km and radians
tdrss.type  = 'Keplerian';

% TDRS East
tdrss.sat(1).epoch = datenum([2008  9 26  0  0   .000]);
tdrss.sat(1).sma  = 42165.3431; %semi-major axis
tdrss.sat(1).ecc  = 0.00026248; %eccentricity
tdrss.sat(1).incl = 1.7350*pi/180; %inclination
tdrss.sat(1).raan = 278.2107*pi/180; %right ascension of the ascending node
tdrss.sat(1).argp = 270.1285*pi/180; %argument of periapse
tdrss.sat(1).mean = 135.9127*pi/180; %mean anomaly

% TDRS West
tdrss.sat(2).epoch = datenum([2008  9 26  0  0   .000]);
tdrss.sat(2).sma  = 42166.4487; %semi-major axis
tdrss.sat(2).ecc  = 0.00030059; %eccentricity
tdrss.sat(2).incl = 8.5680*pi/180; %inclination
tdrss.sat(2).raan = 61.9693*pi/180; %right ascension of the ascending node
tdrss.sat(2).argp = 92.0025*pi/180; %argument of periapse
tdrss.sat(2).mean = 37.3454*pi/180; %mean anomaly

% TDRS Spare
tdrss.sat(3).epoch = datenum([2008  9 26  0  0   .000]);
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

%% Input TDRS Antenna Data
tdrss.sat(1).antenna.maxrange  = 70000; %km
tdrss.sat(1).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(2).antenna.maxrange  = 70000; %km
tdrss.sat(2).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(3).antenna.maxrange  = 70000; %km
tdrss.sat(3).antenna.halfangle = 13*pi/180; %radians

%% Set propagator information
epoch  = datenum([2008  9 26  0  0   .000]);%datenum('1 Oct 2008 09:27:52.832');
tspan  = 0:10:360;
x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec

jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
    % default information for propagating Earth-centric orbits using JAT 
    % force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

%% Measurement Options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useLightTime', false);
measOptions = setOdtbxOptions(measOptions,'EarthAtmMaskRadius',6478.12); %km
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);

%% Run tdrssmeas
[y,H,R] = tdrssmeas(t,x,measOptions);
R
H
y

toc