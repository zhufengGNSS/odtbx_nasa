%% This file is an example for the lunar relay measurement model, lnrmeas.m
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

%% Input Relay States as STKEphem type ephemeris files
relay.type            = 'STKEphem';
relay.sat(1).filename = 'Relay1.e'; %Lunar Relay Satellite #1
relay.sat(2).filename = 'Relay2.e'; %Lunar Relay Satellite #2

%% Set the radius of Moon Obscuration
relay.MoonMaskRadius = 20+JATConstant('meanRadius','Moon')/1000;
%Using the Mean Radius of the Moon plus 20 km to account for lunar mountains

%% Input relay antenna data
relay.sat(1).antenna.maxrange  = 100000; %km
relay.sat(2).antenna.maxrange  = 100000; %km

%% Set propagator information
epoch  = datenum('1 Dec 2016 00:00:00.000');
tspan  = 0:60:7200;
x0     = [-3.10570021720167e+004
    -3.80943874664242e+005
    -1.28680660459686e+005
    1.05814791201683e+000
    1.44219161241440e+000
    6.24467648736737e-001]; %in km & km/sec

jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
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

%% Measurement Options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useLightTime', false);
measOptions = setOdtbxOptions(measOptions,'relay', relay);

%%
[y,H,R] = lnrmeas(t,x,measOptions);
R
H
y

toc