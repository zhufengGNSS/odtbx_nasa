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

%% Input TDRS States as SPEphem type ephemeris files
tdrss.type                     = 'SPEphem';
tdrss.sat(1).filename          = 'TD0_SPephem_2008270.txt'; %TDRS-East
tdrss.sat(2).filename          = 'TD6_SPephem_2008270.txt'; %TDRS-West
tdrss.sat(3).filename          = 'TD4_SPephem_2008270.txt'; %TDRS-Spare
 
%% Input TDRS Antenna Data
tdrss.sat(1).antenna.maxrange  = 70000; %km
tdrss.sat(1).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(2).antenna.maxrange  = 70000; %km
tdrss.sat(2).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(3).antenna.maxrange  = 70000; %km
tdrss.sat(3).antenna.halfangle = 13*pi/180; %radians

%% Set propagator information
epoch  = datenum([2008  9 26  0  0   .000]);
tspan  = 0:10:360;
x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec
options = [];

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