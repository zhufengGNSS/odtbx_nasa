%% OD Toolbox Tutorial 2: Batch Estimation with Multiple Ground Stations
% This tutorial will teach you how to use OD Toolbox functions to simulate
% measurements. In particular, you will learn how to perform the following 
% analysis with ODTBX:
%
% # Transform Keplerian elements to Cartesian inertial states
% # Use ODTBX dynamics functions
% # Use ODTBX measurement-generation functions
% # Manipulate the ODTBX Options structure to specify options for the ODTBX
% estimators and measurement models
% # Use ODTBX plotting routines
%
% This tutorial expands on the following tutorials: <matlab:doc('pancake_tutorial') Pancake>
%
% See also: pancake_tutorial
%
% You can either run this example after setting "echo on" to see output in
% the command line, or do "publish chaser_target_tutorial" to create formatted
% html output in the ODTBX_Examples/html directory.

%% Introduction
% There are 3 ground stations on the Earth, each of which are tracking a 
% satellite that is in a specified two-body orbit. In this tutorial, your goal is
% to estimate the satellite's position and velocity using measurements from
% all ground stations.

%% Initialize Variables
% Suppose that on April 1 2012 at 4:27:46am, the satellite is at the periapse
% of a specified orbit. The date must be specified as a date number, and
% let's choose to specify the orbit by its Keplerian orbital elements.

epoch = datenum('01 April 2012 04:27:46');

mu_Earth   = 398600.4415; % [km^3/s^2] Earth gravitational parameter

% Spacecraft state using Keplerian elements
KOE.sma = 702+6378.1;       % [km] Semimajor Axis
KOE.ecc = 0.0001232;        % [] Eccentricity
KOE.incl = 98.2280*pi/180;  % [rad] Inclination
KOE.raan = 	33.8594*pi/180; % [rad] Ascending Node
KOE.argp = 90.8017*pi/180;  % [rad] Argument of Perigee
KOE.tran = 0;               % [rad] True Anomaly

%% 
% ODTBX requires inertial vectors, and provides the function "kep2cart" to
% create them from Keplerian elements. Since we only want to estimate the
% satellite's position and velocity, our estimation state vector is the
% same as the satellite's state vector.

x0 = kep2cart(KOE, mu_Earth); % Keplerian -> Inertial

%%
% Suppose the ground stations take measurements every 3 minutes for 12
% hours. The corresponding simulation times are,

tspan = 0:180:(12*3600); % Every 3 minutes for 12 hours

%%
% Let's finish up general initialization by defining the 2-body dynamics
% model specified in the introduction. ODTBX provides a dynamics function
% called 'r2bp' that implements the appropriate dynamics, including any
% necessary State-Transition Matrix and Process Noise models.

dynfun = @r2bp;    % Reference to 2-body dynamics function
dynarg = mu_Earth; % Input gravitational parameter for r2bp

%% The OD Toolbox Options Structure
% In order to simulate measurements from ground stations, we need some way
% to specify several properties of the stations themselves. ODTBX provides
% an "Options" structure to organize and facilitate the specification of
% a set of ground stations. While this structure offers many options, they
% all have default values.

% Create an Options structure for ground station measurements with default values
datarg = odtbxOptions('measurement')

%%
% Note that all measurement options seem to be empty! This actually
% indicates that defaults will be used for all empty fields. The default 
% options for ground station measurements specify that ground 
% stations should make 2-way range and rangerate measurements. To see all
% possibilities for measurement options, simply omit the output argument;
% defaults will be indicated in braces.

odtbxOptions('measurement'); % List all measurement options, with defaults in braces

%%
% In addition to these default options, we need to specify our custom requirements.
% First, let's specify the epoch that should be used,

datarg = setOdtbxOptions(datarg, 'epoch', epoch);

%%
% Now we need to specify simulation models for the ground stations.
% ODTBX provides a wrapper function, called 'createGroundStationList', that
% uses JAT to return a list of all stations in the NASA Directory of
% Station Locations (NDOSL).

datarg = setOdtbxOptions(datarg, 'gsList', createGroundStationList());

%%
% Next, let's specify the NDOSL ID for each ground station that we want to
% use. Assume that the three ground stations have IDs 'HBKS', 'USHS', and
% 'USPS'.

datarg = setOdtbxOptions(datarg, 'gsID', {'HBKS','USHS','USPS'});

%%
% Finally, if we know the uncertainties in measurements from each station,
% we can specify them here. The 'rSigma' option is a vector indicating the
% $\sigma$ uncertainty of each selected measurement type, for each station.

datarg = setOdtbxOptions(datarg, 'rSigma',...
    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]); % [km, km/s, km, km/s, km, km/s] sigma uncertainties 

%%
% Here, both measurement sigmas are given for the first ground station, then
% repeated for each additional ground station in the 'gsID' list. Note that 
% if we wanted to use more (or different) ground stations, we would only need to
% modify the 'gsID' and 'rSigma' options accordingly.
%
% Before continuing, let's check that all of our custom measurement options
% have been set correctly.

datarg

%%
% Finally, we need to define the function that will create the simulated
% ground station measurements. ODTBX provides a function for this, called
% 'gsmeas'.

datfun = @gsmeas; % Reference to ground station measurement function

%%
% Now that we've specified options for the measurement function, let's
% specify options for the estimator itself. Conveniently, ODTBX uses the
% same "Options" structure to organize estimator options, and the same
% access functions to create and set these options.

% Create an Options structure for the estimator with default values
estarg = odtbxOptions('estimator');

%%
% As before, the default estimator options can be listed by omitting the
% output argument.

odtbxOptions('estimator'); % List all estimator options, with defaults in braces

%%
% Let's specify that we don't want to include process noise, and that we
% want to update the estimated solution twice.

estarg = setOdtbxOptions(estarg, 'UseProcNoise', false); % No process noise
estarg = setOdtbxOptions(estarg, 'UpdateIterations', 2); % Update solution twice

%% Batch Least Squares Estimation
% Now run the batch estimator, assuming we have a priori state covariance
% knowledge. This step is similar to the estimation from Tutorial 1.

P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2); % Define initial covariance

% Start the batch estimator with specified functions and options
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt]= ...
    estbat(dynfun,datfun,tspan,x0,P0,estarg,dynarg,datarg);

%% Plot Results
% To visualize the estimation results, let's plot estimator error,
% measurement residuals, and variance sandpiles. ODTBX provides the
% 'plot_results' function, which conveniently plots each of these
% quantities for all ground stations and estimated states.

plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt);

%%
% Finally, for good measure, let's plot the estimated orbit of the spacecraft. 
% ODTBX provides the 'EarthOrbitPlot' function to easily plot 3D
% orbits around a texture-mapped Earth.

EarthOrbitPlot(xhat{1});

%% License 
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