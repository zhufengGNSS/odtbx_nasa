%% OD Toolbox Tutorial 2: Batch Estimation with Multiple Ground Stations
%
% This tutorial will teach you how to use OD Toolbox functions to simulate
% measurements. In particular, you will learn how to perform the following 
% analysis with ODTBX:
%
% # Perform coordinate transformations
% # Use ODTBX dynamics functions
% # Use ODTBX measurement-generation functions
% # Manipulate the ODTBX Options structure to specify options for the ODTBX
% estimators and measurement models
% # Use ODTBX plotting routines
%
% This tutorial expands on the following tutorials: <matlab:doc('pancake_tutorial') Tutorial 1>
% See also: pancake_tutorial

%% Introduction
%
% There are 3 ground stations on the Earth, each of which are tracking the
% same satellite that is in a specified two-body orbit. In this tutorial, your goal is
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

% Create an Options structure for ground station measurements with default settings
datarg = odtbxOptions('measurement');

%%
% Here, the default settings for ground station measurements specify that ground 
% stations should make 2-way range and rangerate measurements. In addition
% to these, we need to specify our custom requirements.

%%
% Custom requirement: Measurement Epoch. This one's easy.

datarg = setOdtbxOptions(datarg, 'epoch', epoch);

%%
% Custom requirement: Ground station IDs. Assume that the three ground
% stations have IDs 'HBKS', 'USHS', and 'USPS'.

datarg = setOdtbxOptions(datarg, 'gsID', {'HBKS','USHS','USPS'});

%%
% Custom requirement: Measurement Covariance. NEEDS EXPLANATION

datarg = setOdtbxOptions(datarg, 'rSigma',...
    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);

%%
% Note that if we wanted more (or different) ground stations, we would only
% need to modify the 'gsID' and 'rSigma' options accordingly.
%
% Next we need to specify simulation models for the ground stations themselves.
% ODTBX provides a wrapper function, called 'createGroundStationList', that
% uses JAT to return a list of all stations in the NDOSL.

datarg = setOdtbxOptions(datarg, 'gsList', createGroundStationList());

%%
% Finally, we need to define the function that will create the simulated
% ground station measurements. ODTBX provides a function for this, called
% 'gsmeas'.

datfun = @gsmeas; % Reference to ground station measurement function

%%
% Now that we've specified options for the measurement functions, let's
% specify options for the estimator itself. In this case, suppose that we 
% don't want to use process noise, and we only want to iterate the solution 2 times.

% Create an Options structure for the estimator, with default settings
estarg = odtbxOptions('estimator');

estarg = setOdtbxOptions(estarg, 'UseProcNoise', false); % No process noise
estarg = setOdtbxOptions(estarg, 'UpdateIterations', 2); % Update solution twice

%% Batch Least Squares Estimation
% Now run the batch estimator, assuming we have a priori state covariance
% knowledge.

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