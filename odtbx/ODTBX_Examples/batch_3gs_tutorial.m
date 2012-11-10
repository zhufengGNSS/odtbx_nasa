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
% all tracking stations.

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
% Suppose the tracking stations take measurements every 3 minutes for 12
% hours. The corresponding simulation times are,

tspan = 0:180:(12*3600); % Every 3 minutes for 12 hours


%%
% Let's finish up general initialization by defining the 2-body dynamics
% model specified in the introduction. ODTBX provides a dynamics function
% called 'r2bp' that implements the appropriate dynamics, including any
% necessary State-Transition Matrix and Process Noise models.

dynfun = @r2bp;    % 2-body dynamics with an input gravitational paramter
dynarg = mu_Earth; % Input gravitational parameter for r2bp

%% The OD Toolbox Options Structure
% 

% Define the measurements models for the truth and the estimator.  Here 
% we use the same model for both.
datfun = @gsmeas;
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
%datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS','WHSX'});
%datarg = setOdtbxOptions(datarg,'rSigma',...
%    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS'});
datarg = setOdtbxOptions(datarg,'rSigma',...
    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsList',createGroundStationList());

% Specify that we don't want to use process noise, and only want to iterate 
% the solution 2 times.
estarg = odtbxOptions('estimator');
estarg = setOdtbxOptions(estarg,'UseProcNoise',false);
estarg = setOdtbxOptions(estarg,'UpdateIterations',1);

%% 
% Now run the batch estimator:

P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);

[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt]= ...
    estbat(dynfun,datfun,tspan,x0,P0,estarg,dynarg,datarg);

%% Plot Results
% Plot estimator error, measurment residuals, and variance sandpiles in
% terms of meters. Also plot the true trajectory in 3D against the Earth.

plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt);

EarthOrbitPlot(xhat{1});