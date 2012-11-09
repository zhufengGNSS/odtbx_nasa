%% OD Toolbox Tutorial 2: Batch Estimation with Multiple Ground Stations
%
% This demo shows how to compute observability using the ODTBX function
% estbat. Three ground stations are used to observe the spacecraft, which
% is assumed to operate in a two-body dynamics system.

%% Define Initial State & Covariance
% Define all scenario-specific variables and ODTBX options.

epoch = datenum('01 April 2012 04:27:46');
mu_Earth   = 398600.4415; %km^3/s^2

% Spacecraft state using Keplerian elements
KOE.sma = 702+6378.1;
KOE.ecc = 0.0001232;
KOE.incl = 98.2280*pi/180;
KOE.raan = 	33.8594*pi/180;
KOE.argp = 90.8017*pi/180;
KOE.tran = 0;

%% Define the initial reference state and a priori covariance matrix

x0 = kep2cart(KOE,mu_Earth);
P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);

% Define the times (epoch seconds) to process measurements and estimate
tspan = 0:180:(12*3600); %Every 3 minutes for 12 hours

%% Define estimator-specific options

% Define the dynamics model. Here we use the same model for the truth 
% and the estimator.
dynfun = @r2bp;
dynarg = mu_Earth;

% Define the measurements models for the truth and the estimator.  Here 
% we use the same model for both.
datfun = @gsmeas;
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS','WHSX'});
datarg = setOdtbxOptions(datarg,'rSigma',...
    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsList',createGroundStationList());

% Specify that we don't want to use process noise, and only want to iterate 
% the solution 2 times.
estarg = odtbxOptions('estimator');
estarg = setOdtbxOptions(estarg,'UseProcNoise',false);
estarg = setOdtbxOptions(estarg,'UpdateIterations',1);

%% Now run the batch estimator:

[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt]= ...
    estbat(dynfun,datfun,tspan,x0,P0,estarg,dynarg,datarg);

%% Plot Results
% Plot estimator error, measurment residuals, and variance sandpiles in
% terms of meters. Also plot the true trajectory in 3D against the Earth.

EarthOrbitPlot(xhat{1});

plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt)