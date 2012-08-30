%% Batch (3 Ground Stations) Demo
% This demo shows how to compute observability using the ODTBX function
% estbat. Three ground stations are used to observe the spacecraft, which
% is assumed to operate in a two-body dynamics system.
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
%
%  Original credit for this example goes to K. Getzandanner.
%  
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Phillip Anderson        08/30/2012      Created from original demo
%   & Ravi Mathur                           file demo1.m

echo on
%
% ************************************************************************
% This demo shows how to compute observability using the ODTBX function
% estbat. Three ground stations are used to observe the spacecraft, which
% is assumed to operate in a two-body dynamics system.
% ************************************************************************

%
%% Define Initial State & Covariance
% Define all scenario-specific variables and ODTBX options.
%
pause % Hit RET to continue

epoch = datenum('01 April 2012 04:27:46');
mu_Earth   = 398600.4415; %km^3/s^2

%
% Spacecraft state using Keplerian elements
KOE.sma = 702+6378.1;
KOE.ecc = 0.0001232;
KOE.incl = 98.2280*pi/180;
KOE.raan = 	33.8594*pi/180;
KOE.argp = 90.8017*pi/180;
KOE.tran = 0;

%
% Define the initial reference state and a priori covariance matrix
x0 = kep2cart(KOE,mu_Earth);
P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);

%
% Define the times (epoch seconds) to process measurements and estimate
tspan = 0:180:(12*3600); %Every 3 minutes for 12 hours

%
%% Define estimator-specific options
%
pause % Hit RET to continue

%
% Define the dynamics model. Here we use the same model for the truth 
% and the estimator.
pause % Hit RET to continue

dynfun = @r2bp;
dynarg = mu_Earth;

%
% Define the measurements models for the truth and the estimator.  Here 
% we use the same model for both.
datfun = @gsmeas;
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS','WHSX'});
datarg = setOdtbxOptions(datarg,'rSigma',...
    [1e-2 1e-5 1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsList',createGroundStationList());

%
% Specify that we don't want to use process noise, and only want to iterate 
% the solution 2 times.
estarg = odtbxOptions('estimator');
estarg = setOdtbxOptions(estarg,'UseProcNoise',false);
estarg = setOdtbxOptions(estarg,'UpdateIterations',1);

%
%% Now run (and time) the batch estimator:
%
pause % Hit RET to continue

tic
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt]= ...
    estbat(dynfun,datfun,tspan,x0,P0,estarg,dynarg,datarg);
toc

%
%% Plot Results
% Plot estimator error, measurment residuals, and variance sandpiles in
% terms of meters. Also plot the true trajectory in 3D against the Earth.
pause % Hit RET to continue

EarthOrbitPlot(xhat{1});

pause % Hit RET to continue

plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt)

echo off