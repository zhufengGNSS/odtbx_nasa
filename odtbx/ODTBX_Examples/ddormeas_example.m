%DDORMEAS_EXAMPLE This demonstrates the use of the ddor measurement model.
% The increasing covariance shows that ddor should be combined with another
% measurement type to get a good solution
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

%% initialize
tic; %starts the matlab timer 
format long g; %forces Command Window outputs to be more precise


%% Set the user satellite state
epoch    = datenum('01 Jul 2007 14:00:00.000');
x0 = [1.5e8 0 0 0 5e-2 0]';% km & km/sec
P0 = diag([10 10 10 .1 .1 .1].^2); %initial covariance matrix (km^2 & km^2/sec^2)

%% Set the Tracking Stations
gsID = {'DS15',... Goldstone 34 Meter
        'DS45',... Canberra 34 Meter
        'DS65',... Madrid 34 Meter
        }; %The NASA ID for each ground station
gsSigma = [1e-6 ... (km range sigma) Goldstone 34 Meter
           1e-6 ... (km range sigma) Canberra 34 Meter
           1e-6 ... (km range sigma) Madrid 34 Meter
           ]; %The measurement sigmas for each ground station

%% Set a Tracking Schedule
% Stations(1:3) = (Goldstone, Canberra, Madrid)
% { gs_number1 gs_number2   start_time   stop_time } times in datenum recognized string
TrackSched = {
    1 2 '01 Jul 2007 14:00:00' '01 Jul 2007 18:00:00'
    1 3 '02 Jul 2007 07:30:00' '02 Jul 2007 11:00:00'
    2 1 '02 Jul 2007 14:00:00' '02 Jul 2007 18:00:00'
    1 3 '03 Jul 2007 07:30:00' '03 Jul 2007 11:00:00'
    1 2 '03 Jul 2007 14:00:00' '03 Jul 2007 18:00:00'
    3 1 '04 Jul 2007 07:30:00' '04 Jul 2007 11:00:00'
    };
dT = 600; % time step of the measurement

%% Re-format the tracking data schedule for seconds from epoch
ddor.Sched      = cell2mat(TrackSched(:,1:2)); %1st and 2nd column = ground station numbers
ddor.Sched(:,3) = (datenum(TrackSched(:,3))-epoch)*86400; %3nd column = start times
ddor.Sched(:,4) = (datenum(TrackSched(:,4))-epoch)*86400; %4rd column = end times

%% Estimator options
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0,...
    'OdeSolver',@ode45,'OdeSolvOpts', odeset('reltol',1e-5,'abstol',...
    1e-5,'initialstep',10,'Vectorized','on'));

%% Set up the options and java object for JAT forces
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', '2Body');
jOptions = setOdtbxOptions(jOptions, 'useSolarGravity', true);
jOptions = setOdtbxOptions(jOptions, 'useLunarGravity', true);
jOptions = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', false);
jOptions = setOdtbxOptions(jOptions, 'useAtmosphericDrag', false);
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
% relevant information for propagating Earth-centric orbits using JAT
% force models. 

%% Set the measurement options
gsList  = createGroundStationList; %generates a list of all of the ground stations

% Initialize options structure. In this case it is a measurement structure
measOptions = odtbxOptions('measurement');

% Choose from any or all of these options and enter desired values
measOptions = setOdtbxOptions(measOptions,'gsID',gsID);
measOptions = setOdtbxOptions(measOptions,'gsList',gsList);
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %datenum format
measOptions = setOdtbxOptions(measOptions,'rSigma',gsSigma);
measOptions = setOdtbxOptions(measOptions,'ddor',ddor);

%% Setup the estimator
dynfun = @jatForces_km;
datfun = @ddormeas;
dynarg = jatWorld;
datarg = measOptions;

%% Set the time to cover all tracking intervals
tspan = [];
for n=1:size(ddor.Sched,1)
    tspan = union(tspan, ddor.Sched(n,3):dT:ddor.Sched(n,4));
end

%% Run the sequential estimator with consider covariance analysis
[t,xhat,Phat,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw]=estseq(dynfun,datfun,tspan,x0,P0,eOpts,dynarg,datarg);

%% Create a "truth" propagation
[xprop] = jatForces_km(tspan, x0, dynarg);

%% Display the initial and final covariance values in Phat
Phat0 = unscrunch(Phat{1});
disp('Initial Covariance: ')
disp(Phat0(:,:,1))
disp(' ')
disp('Final Covariance: ')
Phat0(:,:,end)
disp(' ')
disp('The increasing X component of the covariance shows that ddor should ')
disp('be combined with another measurement type to get a good solution')
