echo on
% This file is an example for the ODTBX Sequential Estimator, estseq.m
% It uses the JAT Force Model of the 2-body gravity for the dynamics model 
% and gsmeas function for the measurement model of range and range rate.
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

% Specify the dynamics models for the truth and the estimator.  Here, we 
% are going to use the same models for both.
dynfun.tru = @jatForces_km; 
dynfun.est = dynfun.tru;
%
% Set up the options and java object for JAT forces
% Initialize options structure
jOptions    = odtbxOptions('force');
%
% Choose from any or all of these options and enter desired values
jOptions    = setOdtbxOptions(jOptions, 'epoch', JATConstant('MJDJ2000') );
jOptions    = setOdtbxOptions(jOptions, 'cD', 2.2);
jOptions    = setOdtbxOptions(jOptions, 'cR', 0.7);
jOptions    = setOdtbxOptions(jOptions, 'mass', 1000);
jOptions    = setOdtbxOptions(jOptions, 'draga', 20, 'srpArea', 20);
jOptions    = setOdtbxOptions(jOptions, 'earthGravityModel', '2body');
jOptions    = setOdtbxOptions(jOptions, 'gravDeg', 2, 'gravOrder', 2);
jOptions    = setOdtbxOptions(jOptions, 'useSolarGravity', false);
jOptions    = setOdtbxOptions(jOptions, 'useLunarGravity', false);
jOptions    = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', false);
jOptions    = setOdtbxOptions(jOptions, 'useAtmosphericDrag', false);
jOptions    = setOdtbxOptions(jOptions, 'atmosphereModel', 'HP');
jOptions    = setOdtbxOptions(jOptions, 'nParameterForHPModel', 2);
jOptions    = setOdtbxOptions(jOptions, 'f107Daily', 150);
jOptions    = setOdtbxOptions(jOptions, 'f107Average', 150);
jOptions    = setOdtbxOptions(jOptions, 'ap', 15);
%
dynarg.tru = createJATWorld(jOptions); 
dynarg.est = dynarg.tru;
%
% Define the initial reference state and a priori covariance matrix:
x0 = [6878;0.00;0.00;0.00;0.00;8.339];         % km & km/sec
P0 = diag([1e-2 1e-4 1e-1 1e-4 1e-7 1e-5].^2); % km^2 & km^2/s^2
%
% Set up the solve-for and consider mapping
S = eye(6);                       % Solve-for map - solve for all 6 states
C=[];                             % Consider map - no consider states
%
% The next four lines say to use the same a priori for truth and estimator:
Xnot.Xo = x0;                     % True initial state
Xnot.Xbaro = S*x0;                % Estimator initial state
Pnot.Po = P0;                     % True initial covariance
if isempty(P0)
    Pnot.Pbaro = [];
else
    Pnot.Pbaro = S*P0*S';         % Estimator initial covariance
end
%
% Specify the measurement models for the truth and the estimator.  The
% models will be identical except the measurement noise for the estimator
% is 3 times as bad as the truth.
datfun.tru = @gsmeas;
datfun.est = datfun.tru;
%
% Set the measurement options
epoch   = datenum(2006,12,31,23,59,38.3)
%gsID    = {'ZZOD','DS16','DS46','DS66'};
gsID    = {'ZZOD'};
gsList  = createGroundStationList('DBS_NDOSL_WGS84_Mod_Example.txt');
gsECEF  = getGroundStationInfo(gsList,gsID{1},'ecefPosition',epoch);
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange',true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useRangerate',true);
measOptions = setOdtbxOptions(measOptions,'useDoppler',false);
measOptions = setOdtbxOptions(measOptions,'rSigma',[1e-3 1e-6]);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',0);
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
%
datarg.tru  = measOptions;
datarg.est  = setOdtbxOptions(datarg.tru,'rSigma',3*[1e-3 1e-6]);
%
% Specify the time vector at which the measurements are taken. Here, the 
% filter will run over a 5 minute pass at 1 second intervals:
dT      = 10;
tspan   = 0:dT:300; 
%
% Set the estimator options, for more options see estseq.m:
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'UpdateIterations',2,'MonteCarloCases',10,...
    'OdeSolver',@ode113,'OdeSolvOpts',...
    odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10));
%
% Run and (time) the sequential filter:
tic
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,eflg,Pdy,Pdyt]=estseq(...
    dynfun,datfun,tspan,Xnot,Pnot,eOpts,dynarg,datarg);
toc
%
% Plot estimator error, measurment innovations, and variance sandpiles in
% terms of meters
plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,1e3);
%
echo off


