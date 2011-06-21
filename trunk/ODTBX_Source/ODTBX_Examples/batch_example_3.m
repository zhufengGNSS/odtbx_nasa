function batch_example_3
% Example for the ODTBX Batch Estimator, estbat.m, similar to batch_example_1 with no process noise.
%
% This file is an example for the ODTBX Batch Estimator, estbat.m
% It uses the JAT Force Model of the 2-body gravity for the dynamics model
% but modified with a wrapper function (defined in this file) to set the
% process noise equal to zero.  The measurement model is rrdot3D_1way for 
% range and range rate.
%      This example is similar to batch_example_1 and batch_example_1 except 
% that the process noise is set to zero and the time span was extended to 
% 60,000 sec.
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

echo on
%
% Specify the dynamics models for the truth and the estimator.  Both models
% are identical.
dynfun.tru = @jatForces_km_noQ; % See the subfunction at the end
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
% Specify the measurement models for the truth and the estimator.  Here, we 
% are going to use the same models for both.
datfun.tru = @rrdot3D_1way;
datfun.est = @rrdot3D_1way;
%
% Specify the inputs to the measurement models, which in this case is the 
% measurement noise.  Here, the estimator will assume 3x more measurement 
% noise than the truth:
sig = diag([1e-3,1e-6]); % km & km/s
datarg.tru = sig;
datarg.est = 3*sig;
%
% Specify the time vector at which the measurements are taken. The time
% periods when measurements are available were determined, and tspan is
% specified to cover just those time periods.
dT      = 10;
tspan   = [0:dT:60000];
%
% Set the estimator options, for more options see Estbat.m:
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'UpdateIterations',2,'MonteCarloCases',10,...
    'OdeSolver',@ode113,'OdeSolvOpts',...
    odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10));
%
% Run and (time) the batch filter:
tic
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt]=estbat(dynfun,datfun,...
    tspan,x0,P0,eOpts,dynarg,datarg);  
toc
%
% Plot estimator error, measurment residuals, and variance sandpiles in
% terms of meters
plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,1e3);
%
echo off

%
% ------------------------------------------------------------------------
%
% Dynamics model used in the example above.
function [xDot,A,Q] = jatForces_km_noQ(t,x,jatWorld)
[xDot,A,Q] = jatForces_km(t,x,jatWorld);
Q=zeros(size(Q));

