echo on
% This file is an example for the ODTBX Sequential Estimator, estseq.m
% It uses the restricted 2-body problem 'r2bp' for the dynamics model, and 
% the range and range rate from an Earth station 'rrdot3D_1way' for the 
% measurement model.
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
dynfun.tru = @r2bp; 
dynfun.est = @r2bp;
%
% Specify the inputs to the dynamics models
mu = 3.986004415e5; % km^3'sec^2
dynarg.tru = mu;
dynarg.est = mu;
%
% Define the initial reference state and a priori covariance matrix:
xo = [6878;0.00;0.00;0.00;0.00;8.339]; %km & km/s
Po = diag([1e-2 1e-4 1e-1 1e-4 1e-7 1e-5].^2); % km^2 & km^2/s^2
%
% Set up the solve-for and consider mapping
S = eye(6);                       % Solve-for map - solve for all 6 states
C=[];                         % Consider map - no consider states
% 
% The next four lines say to use the same a priori for truth and estimator:
Xnot.Xo = xo;                 % True initial state
Xnot.Xbaro = S*xo;            % Estimator initial state
Pnot.Po = Po;                 % True initial covariance
if isempty(Po)
    Pnot.Pbaro = [];
else
    Pnot.Pbaro = S*Po*S';         % Estimator initial covariance
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
% Specify the time vector at which the measurements are taken. Here, the 
% filter will run over a 5 minute pass at 10 second intervals:
tspan = 0:10:300;
%
% Set some options for the sequential filter.  We will run 30 Monte Carlo cases 
% and for each case, iterate on the sequential solution 3 times.
%
opts = setOdtbxOptions('MonteCarloCases',10,'UpdateIterations',2);
%
% Now run (and time) the sequential estimator:
%
tic
% % Run estimator without controls
% myest = estnew(dynfun,datfun,tspan,Xnot,Pnot,opts,dynarg,datarg);
% 
% [t,xhat,P,e,Y]=myest.run_estimator();

% Run estimator with controls
mysim = est_control('estnew',dynfun,datfun,tspan,Xnot,Pnot,opts,dynarg,datarg);
mysim.set_controllers();
[t,xhat,P,e,Y] = mysim.run_sim();

toc
%
% Plot estimator error, measurment innovations, and variance sandpiles in
% terms of meters
% plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,1e3);
%
echo off
