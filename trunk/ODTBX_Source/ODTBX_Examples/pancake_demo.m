%% Pancake Demo
%  This demo shows how to compute observability and perform batch least
%  squares estimation with ODTBX. A planar two-body orbit is propagated 
%  and estimated (assuming one ground station) with estbat.
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
%  Original credit for this example goes to K. Getzandanner & R. Carpenter 
%  (24 JUL 2012)
%  
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Phillip Anderson        08/30/2012      Created from original pancake
%   & Ravi Mathur                           demo files

echo on
%% ************************************************************************
%  This demo shows how to compute observability and perform batch least
%  squares estimation with ODTBX. A planar two-body orbit is propagated 
%  and estimated (assuming one ground station) with estbat.
% ************************************************************************

%
%% Define Initial State & Covariance
%
% Define all scenario-specific variables and ODTBX options. Specify the
% spacecraft initial state as well as the ground station location.
%
pause % Hit RET to continue

wP = 2*pi/86400;        % Pancake rotation rate
mu = 3.986e5;           % Pancake gravitational parameter
R0 = [6378 + 500; 0];   % S/C initial position
V0 = [0; 8.339];        % S/C initial velocity
Rs0 = [6378; 0];        % G/S initial state

% Define state vector
x0 = [R0; V0; mu; Rs0];

% Define initial covariance
%P0 = diag([inf inf inf inf inf inf inf]);
P0 = diag([inf inf inf inf 1e-12 1e-6 1e-6]);

%
%% Numerical Integration Using Integ
%
% Integ function takes scenario-specific variables, time span, and
% integration tolerance.
%
pause % Hit RET to continue

% Define integration time span and step size.
tspan = 0:60:86400;

% Set numerical integration tolerances
opts = odeset('reltol',1e-9,'abstol',1e-9);

% Numerically integrate the spacecraft orbit
[t x Phi] = integ(@pancake_dyn,tspan,x0,opts,wP);

%
%% Check Observability
%
% Use ODTBX functions to determine observability of spacecraft.
%
pause % Hit RET to continue

% Generate measurements
y = pancake_dat(t,x,[]);

% Calculate observability grammian
M = observ(@pancake_dyn,@pancake_dat,tspan,x0,[],P0,y,wP,[]);

% Check rank
rank_M = rank(M);
fprintf('Rank of Normal Matrix: %i\n',rank_M)

% Determine if observable
if rank_M < 7
    fprintf('System is not observable!\n')
end

%
%% Batch Least Squares Estimation
%
% Use the ODTBX estbat function to compute estimated states.
%
pause % Hit RET to continue

% Run estbat
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt] = ...
    estbat(@pancake_dyn,@pancake_dat,tspan,x0,P0,[],wP,[]);


%% Output Results
%
% Two ODTBX functions (estval and plot_ominusc) and one generic Matlab
% plotting function are used to show the calculated data.
%
pause % Hit RET to continue

% Plot MC errors and covariance
 estval(t,e,P)

% Plot residuals and Pdy
 plot_ominusc(t,dy,Pdy,Pdyt)

% Plot spacecraft and ground station position
figure;
plot(x(1,:),x(2,:),'b',x(6,:),x(7,:),'g')
title('Pancake Demo')
legend('S/C','Pancake')
xlabel('X [km]')
ylabel('Y [km]')
axis equal

% Output final state and STM
X_f = x(:,end);
Phi_f = Phi(:,:,end);

fprintf('X_f = \n')
disp(X_f)
fprintf('Phi_f = \n')
disp(Phi_f)

echo off

% These hard file references exist because plot_ominusc currently has a bug
% related to creating figures that crashes the 'publish' command. The
% references count on the 5 figures produced by this example being saved
% with the appropriate names. This section can be removed once plot_omniusc
% is fixed.
% 
% <<pancake1.png>>
%
% <<pancake2.png>>
%
% <<pancake3.png>>
%
% <<pancake4.png>>
%
% <<pancake5.png>>
%
