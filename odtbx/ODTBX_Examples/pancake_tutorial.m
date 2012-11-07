%% OD Toolbox Tutorial 1: Planet Pancake
%
% This tutorial will teach you how to get started with using basic OD
% Toolbox functionality. In particular, you will learn how to perform the
% following analysis with ODTBX:
%
% # Numerically integrate a state vector, which includes a satellite's
% position and velocity.
% # Determine the State Transition Matrix $\Phi(t_f, t_0)$ that maps state
% changes at the initial time $t_0$ to state changes at the final time
% $t_f$.
% # Compute the observability of any element(s) of the state vector.
% # Estimate all components of the state vector.
%

%% Introduction
%
% Pancake is a circular planet that exists in a 2D universe with no other gravitational influences.
% Its rotation rate $\omega_p$, radius $r_p$, and gravitational parameter $\mu$ are all known. In
% addition, there is a satellite tracking station on the surface of
% Pancake that tracks a satellite currently in orbit around the planet.
%
% In this tutorial, your goal is to estimate the satellite's position and
% velocity, Pancake's gravitational parameter, and the tracking station's position.

%% Initialize Variables
% Suppose that at some time $t_0$, the tracking station lies on the
% inertial X axis, and the satellite is directly above the tracking station
% at an altitude of 500km with a speed of 8.339km/s and zero flight-path angle.
%
% First, initialize all values (constants, parameters, and variables) associated with Pancake,
% the satellite tracking station, and the orbiting satellite. Note that all
% vectors are given in inertial coordinates.
%

w_p = 2*pi/86400;    % [rad/s] Pancake rotation rate
mu = 3.968e5;        % [km^3/s^2] Pancake gravitational parameter
r_p = 6378;          % [km] Pancake radius
Rs0 = [r_p; 0];      % [km] Tracking station initial position
R0 = [r_p + 500; 0]; % [km] Satellite initial position
V0 = [0; 8.339];      % [km/s] Satellite initial velocity

%%
% Now define the initial value of state vector, which consists of each
% parameter that we want to estimate.

x0 = [R0; V0; mu; Rs0];

%% Numerical Integration
% Before any estimation can be performed, we need the "truth" model for the
% estimation state vector along with the state transition matrix $\phi(t_f,
% t_0)$. To create the truth model, let's use a time span of 1 day, and
% assume that the tracking station will track the satellite every 60
% seconds.

tspan = 0:60:86400; % Integration time span and step size

%%
% ODTBX provides a numerical integrator called "integ", which automatically 
% computes the state transition matrix on request. Like other Matlab
% integrators, integ requires a user-supplied dynamics function.

opts = odeset('reltol', 1e-9, 'abstol', 1e-9); % Numerical integration tolerances

% Numerically integrate the spacecraft orbit
[t x Phi] = integ(@pancake_dyn, tspan, x0, opts, w_p); % w_p will be passed to pancake_dyn()

%%
% where the function pancake_dyn(), defined in pancake_dyn.m, defines
% dynamics for both the system and the state transition matrix.

return;

%% Check Observability
%
% Use ODTBX functions to determine observability of spacecraft.
%

% Generate measurements
y = pancake_dat(t,x,[]);

% Define initial covariance
%P0 = diag([inf inf inf inf inf inf inf]);
P0 = diag([inf inf inf inf 1e-12 1e-6 1e-6]);

% Calculate observability grammian
M = observ(@pancake_dyn,@pancake_dat,tspan,x0,[],P0,y,w_p,[]);

% Check rank
rank_M = rank(M);

% Determine if observable
if rank_M < length(x0)
    fprintf('System is not observable! The observability grammian has rank %i\n', rank_M)
else
    fprintf('System is observable!\n')
end

%
%% Batch Least Squares Estimation
%
% Use the ODTBX estbat function to compute estimated states.
%

% Run estbat
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt] = ...
    estbat(@pancake_dyn,@pancake_dat,tspan,x0,P0,[],w_p,[]);


%% Output Results
%
% Two ODTBX functions (estval and plot_ominusc) and one generic Matlab
% plotting function are used to show the calculated data.
%

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