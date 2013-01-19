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

%% Introduction
%
% Pancake is a circular planet that exists in a 2D universe with no other gravitational influences.
% Its rotation rate $\omega_p$, radius $r_p$, and gravitational parameter $\mu$ are all known. In
% addition, there is a satellite ground station on the surface of
% Pancake that tracks a satellite currently in orbit around the planet.
%
% In this tutorial, your goal is to estimate the satellite's position and
% velocity, Pancake's gravitational parameter, and the ground station's position.

%% Initialize Variables
% Suppose that at some time $t_0$, the ground station lies on the
% inertial X axis, and the satellite is directly above the ground station
% at an altitude of 500km with a speed of 8.339km/s and zero flight-path angle.
%
% First, initialize all values (constants, parameters, and variables) associated with Pancake,
% the satellite ground station, and the orbiting satellite. Note that all
% vectors are given in inertial coordinates.

w_p = 2*pi/86400;    % [rad/s] Pancake rotation rate
mu = 3.986e5;        % [km^3/s^2] Pancake gravitational parameter
r_p = 6378;          % [km] Pancake radius
Rs0 = [r_p; 0];      % [km] ground station initial position
R0 = [r_p + 500; 0]; % [km] Satellite initial position
V0 = [0; 8.339];     % [km/s] Satellite initial velocity

%%
% Now define the initial value of state vector, which consists of each
% parameter that we want to estimate.

x0 = [R0; V0; mu; Rs0];

%% Numerical Integration
% Before any estimation can be performed, we need the "truth" model for the
% estimation state vector along with the state transition matrix $\phi(t_f,
% t_0)$. To create the truth model, let's use a time span of 1 day, and
% assume that the ground station will track the satellite every 60
% seconds.

tspan = 0:60:86400; % Integration time span and step size

%%
% ODTBX provides a numerical integrator called "integ", which automatically 
% computes the state transition matrix on request. Like other Matlab
% integrators, integ requires a user-supplied dynamics function.

opts = odeset('reltol', 1e-9, 'abstol', 1e-9); % Numerical integration tolerances

% Numerically integrate the spacecraft orbit and state transition matrix
[t x Phi] = integ(@pancake_dyn, tspan, x0, opts, w_p); % w_p will be passed to pancake_dyn()

%%
% where the function pancake_dyn(), defined in pancake_dyn.m, defines
% dynamics for both the system and the state transition matrix.

%% Check Observability
% A state's observability indicates how accurately that state can be estimated
% using only measured quantities. In this case, suppose the ground station measures
% the range to the satellite when the satellite is in view (i.e. above the
% local horizon).

y = pancake_dat(t,x,[]); % Generate range measurements

%%
% where the function pancake_dat(), defined in pancake_dat.m, computes the
% range (from the station to the satellite) for each time, as long as the
% satellite is visible at that time. If we have no prior knowlege about how
% the measured range is related to the estimated states, we can indicate
% that via the a priori covariance matrix P0,

P0 = diag([inf inf inf inf inf inf inf]); % Unknown initial covariance

%%
% ODTBX provides a function called "observ", which computes the
% observability gramian of a system given observations that have been made.
% The a priori covariance information can be provided to observ in order to
% improve the observability estimate.

M = observ(@pancake_dyn,@pancake_dat,tspan,x0,[],P0,y,w_p,[]); % Compute observability gramian

%%
% The observability gramian M is 7x7, corresponding to the number of
% estimated states. These states are observable if (and only if) M is full
% rank, so let's check the observability.

rank_M = rank(M); % Get rank of observability gramian

% The system is observable if M is full rank
if rank_M < length(M)
    fprintf('System is not observable! The observability gramian has rank %i\n', rank_M)
else
    fprintf('System is observable!\n')
end

%%
% In this case, the system is unobservable because knowledge of range alone
% cannot be used to differentiate between the station and spacecraft
% positions. However, if we have prior knowledge of how well the range 
% measurement is related to the states, then we can provide this via the a 
% priori covariance matrix P0. Suppose that we know how well the measured
% range relates to the gravitational parameter and station position,

P0 = diag([inf inf inf inf 1e-12 1e-6 1e-6]); % Define initial covariance

%%
% Note that we don't have to provide ALL covariance estimates; we can pick
% and choose depending on our knowlege. Let's re-run the observability
% function and see what the gramian indicates

M = observ(@pancake_dyn,@pancake_dat,tspan,x0,[],P0,y,w_p,[]); % Compute observability gramian
rank_M = rank(M); % Get rank of observability gramian

% The system is observable if M is full rank
if rank_M < length(M)
    fprintf('System is not observable! The observability gramian has rank %i\n', rank_M)
else
    fprintf('System is observable!\n')
end

%% Batch Least Squares Estimation
% Assuming that all states are observable, let's estimate the states at
% all times using a batch least squares estimation method. To do this,
% ODTBX provides the function "estbat".

myest = estbat;

[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt] = ...
    myest.run_estimator(@pancake_dyn,@pancake_dat,tspan,x0,P0,[],w_p,[]);

% [t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt] = ...
%     estbat(@pancake_dyn,@pancake_dat,tspan,x0,P0,[],w_p,[]);

%% Output Results
% It is important to know the errors incurred by an estimation process. In
% particular, we are interested in how well the estimator errors lie within
% their predicted error bounds, and how closely the observed measurements
% match the computed measurements.

%%
% ODTBX provides the function "estval" to plot errors in estimated statesm, relative to
% their predicted error bounds. Each element of the estimated state is
% plotted separately, with the 2-$\sigma$ error bounds overlaid onto the
% plot.

estval(t,e,P);

%%
% ODTBX provides the function "plot_ominusc" to plot errors in the
% measurement, relative to its predicted error bound. 

plot_ominusc(t,dy,Pdy,Pdyt);

%%
% Finally, for good measure, let's plot the orbit of the spacecraft around
% Pancake. We'll use the ground station position to represent the planet's
% surface.

figure;
plot(x(1,:),x(2,:),'b',x(6,:),x(7,:),'g')
title('Pancake Demo')
legend('S/C','Pancake')
xlabel('X [km]')
ylabel('Y [km]')
axis equal