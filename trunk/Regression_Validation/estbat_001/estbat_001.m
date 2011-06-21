% Regression Test Case
% Function estbat_001: Basic Validation Test Case
%    Models:
%        Dynamics: r2bp.m
%        Measurement: rrdot3D.m
%    Initial Conditions:
%        x0 = [6878;0.00;0.00;0.00;0.00;8.339]; km & km/sec
%        P0 = diag([1 1 1 1e-3 1e-3 1e-3].^2); km2 & km2/sec
%        tspan = [0:10:86400]; s
%        mu = 3.986004415e5; km3/sec2
%        sig = diag([1e-3,1e-6]); km & km/sec
%    Options:
%        integrator: ode45
%            odeset('Vectorized','on','MaxStep',10) sec
%        DatVectorized: off
%        DatJtolerance: N/A
%        DatJpattern: N/A
%        UpdateIterations: 1
%        MonteCarloCases: 1
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

function [failed] = estbat_001;
cd estbat_001;
testCaseName = 'estbat_001';
format long;

% The next two variables are inputs to the dynamics and measurement models:
mu = 3.986004415e5; % km^3'sec^2
sig = diag([1e-3,1e-6]); % km & km/s for rrdot3D

% Run the filter for a day with measurements every 10 seconds:
tspan = [0:10:86400]; %Measurement step size and the length of the run are set here

% The next two variables define the initial reference state and a priori
% covariance matrix:
x0 = [6878;0.00;0.00;0.00;0.00;8.339]; % For 3D measurement model (rrdot3D & r2bp) in km & km/sec
P0 = diag([1 1 1 1e-3 1e-3 1e-3].^2); % For 3D measurement model (rrdot3D & r2bp) in km^2 & km^2/sec^2

% Set some options, for more options see Estbat.m:
opts = setOdtbxOptions('UpdateIterations',1,'MonteCarloCases',1,...
    'ValidationCase', 1,...
    'OdeSolvOpts',odeset('Vectorized','on','InitialStep',10));

% Run the batch filter:
[t,xhat,P,e,dy]=estbat(@r2bp, @rrdot3D,tspan,x0,P0,opts,mu,sig);

% This bit of code is used to create the baseline covariance file, it
% should only be used the first time to create a baseline for covariance,
% measurements and state and all code below it should be commented out when it is run.
% ODTBX_Covariance_Time_Zero = (unscrunch(P(:,1)))*10^6;
% ODTBX_Covariance_Time_Last = (unscrunch(P(:,end)))*10^6;
% save covariance_baseline_estbat_001.mat ODTBX_Covariance_Time_Zero ODTBX_Covariance_Time_Last
% failed = 0;

%Call comparison scripts to compare trajectory & measurement
%differences between current run & baseline run
plotting = 0; % Do Not Plot
regressionTest = 1; % Running a regression test
[agreement, measdiff, trajdiff]= Compare_Old_New(plotting,regressionTest, testCaseName); %Will compare current run with baseline run for measurement & trajectory comparison
[covAgreement] = compareCovariance(P, testCaseName); % Will compare covariance against baseline and write a file with different

if ((agreement == 1) || (covAgreement == 1));
    failed = 1;
else
    failed = 0;
end

cd ..   % return to regression directory