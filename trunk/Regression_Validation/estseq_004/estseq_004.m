% Regression Test Case
% estseq_004: Changed EditFlag to Inhibit (NOTE: EditFlag has been removed
% as an option)
%    Models:
%        Dynamics: r2bp.m
%        Measurement: rrdot3D.m
%    Initial Conditions:
%        x0 = [6878;0;0;8.339]; km & km/sec
%        P0 = diag([1e-2 1e-2 1e-2 1e-5 1e-5 1e-5].^2); km2 & km2/sec
%        tspan = [0:10:600]; s
%        mu = 3.986e5; km3/sec2
%        sig = diag([1e-3 1e-6]); km & km/sec
%    Options:
%        integrator: ode45
%            Default Options
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

function [failed] = estseq_004
cd estseq_004
testCaseName = 'estseq_004';
format long

mu = 3.986e5;
sig = diag([1e-3 1e-6]);
tspan = [0:10:600];
x0 = [6878;0;0;0;0;8.339];
P0 = diag([1e-2 1e-2 1e-2 1e-5 1e-5 1e-5].^2);

opts = setOdtbxOptions('UpdateIterations',1,'MonteCarloCases',1,...
    'ValidationCase', 1,...
    'OdeSolvOpts',odeset('Vectorized','on','InitialStep',300));
        
[t,x,P,e,dy,Pdy,arf] = estseq(@r2bp,@rrdot3D,tspan,x0,P0,opts,mu,sig);

% Need to set where the measurements end so we don't have NaNs for Traj &
% Covariance
theEnd = min(find(isnan(P{1}(1,:))))-1;

% This bit of code is used to create the baseline covariance file, it
% should only be used the first time to create a baseline for covariance,
% measurements and state and all code below it should be commented out when it is run.
%  ODTBX_Covariance_Time_Zero = (unscrunch(P{1}(:,1)))*10^6;
%  ODTBX_Covariance_Time_Last = (unscrunch(P{1}(:,theEnd)))*10^6;
%  save covariance_baseline_estseq_004.mat ODTBX_Covariance_Time_Zero ODTBX_Covariance_Time_Last
%  failed = 0;

%Output the trajectory result to a file
trajOutput(t{1}(:,1:theEnd),0,x{1}(:,1:theEnd));



%Call comparison scripts to compare trajectory & measurement
%differences between current run & baseline run
plotting = 0; % Do Not Plot
regressionTest = 1; % Running a regression test
[agreement, measdiff, trajdiff]= Compare_Old_New(plotting,regressionTest,testCaseName); %Will compare current run with baseline run for measurement & trajectory comparison
[covAgreement] = compareCovariance(P{1},testCaseName); % Will compare covariance against baseline and write a file with different

if ((agreement == 1) || (covAgreement == 1))
   failed = 1;
else
   failed = 0;
end

cd ..   % return to regression directory