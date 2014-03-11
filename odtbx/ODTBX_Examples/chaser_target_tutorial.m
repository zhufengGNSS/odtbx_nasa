%% OD Toolbox Tutorial 3: Chaser-Target Estimation with Multiple Dynamics/Measurement Models
% This tutorial will teach you how to combine multiple measurement and
% dynamics models in a single OD Toolbox estimation run. In particular, you 
% will learn how to perform the following analysis with ODTBX:
%
% # Perform coordinate transformations using ODTBX-created direction cosine
% matrices.
% # Use solve-for and consider variables to specify which parameters should
% be solved for during the estimation process.
% # Define separate functions and options for the truth and estimated
% models used by the ODTBX estimators.
% # Use a wrapper function to combine the dynamics and measurement
% functions of multiple spacecraft.
% # Use ODTBX plotting routines.

%% Introduction
% Two spacecraft, a target and a chaser, are in prescribed two-body
% Earth-centered orbits. The chaser spacecraft is equipped with GPS and a
% ranging device that allows for "pose" measurements. In this tutorial,
% your goal is to estimate the target satellite's state using pose
% measurements from the chaser satellite.


%% Initialize Variables
% Suppose that the epoch time and target vehicle states are known, and we
% want to estimate the satellite orbits over 300 seconds.

epoch = datenum('13 DEC 2010 00:00:0.000'); % Epoch time

% Target vehicle initial state
Xst = [ 5426.510869339219;  % [km] X
      -41803.696487556008;  % [km] Y
        -240.913638515624;  % [km] Z
           3.049219473653;  % [km/s] VX
           0.396422336331;  % [km/s] VY
          -0.053147501025]; % [km/s] VZ
      
tspan = 0:10:300; % 300 second estimation time interval

%%
% Assume that the chaser vehicle's initial state is known with respect to
% the target vehicle in a Radial-Intrack-Crosstrack reference frame. ODTBX
% has a function, 'DCM()', which creates a direction cosine matrix for
% certain reference frame transformations. In this case, we want to use the
% 'RIC' frame based on the target spacecraft, so let's first create it.

RIC_IJK = dcm('ric',Xst(1:3),Xst(4:6)); % Create RIC transformation from radius and velocity vectors

%%
% Now let's compute the chaser's state in the IJK frame, using its relative
% state in the RIC frame. We can use the 
rel_RIC = [5e-3; 0; 0; 0; 0; 0];         % [km] and [km/s] Chaser state relative to the target state, in RIC frame
Xst_RIC = blkdiag(RIC_IJK,RIC_IJK)*Xst;  % Convert target state from IJK to RIC frame
Xsc_RIC = Xst_RIC + rel_RIC;             % Compute absolute chaser state in RIC frame

Xsc = blkdiag(RIC_IJK,RIC_IJK)'*Xsc_RIC;

%%
% Given the initial target and chaser states in the IJK frame, we can form
% the full vector of states to be estimated. We can also specify the a-priori 
% covariance matrix (see Tutoral 1 for more info).

x0 = [Xsc;  % Chaser state
      Xst]; % Target state
 
P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6 ...
           1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);


%%
% We can specify which states should be estimated using the "Solve-For"
% matrix. In addition, for states that are not estimated, we can specify
% whether they shoule be considered during estimation using the "Consider"
% matrix. Here, we want to solve for all target and chaser states, so let's
% specify that in the solve-for matrix.

S = eye(12,12); % Solve-for map - solve for all 12 states
C = [];         % Consider map - no consider states

%%
% In previous tutorials, we specified a single initial state and covariance
% matrix, each of which were used as the "truth" and "estimated" initial
% states. Alternatively, ODTBX lets us specify a separate truth and
% estimated initial state and covariance. Although we will use the same
% state and covariance for both, here is how we would specify separate
% values if needed:

Xnot.Xo = x0;         % True initial state
Xnot.Xbaro = S*x0;    % Estimator initial state (note that S is the identity matrix)

Pnot.Po = P0;         % True initial covariance
Pnot.Pbaro = S*P0*S'; % Estimator initial covariance
  
%% Define Estimator Options
% In previous tutorials, we specified a single dynamics and measurement
% model for the truth and estimator. Alternatively, ODTBX lets us specify
% separate truth and estimator models in a method similar to the initial
% state and covariance above. As before, we choose to use the same truth and 
% estimated functions, but here is how we would specify different ones if
% needed:

dynfun = struct;
dynfun.tru = @dualIADyn; % True dynamics model
dynfun.est = @dualIADyn; % Estimator dynamics model

dynopts = []; % No extra options for dynamics model

datfun = struct;
datfun.tru = @dualIADat; % True measurement model
datfun.est = @dualIADat; % Estimator measurement model

%%
% Note that 'dualIADyn' is a wrapper around the individual dynamics of
% the chaser and target, and 'dualIADat' is a wrapper around the individual
% pose and GPS measurements made by the chaser.
%
% Let's specify custom measurement options. See Tutorials 1 & 2 for
% explanations of each option.

datopts = odtbxOptions('measurement');
datopts = setOdtbxOptions(datopts,'epoch',epoch);
datopts = setOdtbxOptions(datopts,'UseRange',true);
datopts = setOdtbxOptions(datopts,'UseRangeRate',false);

rsigma = 10e-3*ones(32,1); % Uncertainties in measurements (See Tutorial 2)
datopts = setOdtbxOptions(datopts,'rsigma',rsigma);

%%
% Since the datopts structure will be sent to the measurement function,
% lets add some parameters specific to the chaser's GPS antenna.

sigma = [0.1e-3 0.1e-3 0.1e-3]./3;
datopts.sig = sigma;
datopts.('AntennaPattern') = {'sensysmeas_ant.txt' 'sensysmeas_ant.txt'}; 
datopts.('AntennaPointing') = [1 -1];
datopts.('RecAcqThresh') = 19;
datopts.('YumaFile') = 'Yuma_590.txt';

%% Estimator Options
% Now let's specify options for the estimator itself. We want to disable
% measurement editing, which is done through the 'EditFlag' option.
% Had we instead wanted to accept or reject measurements based on a
% "3-sigma edit" rule, we could have passed that to ODTBX through the
% 'EditRatio' option.

opts = setOdtbxOptions('MonteCarloCases',1,'UpdateIterations',1);
opts = setOdtbxOptions(opts,'EditFlag',2*ones(35,1));  % Disable measurement editing
%opts = setOdtbxOptions(opts,'EditRatio',9*ones(35,1)); % 3-sigma edit rule for measurements


%% Sequential Estimation
% Let's use a sequential estimator to estimate the target and chaser
% spacecraft states. The estimator will automatically check for multiple
% dynamics and measurement functions and their associated options.

[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,sigsa,eflag,Pdy,Pdyt] = ...
    estseq(dynfun,datfun,tspan,Xnot,Pnot,opts,dynopts,datopts,S,C);

%% Plot Estimation Errors
% Given the estimator output, it's always of interest to plot the errors in
% the estimation process. However, the estimator outputs are all in the IJK
% frame; we first need to transform them to the RIC frame.

% The estimator outputs are cell arrays, so loop over each element of the cells
for i=length(xhat):-1:1
    % Correct the estimated states using the estimated errors
    x_tru{i} = xhat{i} - e{i};
    
    % Compute the relative position of the chaser with respect to the target
    x_rel{i} = x_tru{i}(1:3,:)-x_tru{i}(7:9,:);
    
    % Create the rotation matrix from IJK to RIC coordinates of the target
    C_IR{i} = dcm('ric',x_tru{i}(7:9,:),x_tru{i}(10:12,:));
end

S = [eye(6,6) -eye(6,6)];

% The output covariances and sensitivities also need to be rotated to the
% RIC frame. Let's build the needed rotation matrix for this.
% Rotate the stae into RIC
% Build rotation matrix that will rotate errors & covariances
% absolute/relative covariance in ric frame
% absolve/relative sensitivities in ric frame
% look at sensmos
for j=length(xhat):-1:1
    for i=size(x_tru{1},2):-1:1
        x_rel_ric{j}(:,i) = C_IR{j}(:,:,i)'*x_rel{j}(:,i);
        A = blkdiag(C_IR{j}(:,:,i),C_IR{j}(:,:,i),C_IR{j}(:,:,i),C_IR{j}(:,:,i));
        
        e_ric{j}(:,i) = A*e{j}(:,i);
        e_ric_rel{j}(:,i) = S*e_ric{j}(:,i);
        
        temp = unscrunch(P{j}(:,i));
        temp = A*temp*A';
        
        temp2 = S*temp*S';
        
        P_ric{j}(:,i) = scrunch(temp);
        P_ric_rel{j}(:,i) = scrunch(temp2);
        
        sigsa_ric(:,:,i) = A*sigsa(:,:,i)*A';
        sigsa_rel_ric(:,:,i) = S*sigsa_ric(:,:,i);
    end
end

%% 
% ODTBX provides the 'estval' function to plot estimation errors

estval(t,e_ric_rel,P_ric_rel);

%% Plot Sensitivity Mosaic
% We plot the sensitivities of the estimated data using the ODTBX 'sensmos'
% function, then customize the plot labels.

figure;
sensmos(t{end}, sigsa_rel_ric);
tickLabel = {'R'; 'I'; 'C'; 'VR'; 'VI'; 'VC'};
set(gca,'xticklabel',tickLabel,'yticklabel',tickLabel);

%% Plot Variance Sandpiles
% ODTBX provides the 'varpiles' function to plot individual contributions
% to a solved-for variable's total variance.

for i=size(Pa,3):-1:1
    Par(:,:,i) = S*Pa(:,:,i)*S';
    Pvr(:,:,i) = S*Pv(:,:,i)*S';
    Pwr(:,:,i) = S*Pw(:,:,i)*S';
    Phatar(:,:,i) = S*Phata(:,:,i)*S';
    Phatvr(:,:,i) = S*Phatv(:,:,i)*S';
    Phatwr(:,:,i) = S*Phatw(:,:,i)*S';
end

for i=size(Pa,3):-1:1
    Pat(:,:,i) = trace(Par(1:3,1:3,i));
    Pvt(:,:,i) = trace(Pvr(1:3,1:3,i));
    Pwt(:,:,i) = trace(Pwr(1:3,1:3,i));
    Phatat(:,:,i) = trace(Phatar(1:3,1:3,i));
    Phatvt(:,:,i) = trace(Phatvr(1:3,1:3,i));
    Phatwt(:,:,i) = trace(Phatwr(1:3,1:3,i));
end

dPat = Pat - Phatat;
dPvt = Pvt - Phatvt;
dPwt = Pwt - Phatwt;

Phatlin = Phatat + Phatvt + Phatwt;
P = Pat + Pvt + Pwt;
for k = 1
    figure
    varpiles(t{1},dPat(k,k,:),dPvt(k,k,:),dPwt(k,k,:),...
        Pat(k,k,:),Pvt(k,k,:),Pwt(k,k,:),...
        Phatat(k,k,:),Phatvt(k,k,:),Phatwt(k,k,:),...
        P(k,k,:),Phatlin(k,k,:));
end
titlestr = 'Relative Position Variance Sandpile [km^2]';
title(titlestr);
xlabel('Elapsed Time [s]');

%% License 
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