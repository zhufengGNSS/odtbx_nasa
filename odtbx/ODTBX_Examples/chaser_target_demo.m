% Demo showing the estimation of relative position in a chaser-target situation.

%% Chaser-Target Demo
% This demo shows the estimation of relative position in a chaser-target
% situation, using measurements and gps on the chaser.
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
%  Original credit for this example goes to K. Getzandanner.
%  
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Phillip Anderson        08/30/2012      Created from original demo
%   & Ravi Mathur                           file demo2.m

echo on
%
% ************************************************************************
% This demo shows the estimation of relative position in a chaser-target
% situation, using measurements and gps on the chaser.
% ************************************************************************

%
%% Define Initial State & Covariance
% Define all scenario-specific variables and ODTBX options.
%

% Define the Target Vehicle state and convert to a direction cosine matrix.
Xss = [5426.510869339219;
       -41803.696487556008;
       -240.913638515624;
       3.049219473653;
       0.396422336331;
       -0.053147501025];
   
DCM = dcm('ric',Xss(1:3),Xss(4:6)); % ODTBX fcn to make DC matrix

% Define the relative Chase Vehicle state.
rel_RIC = [5e-3 0 0 0 0 0]';
Xss_RIC = blkdiag(DCM,DCM)'*Xss;
Xsc_RIC = Xss_RIC + rel_RIC;

Xsc = blkdiag(DCM,DCM)*Xsc_RIC;

% Combine chase and target vehicle states to create the initial reference
% state. Also define the a priori covariance matrix.
x0 = [Xsc;
      Xss];
 
P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6 ...
           1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);


% Set up the solve-for and consider mapping
I = eye(6,6);   
O = zeros(6,6); 

S = eye(12,12); % Solve-for map - solve for all 12 states
C = [];         % Consider map - no consider states

% Use the same a priori for both truth and estimator
Xnot.Xo = x0;         % True initial state
Xnot.Xbaro = S*x0;    % Estimator initial state

Pnot.Po = P0;         % True initial covariance
Pnot.Pbaro = S*P0*S'; % Estimator initial covariance

%
% Define the times and epoch to process measurements and estimate
tspan = 0:10:300;
epoch = datenum('13 DEC 2010 00:00:0.000');
  
%
%% Define estimator-specific options
%

% Use the same dynamics model for the truth and the estimator.
dynfun.tru = @dualIADyn;
dynfun.est = @dualIADyn;

dynopts = [];

%
% Define the data model. Here we use the same model for the truth 
% and the estimator.

datfun.tru = @dualIADat;
datfun.est = @dualIADat;

% Measurement uncertainty
rsigma = 10e-3*ones(32,1);
sigma = [0.1e-3 0.1e-3 0.1e-3]./3;

% Measurement options
measopts = odtbxOptions('measurement');
measopts = setOdtbxOptions(measopts,'epoch',epoch);
measopts = setOdtbxOptions(measopts,'UseRange',true);
measopts = setOdtbxOptions(measopts,'UseRangeRate',false);
measopts = setOdtbxOptions(measopts,'rsigma',rsigma);

measopts.('AntennaPattern') = {'sensysmeas_ant.txt' 'sensysmeas_ant.txt'}; 
measopts.('AntennaPointing') = [1 -1];
measopts.('RecAcqThresh') = 19;
measopts.('YumaFile') = 'Yuma_590.txt';
measopts.sig = sigma;

%
%% Estimator Options

opts = setOdtbxOptions('MonteCarloCases',1,'UpdateIterations',1);
opts = setOdtbxOptions(opts,'EditFlag',2*ones(35,1));
opts = setOdtbxOptions(opts,'EditRatio',9*ones(35,1));
% opts = setOdtbxOptions(opts,'Refint',-3);


%
%% Time the Sequential Estimator

tic

% myest = estseq(dynfun,datfun,tspan,Xnot,Pnot,opts,dynopts,measopts,S,C);
% 
% [t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,sigsa,eflag,Pdy,Pdyt] = myest.run_estimator()

mysim = est_control('estnew',dynfun,datfun,tspan,Xnot,Pnot,opts,dynopts,measopts);
[t,xhat,P,e,Y] = mysim.run_sim();

toc

%
%% Rotate to RIC

for i=length(xhat):-1:1
    x_tru{i} = xhat{i} - e{i};
    x_rel{i} = x_tru{i}(1:3,:)-x_tru{i}(7:9,:);
    C_IR{i} = dcm('ric',x_tru{i}(7:9,:),x_tru{i}(10:12,:));
end

S = [eye(6,6) -eye(6,6)];

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

%
%% Plot estimation errors using estval

estval(t,e_ric_rel,P_ric_rel)

%
%% Plot sensitivity mosaic

ns = 6;
n = 12;
lentr = size(sigsa,3);

lastfig = ns;
figure
hold off
pcolor(eye(ns+1,ns)*10*log10(abs(sigsa_rel_ric(:,:,end)))*eye(n,n+1))
tickLabel = {'R'; 'I'; 'C'; 'VR'; 'VI'; 'VC'};
set(gca,'xtick',(1:n)+0.5,'ytick',(1:ns)+0.5,...
    'xticklabel',tickLabel,'yticklabel',tickLabel)
xlabel('{\it A Priori} State Index')
ylabel('Relative State Index')
set(gca,'ydir','rev','xaxisloc','top')
axis equal
colorbar
title('Sensitivity Mosaic')
hold on

%
%% Plot Varpiles

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
title(titlestr)
xlabel('Elapsed Time [s]')

echo off