function XNAV_only_LEO
%% LEO Orbit Determination Simulation of XNAV
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

% Original Author - Kevin Berry
% Dennis Woodfork - edited scenario for XNAV
% 2010/08/16 Sun Hur-Diaz - XNAV-only version with LEO example with ISS

%% initialize
clear;
clc;
tic; %starts the matlab timer
format long g; %forces Command Window outputs to be more precise

%% Common user modified variables
% International Space Station (ISS) states
% Source http://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html
% ISS  TRAJECTORY DATA
%  
%    Lift off time (UTC)  :  N/A
%    Area (sq ft)  :  16953.16
%    Drag Coefficient (Cd)  :   2.00
%    Monthly MSFC 50% solar flux (F10.7-jansky)  :   74.8
%    Monthly MSFC 50% earth geomagnetic index (Kp)  :   0.74
%    ET - UTC  (sec) :   66.18
%    UT1 - UTC (sec) :    0.00
%
%     Coasting Arc #1 (Orbit 3217)
%     ---------------------------------------
%  
%     Vector Time (GMT): 2010/223/12:00:00.000 (Aug. 11)
%     Vector Time (MET): N/A
%     Weight (LBS)     : 828730.9
%
%               M50 Cartesian                         J2K Cartesian  
%     -----------------------------------       --------------------------------
%      X    =       17242850.30                 X    =         5263998.56
%      Y    =       -7583853.79  feet           Y    =        -2252740.29  meter
%      Z    =       11535528.63                 Z    =         3541580.75
%      XDOT =       -804.740668                 XDOT =        -336.973421
%      YDOT =      20539.607630  feet/sec       YDOT =        6257.216826  meter/sec
%      ZDOT =      14656.538506                 ZDOT =        4465.898658
%
wt_lb = 828730.9;
area_sqft = 16953.16;
EpochStr = '11 Aug 2010 12:00:00';
epoch    = datenum(EpochStr);
x0.Xbaro = [5263998.56;-2252740.29;3541580.75;-336.973421;6257.216826;4465.898658]*1e-3; % km & km/sec
P0.Pbaro = diag([1e-1 1e-1 1e-1 1e-4 1e-4 1e-4].^2); %initial covariance matrix (km^2 & km^2/sec^2)
cR = 1.8; % initial Coefficient of Reflectivity
mu = 3.986004415e5; % km^3/sec^2
Pis = sqrt((x0.Xbaro(1))^2+(x0.Xbaro(2))^2+(x0.Xbaro(3))^2);
T_orb = 2*pi*sqrt(Pis^3/mu);
fprintf('\n Orbit period is %10.8g sec.\n',T_orb)
mass = wt_lb*0.45359237;
area_sqm = area_sqft * 0.09290304;

UpdateIterations = 1; % The number of times to iterate the measurement update
MonteCarloCases  = 1; % Number of Monte Carlo runs to be performed

%% Set the Consider states
x0.Xo    = x0.Xbaro;
P0.Po    = P0.Pbaro;

% Set the S and C matrices, which let estseq know what is a consider state
 S=eye(6); 
%  C=[]; % Empty matrix is the default

%% Estimator options
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'UpdateIterations',UpdateIterations,...
    'MonteCarloCases',MonteCarloCases,...
    'refint',3,...
    'MonteCarloSeed',0,...%    'EditFlag',2,...
    'OdeSolver',@ode45,...
    'OdeSolvOpts', odeset('reltol',1e-7,'abstol',1e-7));

%% Set up the options and java object for JAT forces
jOptions    = odtbxOptions('force');
jOptions    = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jOptions    = setOdtbxOptions(jOptions, 'cD', 2.0); %coeff of drag
jOptions    = setOdtbxOptions(jOptions, 'cR', cR); %coeff of reflectivity
jOptions    = setOdtbxOptions(jOptions, 'mass', mass ); %kg
jOptions    = setOdtbxOptions(jOptions, 'dragArea', area_sqm  , 'srpArea', area_sqm); %m^2
jOptions    = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM2');
jOptions    = setOdtbxOptions(jOptions, 'gravDeg', 8, 'gravOrder', 8); %max 20x20
jOptions    = setOdtbxOptions(jOptions, 'useSolarGravity', true);
jOptions    = setOdtbxOptions(jOptions, 'useLunarGravity', true);
jOptions    = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', true);
jOptions    = setOdtbxOptions(jOptions, 'useAtmosphericDrag', true);
jOptions    = setOdtbxOptions(jOptions, 'atmosphereModel', 'HP');
jOptions    = setOdtbxOptions(jOptions, 'nParam', 4.3); %inclination param for HP model

% Create a java object that stores the relevant information for 
% propagating Earth-centric orbits using JAT force models.

jatWorldtru = createJATWorld(jOptions); 
jatWorldest = createJATWorld(jOptions); %creates a java object that stores the
% 
%% Set the measurement options
 measOptions = odtbxOptions('measurement');

% The .tru values are the truth dynamics and measurements
% The .est values are the estimated dynamics and measurements
dynfun.tru = @jatForces_km_tru; % See function below
dynfun.est = @jatForces_km_est; % See function below
datfun.tru = @xnavmeas;
datfun.est = @xnavmeas;
dynarg.tru = jatWorldtru;
dynarg.est = jatWorldest;

t_obs       = 1e5; % sec, duration of each observation used in noise model       

% Retrieive pulsar data and create fps options structure and add fps 
% specific fields
pulsar_id = {'B0531+21','B1821-24','B1937+21',};
pulsar_database = 'pulsars.txt';
[rso_eci,C_1] = getPulsarData(pulsar_database,pulsar_id);

xnavOpts.rso_eci = rso_eci;
xnavOpts.C_1 = C_1;
xnavOpts.obs_time = t_obs;
xnavOpts.epoch = epoch;
xnavOpts.num_meas = 3;
xnavOpts.sigma_r = [];
xnavOpts.sigma_rr = [];
xnavOpts.useRangeRate = true;

measOptions = setOdtbxOptions(measOptions,'xnav',xnavOpts);

%% Setting up the measurement time vector and schedule
%
% The measurements from the three pulsars are staggered assuming a single 
% detector is time-shared.  We force measurements from only one pulsar at 
% a time by specifying schedules that are exclusive of each other.  This 
% is done by making each schedule duration very short.

t_interval  = 3e5; % sec, time between measurements for a given pulsar
tfinal      = 3*t_interval;

% The start times of the measurements for the three pulsars
tstart = [0 1e5 2e5];

% Set the measurement times of each pulsar
Pulsar1_times = tstart(1)+0:t_interval:tfinal-t_interval;
Pulsar2_times = tstart(2)+0:t_interval:tfinal-t_interval;
Pulsar3_times = tstart(3)+0:t_interval:tfinal-t_interval;

% Concatenate and sort to form the measurement time vecor
tspan = sort([Pulsar1_times(:);Pulsar2_times(:);Pulsar3_times(:)]);

% Set up measurement schedules for each pulsar.  The end times are one
% second later than the measurement times to ensure no overlap
Sched1 = [1*ones(length(Pulsar1_times),1) Pulsar1_times(:) Pulsar1_times(:)+1];
Sched2 = [2*ones(length(Pulsar2_times),1) Pulsar2_times(:) Pulsar2_times(:)+1];
Sched3 = [3*ones(length(Pulsar3_times),1) Pulsar3_times(:) Pulsar3_times(:)+1];
Sched = [Sched1;Sched2;Sched3];
measOptions = setOdtbxOptions(measOptions,'Schedule',Sched); % tracking schedule

datarg.tru = measOptions;
datarg.est = datarg.tru;

% Setting different measurement covariance for estimator
% -------------------------------------------------------
t1 = 0;
% Set a temporary blank Schedule to get sigmas for all pulsars
measOptions_nosched = setOdtbxOptions(measOptions,'Schedule',[]);
[~, ~, R1] = xnavmeas(t1,x0.Xbaro,measOptions_nosched);
sigma_R1 = diag(sqrt(R1));
if xnavOpts.useRangeRate == true
    sigma_r = sigma_R1(1:2:end);
    sigma_rr = sigma_R1(2:2:end);
else
    sigma_r = sigma_R1;
end
xnavOpts.sigma_r = .1*sigma_r;
xnavOpts.sigma_rr = .1*sigma_rr;
datarg.est  = setOdtbxOptions(datarg.est,'xnav',xnavOpts);

%% Run the sequential estimator 

[t,xhat,Phat,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,~,eflag,Pdy]=estseq(dynfun,datfun,tspan,x0,P0,eOpts,dynarg,datarg);
[Phat1,e1,Pa1,Pv1,Pw1,Phata1,Phatv1,Phatw1]=rot2ric(t,xhat,Phat,e,Pa,Pv,Pw,S,Phata,Phatv,Phatw);

runtime = toc %#ok<NASGU,NOPRT> % End the Matlab timer

disp('You may ''save'' the results of the filter simulation at this time.')
disp('Type ''return'' to continue with plotting.')
keyboard
%% Plot the Ensemble Statistics, Sandpiles
for k = length(Pa1):-1:1,
    SPaSt(:,:,k) = S*Pa1(:,:,k)*S'; 
    SPvSt(:,:,k) = S*Pv1(:,:,k)*S'; 
    SPwSt(:,:,k) = S*Pw1(:,:,k)*S'; 
end
SPSt = SPaSt + SPvSt + SPwSt;
Ptrue = scrunch(SPSt);
estval(t,e1,Phat1,Ptrue,0)
dPa = SPaSt - Phata1;
dPv = SPvSt - Phatv1;
dPw = SPwSt - Phatw1;
Phatlin = Phata1 + Phatv1 + Phatw1;
for k = 1:size(S,1),
    figure
    varpiles(t{1},dPa(k,k,:),dPv(k,k,:),dPw(k,k,:),...
        SPaSt(k,k,:),SPvSt(k,k,:),SPwSt(k,k,:),...
        Phata1(k,k,:),Phatv1(k,k,:),Phatw1(k,k,:),SPSt(k,k,:),Phatlin(k,k,:));
end

%% Plot 3 Sigma Position and Velocity Errors
PosErr=zeros(MonteCarloCases,length(t{1}));
VelErr=zeros(MonteCarloCases,length(t{1}));
for n=1:MonteCarloCases
PosErr(n,:) = sqrt([1 1 1 0 0 0]*e{n}.^2);
VelErr(n,:) = sqrt([0 0 0 1 1 1]*e{n}.^2);
end
figure; hold on;
title('RSS Position Errors')
H1 = plot(t{1}/3600,mean(PosErr)+3*std(PosErr),'r','LineWidth',3);
H2 = plot(t{1}/3600,mean(PosErr),'b','LineWidth',3);
H3 = plot(t{1}/3600,PosErr','k'); %plots all of the position error values (km)
xlabel('Hours from Epoch');ylabel('km')
legend([H1(1),H2(1),H3(1)],'+3sigma Error','Mean Error','Raw Errors','Location','best');
hold off
figure; hold on;
title('RSS Velocity Errors')
H1 = plot(t{1}/3600,mean(VelErr)+3*std(VelErr),'r','LineWidth',3);
H2 = plot(t{1}/3600,mean(VelErr),'b','LineWidth',3);
H3 = plot(t{1}/3600,VelErr','k'); %plots all of the velocity error values (km/s)
xlabel('Hours from Epoch');ylabel('km/s')
legend([H1(1),H2(1),H3(1)],'+3sigma Error','Mean Error','Raw Errors','Location','best');
hold off

%% Plot all xhats with the Earth
figure; hold on;
for n=1:MonteCarloCases
    xhat_n = xhat{n};
    plot3(xhat_n(1,:),xhat_n(2,:),xhat_n(3,:),'.');
end
xlabel('X (Km)');ylabel('Y (Km)');zlabel('Z (Km)');
title('Estimated Trajectories');
[X,Y,Z] = sphere(20);surf(6378*X,6378*Y,6378*Z);
colormap( [(0) (.7) (1)]);axis equal;

%% Plot the Range and Range Rate Residuals for each Monte Carlo run
for n=1:MonteCarloCases
    t_n    = t{n};
    dy_n   = dy{n};
    Puls_lab = {'Pulsar 1';'Pulsar 2';'Pulsar 3';'Pulsar 4'};
    if xnavOpts.useRangeRate == true
        figure;
        plot(t_n/3600,dy_n(1:2:end,:)','Linestyle','none','Marker','.','MarkerSize',15);
        legend(Puls_lab(1:xnavOpts.num_meas),'Location','southOutside');
        xlabel('Hours');ylabel('km');title('Range residuals');
        figure;
        plot(t_n/3600,dy_n(2:2:end,:)','Linestyle','none','Marker','.','MarkerSize',15);
        legend(Puls_lab(1:xnavOpts.num_meas),'Location','southOutside');
        xlabel('Hours');ylabel('km/sec');title('Range Rate residuals');
    else
        figure;
        plot(t_n/3600,dy_n','Linestyle','none','Marker','.','MarkerSize',15);
        legend(Puls_lab(1:xnavOpts.num_meas),'Location','southOutside');
        xlabel('Hours');ylabel('km');title('Range residuals');
    end
end

%% Plot measurement innovations
plot_ominusc(t,dy,Pdy,[],[],[],eflag)

%% Wrappers for the dynamics functions to allow changing th eprocess noise
function [xDot,A,Q] = jatForces_km_est(t,x,jatWorld)
[xDot,A,Q] = jatForces_km(t,x,jatWorld);
Q=Q*5e10;

function [xDot,A,Q] = jatForces_km_tru(t,x,jatWorld)
[xDot,A,Q] = jatForces_km(t,x,jatWorld);
Q=Q*1e0;

%% Transforms from ECI to Radial-Intrack-Crosstrack
function [varargout]=rot2ric(varargin)
% [P,e,Pa,Pv,Pw,Phata,Phatv,Phatw]=rot2ric(t,x,P,e,Pa,Pv,Pw,S,Phata,Phatv,Phatw)
t=varargin{1};
x=varargin{2};
P=varargin{3};
if nargin >= 4
    e=varargin{4};
    if nargin >= 5
        Pa=varargin{5};
        Pv=varargin{6};
        Pw=varargin{7};
        S=varargin{8};
        if nargin >= 9
            Phata=varargin{9};
            Phatv=varargin{10};
            Phatw=varargin{11};
        end
    end
end

if exist('S','var') && isempty(S)
    S=eye(size(Pa));
end

for el = 1:length(t),
    for k = 1:length(t{el}),
        M = dcm('ric',x{el}(1:3,k),x{el}(4:6,k));
        M = blkdiag(M,M);
        phatel = unscrunch(P{el}(:,k));
        phatel(1:6,1:6) = M*phatel(1:6,1:6)*M';
        phatel = (phatel + phatel')/2;
        P{el}(:,k) = scrunch(phatel);%*1e6;
        if nargin >= 4
            e{el}(1:6,k)=M*e{el}(1:6,k);
            if nargin >= 5
                Pa(:,:,k) = S*Pa(:,:,k)*S'; 
                Pv(:,:,k) = S*Pv(:,:,k)*S'; 
                Pw(:,:,k) = S*Pw(:,:,k)*S';
                Pa(1:6,1:6,k) = M*Pa(1:6,1:6,k)*M';
                Pv(1:6,1:6,k) = M*Pv(1:6,1:6,k)*M';
                Pw(1:6,1:6,k) = M*Pw(1:6,1:6,k)*M';
                if nargin >=9
                    Phata(1:6,1:6,k) = M*Phata(1:6,1:6,k)*M';
                    Phatv(1:6,1:6,k) = M*Phatv(1:6,1:6,k)*M';
                    Phatw(1:6,1:6,k) = M*Phatw(1:6,1:6,k)*M';
                end
            end
        end
    end
end

if nargout>=1
    varargout{1}=P;
    if nargout>=2
        varargout{2}=e;
        if nargout >= 3
            varargout{3}=Pa;
            varargout{4}=Pv;
            varargout{5}=Pw;
            if nargout >= 6
                varargout{6}=Phata;
                varargout{7}=Phatv;
                varargout{8}=Phatw;
            end
        end
    end
end
