function fail = restartRegression
% RESTARTREGRESSION  Restart Record and Impulsive Burn Regression Test.
%
%   [fail] = RESTARTREGRESSION regression test to verify results from the
%   restart record implementation in ESTSEQ and the IMPULSIVEBURN function.
%   The test compares data from a simple Hohmann Transfer example.
%
%   keyword: Impulsive, Maneuver, Restart, Regression
%
%   See also
%      Estimation:      ESTSEQ
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

% Kenneth M. Getzandanner
% NASA Goddard Space Flight Center
%
% Modification History
% ---------------------
%

%% Initialize Variables

% Initialize fail to "pass"
fail = 0;

% Define initial orbit it Keplarian Elements
KOE.sma = 6778;
KOE.ecc = 0;
KOE.incl = 0;
KOE.raan = 0;
KOE.argp = 0;
KOE.tran = 0;

% Define dynamics and measurment model functions/arguments
dynfun = @r2bp;
datfun = @rwdat;
dynarg = [];
datarg = 1e-3;

% Set up ODTBX options structure
opts = setOdtbxOptions('MonteCarloCases',1,'UpdateIterations',1);
opts = setOdtbxOptions(opts,'EditFlag',1*ones(14,1));
opts = setOdtbxOptions(opts,'EditRatio',9*ones(14,1));
opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);

% Convert KOE to cartesian for initial conditions
x0 = kep2cart(KOE,3.986004415e5);

% Assign initial covariance
P0 = 1e-3*eye(6);

% Define orbit segment durations
T1 = 5555;
TH = 8758.47;
T2 = 33315;

% Define timestep
dt = 600;

% Define impulsive maneuvers
dv1 = 1.832;
dv2 = 1.342;

% Set up time arrays
tspan1 = 0:dt:T1;
tspan2 = tspan1(end):dt:(tspan1(end)+TH);
tspan3 = tspan2(end):dt:(tspan2(end)+T2);

%% Initial Orbit
%
[t1,xhat1,P1,e1,dy1,Pa1,Pv1,Pw1,Phata1,Phatv1,Phatw1,~,~,~,...
    ~,Pm1,Phatm1,restartRecord2]= estseq(dynfun,datfun,tspan1,x0,...
    P0,opts,dynarg,datarg);

% Assign DV values
dv(:,1) = unit(xhat1{1}(4:6,end))*dv1;
dv(:,2) = dv(:,1);

% Perform Impulsive Maneuver
restartRecord2 = impulsiveBurn(restartRecord2,dv,0,0);

%% Transfer Orbit
%
[t2,xhat2,P2,e2,dy2,Pa2,Pv2,Pw2,Phata2,Phatv2,Phatw2,~,~,~,...
    ~,Pm2,Phatm2,restartRecord3]= estseq(restartRecord2,tspan2);

% Assign DV values
dv(:,1) = unit(xhat2{1}(4:6,end))*dv2;
dv(:,2) = dv(:,1);

% Perform Impulsive Maneuver
restartRecord3 = impulsiveBurn(restartRecord3,dv,0,0);

%% Final Orbit
%
[t3,xhat3,P3,e3,dy3,Pa3,Pv3,Pw3,Phata3,Phatv3,Phatw3,~,~,~,...
    ~,Pm3,Phatm3]= estseq(restartRecord3,tspan3);

%% Plot Trajectory
%
plot(xhat1{1}(1,:),xhat1{1}(2,:),'b',...
    xhat2{1}(1,:),xhat2{1}(2,:),'r',...
    xhat3{1}(1,:),xhat3{1}(2,:),'b')

axis equal
title('Hohmann Transfer Example')
xlabel('X [km]')
ylabel('Y [km]')

%% Compare Results
%
t = [t1{1} t2{1} t3{1}];
e = [e1{1} e2{1} e3{1}];
P = [P1{1} P2{1} P3{1}];
xhat = [xhat1{1} xhat2{1} xhat3{1}];
dy = [dy1{1} dy2{1} dy3{1}];

Pa = cat(3,Pa1,Pa2,Pa3);
Pv = cat(3,Pv1,Pv2,Pv3);
Pw = cat(3,Pw1,Pw2,Pw3);
Pm = cat(3,Pm1,Pm2,Pm3);

P_LC = Pa+Pv+Pw+Pm;

Phata = cat(3,Phata1,Phata2,Phata3);
Phatv = cat(3,Phatv1,Phatv2,Phatv3);
Phatw = cat(3,Phatw1,Phatw2,Phatw3);
Phatm = cat(3,Phatm1,Phatm2,Phatm3);

Phat_LC = Phata+Phatv+Phatw+Phatm;

% Uncomment to save regression test data:
% t_test = t;
% e_test = e;
% xhat_test = xhat;
% P_test = P;
% P_LC_test = P_LC;
% Phat_LC_test = Phat_LC;
% dy_test = dy;
% 
% save('restartRecordData.mat','*test')

load('restartRecordData.mat')

dmax = 3e-8; % Note: set at 3e-8 to obtain passing results on 32-bit 
             % platforms.  5e-9 passed on 64 bit platforms.

if(any(t-t_test))
    disp('Restart Regression Test Failed: t')
    fail = 1;
end

if(max(max(abs(xhat-xhat_test)))>dmax)
    disp('Restart Regression Test Failed: xhat')
    fail = 1;
end

if(max(max(abs(e-e_test)))>dmax)
    disp('Restart Regression Test Failed: e')
    fail = 1;
end

if(max(max(abs(P-P_test)))>dmax)
    disp('Restart Regression Test Failed: P')
    fail = 1;
end

if(max(max(abs(P_LC-P_LC_test))))
    disp('Restart Regression Test Failed: P (LinCov)')
    fail = 1;
end

if(max(max(abs(Phat_LC-Phat_LC_test))))
    disp('Restart Regression Test Failed: Phat (LinCov)')
    fail = 1;
end

if(max(max(abs(dy-dy_test)))>dmax)
    disp('Restart Regression Test Failed: dy')
    fail = 1;
end

end

%% Simple Measurement Function
% Returns the state at time t
function [Y,H,R] = rwdat(t,X,r)
Y = X;
H = repmat(eye(6),[1,1,length(t)]);
R = r*repmat(eye(6),[1,1,length(t)]);
end