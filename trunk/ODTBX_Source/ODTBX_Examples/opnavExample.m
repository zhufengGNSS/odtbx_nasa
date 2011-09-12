function varargout = opnavExample(flag)
% OPNAVEXAMPLE  Example to demonstrate the functionality of OPNAVMEAS and perform regression tests.
%
% OPNAVEXAMPLE() runs OPNAVMEASEXAMPLE in demo mode.
%
% FAIL = OPNAVEXAMPLE() performs a regression test to validate OPNAVMEAS 
% based on stored measurement data and numerical calculation of measurement
% partials.
%
% keyword: measurement
% See also OPNAVMEAS, CAMERA, ATTITUDE, BODY
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

%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Kenneth Getzandanner    03/30/2011      Original opnavexample.m

% Set mode to test or demo
if nargout > 0
    
    % Running in test mode
    testMode = true;
    fail = 0;
    
    % Set the regression test tolerance
    tol = 2e-6;
    
    % Set the regression test measurement partials threshold
    thresh = 1e-20;
else
    % Running in demo mode
    testMode = false;
end

if nargin>0
    useCCD = strcmp(flag,'CCD');
else
    useCCD = 0;
end

%% BODY CLASS
% Initialize a body object to describe the observed target body

% Load relavent SPICE files to obtain planetary/asteroid states
cspice_furnsh(which('wld7802.15'))
cspice_furnsh(which('de421.bsp'))
cspice_furnsh(which('naif0009.tls'))

% Body Epoch
start = '09 MAR 2023 20:00:0.000';
epoch = cspice_str2et(start);

% SPICE ID for the asteroid 4660 Nereus
SPID = '2004660';

% Body Gravitational Parameter
GM = 1.77536e-9;

% Define the Body Semi-Axes
a = 0.2;
b = 0.15;
c = 0.1;

SA = [a b c];

% Define the Body Spin Rate & Orientation
RA = 0;
DEC = 0;
PRA = 0;
w = 2*pi/(18*3600);

SPIN = [RA DEC PRA w];

% Generate Body-Fixed Spherical Landmark Coordinates
lmk = [0 0 0 0 0 0 45 45 45 -45 -45 -45 80 -80;
       0 60 120 180 -120 -60 30 150 -90 90 -150 -30 90 -90]';

lmk = lmk*pi/180;

% Create an Instance of the Body Class
asteroid = Body(SA,SPIN,lmk,epoch,SPID,GM);

%% ATTITUDE CLASS
% Initialize an attitude object with tspan

tspan = 0:600:(18*3600);

att = Attitude(tspan);

%% INITIAL STATE
% Define the spacecraft initial state and covariance

x_kep.sma  = 0.5; %km
x_kep.ecc  = 0.01; %unitless
x_kep.incl = 70*(pi/180); %radians
x_kep.raan = 250*(pi/180); %radians
x_kep.argp = 270*(pi/180); %radians
x_kep.tran = 180*(pi/180); %radians

x0 = kep2cart(x_kep,asteroid.GM); %km & km/sec
r = x0(1:3,1);

%% CAMERA CLASS
% Initialize a camera object for measurement generation

% Define the inertial position of the camera relative to the body
Rc = r(:,1);

% Define camera field of view and focal length
FOV = 50; % degrees
f = 30e-6;

% Create an Instance of the camera class

if useCCD
    E = [1e6 1 1];
    Res = [720 720];
    
    % Create an Instance of the camera class
    ocam = CCDCamera(asteroid,Rc,eye(3,3),f,FOV*pi/180,false,E,Res);
else
    ocam = Camera(asteroid,Rc,eye(3,3),f,FOV*pi/180,false);
end

if testMode
    if useCCD
        x0 = [x0; ocam.f; ocam.Kx; ocam.Ky; ocam.s0; ocam.l0; ...
            ocam.Er; ocam.Em; ocam.En; 1e-6*ones(3,1);1e-7*ones(2,1)];
    else
        x0 = [x0; ocam.f; 1e-6*ones(3,1); 1e-7*ones(2,1)];
    end
else
    P0 = diag([50e-3 50e-3 50e-3 1e-5 1e-5 1e-5].^2);
end
%% ESTIMATOR
% Define estimator options and run estseq

% Define the dynamics function (restricted 2-body)
dynfun = @r2bp_wrap;

% Define the measurement functions
datfun.tru = @opnavmeas_tru;
datfun.est = @opnavmeas;

% Define measurement options
measopts = odtbxOptions('measurement');
measopts = setOdtbxOptions(measopts,'Camera',ocam);
measopts = setOdtbxOptions(measopts,'Attitude',att);

if testMode
    measopts = setOdtbxOptions(measopts,'CalibrationFlag', 1);
else
    measopts = setOdtbxOptions(measopts,'CalibrationFlag', 0);
end

if useCCD
    measopts = setOdtbxOptions(measopts,'OpticalSigma',1);
else
    measopts = setOdtbxOptions(measopts,'OpticalSigma',6e-10);
end

% Calculate the number of measurements
n = 2*size(lmk,1);

% Set up estimator options
opts = setOdtbxOptions('MonteCarloCases',1,'UpdateIterations',1);
opts = setOdtbxOptions(opts,'EditFlag',2*ones(n,1));
opts = setOdtbxOptions(opts,'EditRatio',9*ones(n,1));

%% TEST MODE
if testMode
    
    % Integrate trajectory
    [~,x] = integ(dynfun,tspan,x0,[],GM);
    
    options = odtbxOptions('estimator');
    options = setOdtbxOptions(options,'DatJTolerance',1e-12);
    
    % Generate measurements
    [y H] = opnavmeas(tspan,x,measopts);
    
    % Uncomment to generate regression data:
%     if useCCD
%         yref = y;
%         Href = H;  
%         save('opnavRegressionCCD','yref','Href')
%     else
%         yref = y;
%         Href = H;
%         save('opnavRegression','yref','Href')
%     end
    
    % Load test data
    if useCCD
        load opnavRegressionCCD
    else
        load opnavRegression
    end
    
    % Calculate differences
    Hdiff = Href-H;
    ydiff = yref-y;
    
    makenan = isnan(y);
    ydiff(makenan) = NaN;
    ypd = ydiff./max(abs(yref),abs(y));
    
    for i=1:size(H,3)
        Hdiff(makenan(:,i),:,i) = NaN;
    end
    Hpdiff = Hdiff./((H+Href)./2);
    
    % Check analytic partials against tolerance
    if max(max(max(Hpdiff(~isnan(Hpdiff)))))>tol
        fail = 1;
        disp(['Opnavmeas regression test failed! Analytic partials ' ...
            'do not match to desired tolerance!'])
    end
    
    % Check measurements against tolerance
    if max(max(ypd))>tol
        fail = 1;
        disp(['Opnavmeas regression test failed! Measurements ' ...
            'do not match to desired tolerance!'])
    end
    
    % Output test results
    varargout{1} = fail;
    

%% DEMO MODE
else
    % Run the sequential estimator
    [t1,~,P1,e1,dy1,Pa,Pv,Pw,Phata,Phatv,Phatw,~,~,Pdy,Pdyt] =...
        estseq(dynfun,datfun,tspan,x0,P0,opts,GM,measopts);
    
    % Plot results
    plot_results(t1,P1,e1,dy1,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt)
    
    % Plot target body
    figure(gcf+1)
    plotLmk(asteroid)
    
    % Plot camera frame
    figure(gcf+1)
    plotCameraFrame(ocam,0)
end

end

%% OPNAVMEAS_TRU
% Example wrapper function to generate an attitude profile based on the
% "true" spacecraft state
function [y H R] = opnavmeas_tru(t,x,options)

% Get the attitude object
att = getOdtbxOptions(options,'Attitude',[]);

% Set nadir pointing for all measurement times
nadir(att,t,x(1:3,:));

% Call opnavmeas
[y H R] = opnavmeas(t,x,options);

end

%% R2BP_WRAP
% Wrapper function to include derivatives of camera calibration parameters
function [xdot A Q] = r2bp_wrap(t,x,options)

n = size(x,1);
xdot = zeros(n,1);
A = zeros(n,n);
Q = zeros(n,n);

[xdotr Ar Qr] = r2bp(t,x,options);

xdot(1:6,1) = xdotr;
A(1:6,1:6) = Ar;
Q(1:6,1:6) = Qr;

end