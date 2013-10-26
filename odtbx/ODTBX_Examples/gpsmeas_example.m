% An example script that demonstrates the use of the GPS measurement model, gpsmeas.
%
% Calculates and plots range and range rate measurements at each time point
% for a GPS receiver in a specified Earth-Centered, Inertial (ECI)
% trajectory from each GPS satellite in a given GPS almanac.  All of the
% parameters of the transmit-receive system link budget are specified
% as input options to gpsmeas.
%
% See also: gpsmeas, jatDCM, JATConstant
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

close all


%% Specifiy input files

% Example YUMA file describing the GPS constellation
yuma_file = 'Yuma1134_example.txt';

% Specify the appropriate antenna pattern settings.
ant_pat{1} = 'sensysmeas_ant.txt';

% Specify the appropriate antenna pointing profile.  [1] is a zenith
% pointing antenna profile.  See gpsmeas.m's AntennaPointing option.
ant_point = 1;

%% Configure the receiver's trajectory and simulation time

% Load a trajectory file.  This not only provides data on the satellite
% position, velocity, and time, it contains additional data on the transmit
% and receive systems.
truth = load('gpsmeas_testtraj.mat');

%  Provide the epoch as a datenum
% epoch = datenum(truth.sim_start_utc);
epoch = datenum([2013 1 1 0 0 0]);

% Create the time vector using the same time intervals, referenced from the
% first time
t = truth.time - truth.time(1);


% Set up satellite state in Km and Km/s
% This truth data is in the Earth-Centered, Earth-Fixed (ECEF) frame.
xecef(1:3,:) = truth.sat_pos/1000;
xecef(4:6,:) = truth.sat_vel/1000;


% Convert true position and velocity from ECEF to Earth-Centered, Inertial
% (ECI).
ecef2eci  = jatDCM('ecef2eci', (t/86400)+epoch);    % Generate rotation matrices (3x3xn)
w         = [0; 0; JATConstant('wEarth')];          % Rotation rate of earth

% Create temporary vectors to hold the state
tmp_sat_pos = zeros(3, length(t));
tmp_sat_vel = zeros(3, length(t));

% Perform the rotations to the ECI frame
for n=1:length(t)
    tmp_sat_pos(1:3,n) = ecef2eci(:,:,n)*xecef(1:3,n);
    tmp_sat_vel(1:3,n) = ecef2eci(:,:,n)*(cross(w,xecef(1:3,n))+xecef(4:6,n));
end

% Save off the new trajectory, place position and velocity into the same
% ECI state vector.
x(1:3,:) = tmp_sat_pos; % km, eci
x(4:6,:) = tmp_sat_vel; % km/s, eci

%% Set the options for the call

% gpsmeas uses atmosphere height above earth for 'AtmosphereMask' below,
% gpstools uses a total distance from the earth's center.  This is how
% gpsmeas defines the earth radius and adds it to 'AtmosphereMask'.

EARTH_RADIUS = JATConstant('rEarth','WGS84') / 1000;  % km Equatorial radius of Earth

% measOptions is the structure that specifies the options required by
% gpsmeas.m.  See gpsmeas.m for information on each of these options.
% Note, only two options are not set in this example.
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions, 'epoch', epoch);                  % Time associated with start of simulation, datenum format
measOptions = setOdtbxOptions(measOptions, 'useRange', true);                   % Compute range measurements
measOptions = setOdtbxOptions(measOptions, 'useRangeRate', true);                   % Compute range rate measurements
measOptions = setOdtbxOptions(measOptions, 'useLightTime', true);               % Move GPS Position to account for travel time

% The rSigma parameter is not demonstrated in this example, this parameter
% is only required when calculating the measurement covariance, R, as in:
% [y, H, R] = gpsmeas(t, x, measOptions);
linkbudget.GPSBand           = 'L1';                   % See truth.freq
linkbudget.YumaFile          = yuma_file;              % Input YUMA file
% The Rotation2ECI parameter is not demonstrated in this example, this
% parameter is only required when providing satellite state in another
% coordinate system other than ECI.
linkbudget.AntennaPointing   = ant_point;              % Specify attitude profile for each antenna
linkbudget.AntennaPattern    = ant_pat;                % Specify receive antenna pattern for each antenna
linkbudget.AntennaMask       = truth.rcv_ant_mask;     % Cut off angle for the receive antenna
linkbudget.AtmosphereMask    = truth.r_mask/1000 - EARTH_RADIUS; % Mask altitude, km
linkbudget.NoiseTemp         = truth.Ts;               % System noise temp of receive antenna (K)
linkbudget.AtmAttenuation    = truth.Ae;               % Attenuation due to atmosphere, should be negative (dB)
linkbudget.TransPowerLevel   = truth.sv_power;         % Transmitter power level (1=min, 2=typical, 3=max)
linkbudget.TransPowerOffset  = truth.xmit_power_offset;% Global transmitter power offset (dB)
linkbudget.GPSBlock          = truth.sv_block;         % GPS Satellite Block  (1-II/IIA, 2-IIR, 3-IIR-M, 4-IIF)
linkbudget.TransAntMask      = truth.xmit_ant_mask;    % Cut off angle for the transmit antenna (rad)
linkbudget.ReceiverNoise     = truth.Nf;               % Noise figure of receiver/LNA (dB)
linkbudget.RecConversionLoss = truth.L;                % Receiver implementation, A/D conversion losses (dB)
linkbudget.SystemLoss        = truth.As;               % System losses, in front of LNA (dB)
linkbudget.LNAGain           = truth.Ga;               % Gain provided by the LNA (dB)
linkbudget.CableLoss         = truth.Ac;               % Cable losses after LNA (dB)
linkbudget.RecAcqThresh      = truth.CN0_lim;          % Receiver acquisition threshold (dB-Hz)
linkbudget.RecTrackThresh    = truth.CN0_lim;          % Receiver tracking threshold (dB-Hz)
linkbudget.DynamicTrackRange = truth.dyn_range;        % Receiver dynamic range (dB)
linkbudget.PrecnNutnExpire   = 0.1;                    % Length of time an Earth rotation parameter set is valid (days)

measOptions = setOdtbxOptions(measOptions, 'linkbudget', linkbudget);



%% GPSMeas function call, calculating measurements only

% [y,H,R,AntLB,dtsv] = gpsmeas(t, x, measOptions);

% When different GPS Satellite Block Types are used per PRN, gpsmeas()
% needs to be called satellite by satellite.  
y=[];
qatt=[];
for sv=1:32
    GeneratingPRN=sv
    [ym,H,R,AntLB,dtsv] = gpsmeas(t, x, measOptions,qatt,sv);
    y=[y;ym];
end


%% Plot the results
% The output measurements, y, are controlled by the 'useRange' and
% 'useRangeRate' input options.  If both are enabled, then the measurements
% are interleaved as column vectors of range then range rate for each G
% satellite at a given time.  If there are 32 satellites then there will be
% 64 outputs.  Also note that NaN are possible if there are geometric or
% signal constraints.
r = 1:2:64; % Pick out the range measurements from the data
d = 2:2:64; % Pick out the range rate measurements from the data

figure;
subplot(2,1,1),plot(y(r,:)','.');
title('gpsmeas() Range Measurements  (per GPS satellite)');
ylabel('Range (Km)');
subplot(2,1,2),plot(y(d,:)','.');
title('gpsmeas() Range Rate Measurements (per GPS satellite)');
ylabel('Range Rate (Km/s');




