% Support script to gps_antenna_tools_demo.m
%
% gps_antenna_tools_support.m
%
% NOTE: This script is primarily for generating example GPS raw measurement
% data similar to what would be generated from a GPS receiver.  Users
% running the gps_antenna_tools_demo are not expected to perform any of
% these steps when they are using the GPS Enhancements toos to analyze
% antenna gains and the GPS link budget.
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

%% We specify our input data files
%
% First, the data file with most of our input data for receiver states,
% etc.
truth_file = 'gps_ant_tools.mat';

% We'll have a zenith-pointing antenna for 
ant_point = 1;
ant_rx_pat = '{truth.rec_pattern}'; % 3D antenna

% Load the truth data file.
truth = load(truth_file);

%% set up time
% gpstools time is GPS seconds
% gpsmeas time is seconds from epoch, with an epoch in UTC
epoch = datenum(truth.sim_start_utc);

% create the time vector using the same time intervals, referenced from the
% first time
t = truth.time - truth.time(1);

%% set up satellite state
% gpstools state is in m (ecef for leo case)
% gpsmeas state is in km, default to eci unless a rotation function is
% supplied
xecef = zeros(6,size(truth.sat_pos,2));
xecef(1:3,:) = truth.sat_pos/1000;
xecef(4:6,:) = truth.sat_vel/1000;

%% Convert true position and velocity from ECEF to ECI
ecef2eci  = jatDCM('ecef2eci', (t/86400)+epoch);
w         = [0; 0; JATConstant('wEarth')];
% change name, added tmp_
tmp_sat_pos = zeros(3, length(t));
tmp_sat_vel = zeros(3, length(t));
for n=1:length(t)
    tmp_sat_pos(1:3,n) = ecef2eci(:,:,n)*xecef(1:3,n);
    tmp_sat_vel(1:3,n) = ecef2eci(:,:,n)*(cross(w,xecef(1:3,n))+xecef(4:6,n));
end
x = zeros(6,length(t));
x(1:3,:) = tmp_sat_pos; % km, eci
x(4:6,:) = tmp_sat_vel; % km/s, eci
% (end of state conversion)

%% Compute attitude quaternion assuming zenith pointing
% Note that the antenna z-axis should be pointing zenith.  The body to
% antenna rotation is done through the measOptions belows with the
% parameter AntennaOrientation.
% The following is used only if the receiver antenna is 2-D.
R = dcm('ric', x(1:3,:),x(4:6,:)); % Rotation from ECI to RIC
qatt = dcm2q(R);

%% setup for gpsmeas

% gpsmeas uses atmosphere height above earth for 'AtmosphereMask' below,
% gpstools uses a total distance from the earth's center.  This is how
% gpsmeas defines the earth radius and adds it to 'AtmosphereMask'.
EARTH_RADIUS = JATConstant('rEarth','WGS84') / 1000;  % km Equatorial radius of Earth

measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions, ...
    'epoch', epoch, ...
    'useRange', true, ...
    'useRangeRate', false, ...
    'useDoppler', false, ...
    'YumaFile', yuma_file, ...
    'AntennaPointing', ant_point, ... % override of truth.pointing_ref
    'useIonosphere', false, ...
    'useTroposphere', false, ...
    'PrecnNutnExpire', 0.1, ... % days
    'AntennaOrientation', dcm('ax2', pi/2));

%rSigma (use default, no analogue in gpsdef)
%Rotation2ECI (use default, no analogue in gpsdef)

link_budget.GPSBand           = 'L1';                             % see truth.freq
link_budget.AntennaPattern    = eval(ant_rx_pat);
link_budget.RXAntennaMask     = truth.rcv_ant_mask;
link_budget.AtmosphereMask    = truth.r_mask/1000 - EARTH_RADIUS; % km
link_budget.NoiseTemp         = truth.Ts;                         % K
link_budget.AtmAttenuation    = truth.Ae;                         % dB
link_budget.TransPowerLevel   = truth.sv_power;                   % enum
link_budget.TransPowerOffset  = truth.xmit_power_offset;          % truth units? gpsmeas: dB
link_budget.GPSBlock          = 6;                   % use a GPS 2D antenna with L1, see gpsmeas
link_budget.TXAntennaMask     = truth.xmit_ant_mask;              % rad
link_budget.ReceiverNoise     = truth.Nf;                         % dB
link_budget.RecConversionLoss = truth.L;                          % dB
link_budget.SystemLoss        = truth.As;                         % dB
link_budget.LNAGain           = truth.Ga;                         % dB
link_budget.CableLoss         = truth.Ac;                         % dB
link_budget.RecAcqThresh      = truth.CN0_lim;                    % dB-Hz
link_budget.RecTrackThresh    = truth.CN0_lim;                    % dB-Hz
link_budget.DynamicTrackRange = truth.dyn_range;                  % dB
link_budget.TX_AntennaPointing= -1; % 1 for zenith pointing, -1 for nadir pointing
measOptions = setOdtbxOptions(measOptions, 'linkbudget', link_budget);
    
%% Call gpsmeas.
% Note we are using the qatt as well as the AntennaOrientation with our 2D
% antenna models.
fprintf(1,'...Calling gpsmeas to generate data, please wait...\n');
[y,~,~,AntLB] = gpsmeas(t, x, measOptions, qatt);

%% Fabricate GPS receiver inputs: 
% gather gpsmeas outputs into a gps data struct
gps_meas = makeGpsData();

prnind = 0;

for i = 1:32

    % only use the data that is visible:
    visind = find(AntLB{1}.HCN0(i,:) > truth.CN0_lim & ~isnan(y(1,:)) & ~isnan(y(2,:))); % visible indices
    
    if ~isempty(visind)
        prnind = prnind + 1;

        times_in_gps = (convertTime('GPS','UTC',epoch) + t(visind)/86400);
        gps_meas.PRN_data{prnind}.epoch          = times_in_gps;
        gps_meas.PRN_data{prnind}.raw_SNR        = AntLB{1}.HCN0(i,visind);
        % Note, the gpsmeas output range is not truly pseudorange (there are no
        % clock errors) but it is representative enough for this test. Pseudorange
        % is not used in the calculations that follow.
        gps_meas.PRN_data{prnind}.pseudorange    = y(1, visind); % range, km
        gps_meas.PRN_data{prnind}.doppler        = y(2, visind); % doppler
        gps_meas.PRN_data{prnind}.phase          = zeros(1,length(visind));
        
        gps_meas.GPS_PRN(prnind) = i;
    end
end
