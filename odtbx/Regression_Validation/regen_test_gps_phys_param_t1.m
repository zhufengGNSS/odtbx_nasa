% A script to regenerate the t1 data used in test_gps_phys_params.m.
% This uses gpsmeas and one of its regression test data files to synthesize
% test data.  Results are stored in the file: test_gps_phys_param_dat1_NEW.
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

% constants
r2d=180/pi; % rad to deg
d2r=pi/180; % deg to rad

%% input data
truth_file = 'RESULTS_2D_user_and_GPS_antenna.mat'; % both 3D antennas
yuma_file = 'Yuma1134.txt';

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

measOptions.('epoch')             = epoch;
measOptions.('useRange')          = true;
measOptions.('useRangeRate')      = false;
measOptions.('useDoppler')        = true;
%rSigma (use default, no analogue in gpsdef)
measOptions.('GPSBand')           = 'L1';
measOptions.('YumaFile')          = yuma_file;
%Rotation2ECI (use default, no analogue in gpsdef)
measOptions.('AntennaPointing')   = ant_point;
measOptions.('AntennaPattern')    = eval(ant_rx_pat);
measOptions.('AntennaMask')       = truth.rcv_ant_mask;
measOptions.('AtmosphereMask')    = truth.r_mask/1000 - EARTH_RADIUS; % km
measOptions.('useiono')           = false;
measOptions.('usetropo')          = false;
measOptions.('NoiseTemp')         = truth.Ts;                         % K
measOptions.('AtmAttenuation')    = truth.Ae;                         % dB
measOptions.('TransPowerLevel')   = truth.sv_power;                   % enum
measOptions.('TransPowerOffset')  = truth.xmit_power_offset;          % truth units? gpsmeas: dB
measOptions.('GPSBlock')          = 6;                                % use a GPS 2D antenna with L1, see gpsmeas
measOptions.('TransAntMask')      = truth.xmit_ant_mask;              % rad
measOptions.('ReceiverNoise')     = truth.Nf;                         % dB
measOptions.('RecConversionLoss') = truth.L;                          % dB
measOptions.('SystemLoss')        = truth.As;                         % dB
measOptions.('LNAGain')           = truth.Ga;                         % dB
measOptions.('CableLoss')         = truth.Ac;                         % dB
measOptions.('RecAcqThresh')      = truth.CN0_lim;                    % dB-Hz
measOptions.('RecTrackThresh')    = truth.CN0_lim;                    % dB-Hz
measOptions.('DynamicTrackRange') = truth.dyn_range;                  % dB
measOptions.('PrecnNutnExpire')   = 0.1;                                % days
measOptions.('AntennaOrientation') = dcm('ax2',pi/2);

%% Call gpsmeas.
% Note we are using the qatt as well as the AntennaOrientation with our 2D
% antenna models.
fprintf(1,'...Calling gpsmeas...\n');
[y,H,~,AntLB] = gpsmeas(t, x, measOptions, qatt);

%% Fabricate GPS receiver inputs: 
% gather gpsmeas outputs into a gps data struct
gps_meas = makeGpsData();

% use satellite #1
gps_meas.GPS_PRN = 1;

% only use the data that is visible:
visind = find(AntLB{1}.HCN0(1,:) > truth.CN0_lim); % visible indices
times_in_gps = (convertTime('GPS','UTC',epoch) + t(visind)/86400);
meas_time_bias = 0; % sec
gps_meas.PRN_data{1}.epoch          = times_in_gps + meas_time_bias/86400; % bias meas times by 0s
gps_meas.PRN_data{1}.raw_SNR        = AntLB{1}.HCN0(gps_meas.GPS_PRN,visind);
% Note, the gpsmeas output range is not truly pseudorange (there are no
% clock errors) but it is representative enough for this test. Pseudorange
% is not used in the calculations that follow.
gps_meas.PRN_data{1}.pseudorange    = y(1, visind); % range, km
gps_meas.PRN_data{1}.doppler        = y(2, visind); % doppler
gps_meas.PRN_data{1}.phase          = zeros(1,length(visind));

% note the metadata for traceability
gps_meas.RX_meta.RX_ID              = -99; %user-defined reveiver system identifier
gps_meas.RX_meta.meas_file          = '_synthetic_no_file';
gps_meas.RX_meta.obs_metadata{1}    = 'Synthetic measurement data from gpsenh_vs_gpsmeas.m';

% note the state source
RX_state_source = 'gpsenh_vs_gpsmeas testing';

% let's identify this GPS SV transmitter:
TX_ID{1} = makeGpsTXID(345, 1, 1);

%% call gps_phys_params
fprintf(1,'...Calling gps_phys_params...\n');
phys_param = gps_phys_params(gps_meas, RX_state_source, -99, epoch, ...
    t, x, qatt, [], getOdtbxOptions(measOptions,'AntennaOrientation'), ...
    yuma_file, TX_ID, []);

%% save the regression data
save test_gps_phys_param_dat1_NEW epoch gps_meas qatt t x phys_param
