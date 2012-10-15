function fail = gpsenh_vs_gpsmeas_test()
% An end-to-end comparison test between gpsmeas and the GPS Enhancement
% Toolset tools: gps_phys_param, gps_est_cno, and gps_gain.
%
% Note: uses one of the datasets from the gpsmeas regression test.
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
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur             08/27/2012      Rename to conform to new
%                                           regression test format 

fprintf(1,'Starting gpsmeas vs GPS Enhancement tools test...\n');
fail = 0;

% constants
% r2d=180/pi; % rad to deg
% d2r=pi/180; % deg to rad

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
[y,~,~,AntLB] = gpsmeas(t, x, measOptions, qatt);

%% Fabricate GPS receiver inputs: 
% gather gpsmeas outputs into a gps data struct
gps_meas = makeGpsData();

% use satellite #1
gps_meas.GPS_PRN = 1;

% only use the data that is visible:
visind = find(AntLB{1}.HCN0(1,:) > truth.CN0_lim); % visible indices
times_in_gps = (convertTime('GPS','UTC',epoch) + t(visind)/86400);
meas_time_bias = 0; % sec (if you bias the times then the exact comparisons below won't work but the plots will)
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

%% the antenna patterns
rx_pat = load(truth.rec_pattern);
xmit_pattern = 'GPSIIA_L1_3Dantenna.txt';	 % from gpsmeas, with above inputs
tx_pat = load(xmit_pattern);

%% Calculate the estimated CN0 from these physical parameters:
RX_link.Nf = truth.Nf;
RX_link.L = truth.L;
RX_link.freq = 1575.42e6;    % GPS L1 frequency, Hz
RX_link.Ts = truth.Ts;
RX_link.As = truth.As;
RX_link.Ae = truth.Ae;
TX_link.P_sv = 14.9; % reflects the truth.sv_power, xmit_power_offset, and block in gpsmeas

%% Check some of the physical params
% some measurements may have been discarded, check via epoch first to find
% the indices into the gps data that match the physical params data
ppe_gps = convertTime('GPS','UTC',phys_param{1}.epoch);
f = FilterGpsExplicitTimes(ppe_gps, 0.1/86400); % 0.1 sec tolerance
gpsdatind = f.getInd(gps_meas, 1); % gps indices of matching phys param data for 1st PRN

if length(gpsdatind) ~= length(phys_param{1}.epoch)
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: time index mismatch!');
end

rangediff = phys_param{1}.range-gps_meas.PRN_data{1}.pseudorange(gpsdatind);
if any(abs(rangediff) > 1e-4) %km
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: range from gps_phys_params fails to match gpsmeas pseudorange!');
else
    fprintf(1,'Passed physical param range comparison.\n');
end

c = JATConstant('c') / 1000; % km/sec speed of light
rrdiff = phys_param{1}.range_rate-(-gps_meas.PRN_data{1}.doppler(gpsdatind)/RX_link.freq*c);
if any(abs(rrdiff) > 1e-6) %km/s
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: range rate from gps_phys_params fails to match gpsmeas doppler!');
else
    fprintf(1,'Passed physical param range rate comparison.\n');
end

%% Call the CN0 estimator:
fprintf(1,'...Calling gps_est_cno...\n');
CN0_est = gps_est_cno(phys_param{1}, RX_link, TX_link, rx_pat, tx_pat);

%% check the resulting estimated CN0
% find the exact closest times gpsmeas used for comparison by simply
% removing the time bias
CN0diff = CN0_est(2,:)-AntLB{1}.HCN0(gps_meas.GPS_PRN,gpsdatind);

if abs(CN0diff) > 1e-8
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: estimated CN0 from gps_est_cno fails to match gpsmeas!');
else
    fprintf(1,'Passed estimated CN0 comparison.\n');
end

%% Gain calculations
% Now that we have the physical parameters, we have seen the estimated CN0,
% and we have fabricated receiver measurements, lets use these to back out
% 2D antenna gain values.
fprintf(1,'...Calling gps_gain...\n');
tx_gain = gps_gain(gps_meas, phys_param{1}, RX_link, TX_link, rx_pat, 1);
rx_gain = gps_gain(gps_meas, phys_param{1}, RX_link, TX_link, tx_pat, 2);

% %% plot the true receive rx_pattern
% hr = figure;surf(rx_pat(1,2:end),rx_pat(2:end,1),rx_pat(2:end,2:end),'LineStyle','none');
% xlabel('Az (deg)');
% ylabel('El (deg)');
% zlabel('Gain (dB)');
% title('Receive antenna gain');
% hold on;
% plot3(rx_gain(3,:)*r2d, rx_gain(4,:)*r2d, rx_gain(2,:),'.');
% hold off;
% view(-113,34);

%% Compare receive gain data
Ar_actual = interp2(rx_pat(1,2:end)*pi/180,rx_pat(2:end,1)*pi/180,rx_pat(2:end,2:end),rx_gain(3,:),rx_gain(4,:),'spline');
if max(abs(rx_gain(2,:)-Ar_actual)) > 1e-6 % dB, started at 1e-7
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: estimated receive gain from gps_gain fails to match actual pattern!\n');
else
    fprintf(1,'Passed estimated rx gain comparison.\n');
end


% %% plot the true transmit rx_pattern
% ht = figure;surf(tx_pat(1,2:end),tx_pat(2:end,1),tx_pat(2:end,2:end),'LineStyle','none');
% xlabel('Az (deg)');
% ylabel('El (deg)');
% zlabel('Gain (dB)');
% title('Transmit antenna gain');
% hold on;
% plot3(tx_gain(3,:)*r2d, tx_gain(4,:)*r2d, tx_gain(2,:),'.');
% hold off;
% view(-113,34);

%% Compare transmit gain data
At_actual = interp2(tx_pat(1,2:end)*pi/180,tx_pat(2:end,1)*pi/180,tx_pat(2:end,2:end),tx_gain(3,:),tx_gain(4,:),'spline');
if max(abs(tx_gain(2,:)-At_actual)) > 1e-6 % dB, started at 1e-7
    fail = 1;
    fprintf(1,'end2end_gpsenh_vs_gpsmeas: estimated transmit gain from gps_gain fails to match actual pattern!\n');
else
    fprintf(1,'Passed estimated tx gain comparison.\n');
end