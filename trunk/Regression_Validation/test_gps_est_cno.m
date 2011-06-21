function fail = test_gps_est_cno()
%
% Regression and unit test for gps_est_cno.m
% Note, the heavy lifting of gpslinkbuget is handled by that test.
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
d2r = pi/180; % degrees to radians
r2d = 180/pi; % degrees to radians
fail = 0;
GAIN_TOL = 1e-12; % 2D interp from gpslinkbudget, two patterns

%% TEST 1 - nominal test
RX_link.Nf = 0;
RX_link.L = 0;
RX_link.freq = 1575.42e6;    % GPS L1 frequency, Hz
RX_link.Ts = 300; % K
RX_link.As = 0;
RX_link.Ae = 0;

TX_link.P_sv = 1; % reflects the truth.sv_power, xmit_power_offset, and block in gpsmeas

real_2d = zeros(10,10);
real_2d(1,2:end) = ([0 5 10 30 60 90 120 150 180]*2)'; %az
real_2d(2:end,1) = [0 5 10 30 60 90 120 150 180]; % el
for i = 2:10
    for j = 2:10
        real_2d(i,j) = real_2d(1,i) + real_2d(j,1);
    end
end

rx_el = (0:10:180)*d2r;
rx_az = (0:20:360)*d2r;
tx_el = rx_el; % totally fake, but numerically useful
tx_az = rx_az; % totally fake, but numerically useful
los_mag = ones(1,19)*20000; % static values
epoch = datenum('19-May-2000');
t = epoch + (1:length(rx_el))/86400;

tx_id = makeGpsTXID(234, 31, 1);

phys_param = struct('epoch',t,...
    'TX_az',tx_az,'TX_el',tx_el,'RX_az',rx_az,'RX_el',rx_el,...
    'range',los_mag,...
    'range_rate',los_mag,... % bogus, but not a factor in gps_gain.m
    'GPS_yaw',zeros(1,length(rx_el)),...
    'meta',[]);
phys_param.meta = struct('RX_meta',[],...
    'RX_state_source',-99,...
    'TX_ID',tx_id,'TX_state_source',[],'gen_date',now);
phys_param.meta.RX_meta.RX_ID              = -99; %user-defined reveiver system identifier
phys_param.meta.RX_meta.meas_file          = '_synthetic_no_file';
phys_param.meta.RX_meta.obs_metadata{1}    = 'Synthetic measurement data from test_gps_est_cno.m';

CN0 = gps_est_cno(phys_param, RX_link, TX_link, real_2d, real_2d); % no filters

% test the values
if any(abs(CN0(1,:)-t) > 1e-15)
    fprintf(1,'Failed CN0 time check, test 1.\n');
    fail = 1;
end


if any(abs(CN0(2,:)-CN0(2,1)-(rx_el*r2d)-(rx_az*r2d)-(tx_el*r2d)-(tx_az*r2d)) > GAIN_TOL)
    fprintf(1,'Failed tx_gain gain check, test 1.\n');
    fail = 1;
end

% link budget check for baseline value
C = JATConstant('c') / 1000;  % km/s Speed of light
Ad_test = 20.*log10((C/RX_link.freq)./(4*pi.*los_mag(1)));
AP_test = TX_link.P_sv + Ad_test + RX_link.Ae;
RP_test = AP_test + RX_link.As;
CN0_test = RP_test - (10*log10(RX_link.Ts)) + 228.6 + RX_link.Nf + RX_link.L;

if any(abs(CN0(2,1)-CN0_test) > GAIN_TOL)
    fprintf(1,'Failed rx_gain CN0 gain check, test 1.\n');
    fail = 1;
end

%% TEST 2 - missing arg checks

% arguments: gps_meas, phys_param, RX_link, TX_link, pattern, ant)
% Missing args or bad ant value
arglist = [...
    0       1           1       1       1;...    
    1       0           1       1       1;...
    1       1           0       1       1;...
    1       1           1       0       1;...
    1       1           1       1       0;...
    ];
for i = 1:size(arglist,1)
    if arglist(i,1); arg1=phys_param; else arg1=[];end
    if arglist(i,2); arg2=RX_link; else arg2=[];end
    if arglist(i,3); arg3=TX_link; else arg3=[];end
    if arglist(i,4); arg4=real_2d; else arg4=[];end
    if arglist(i,5); arg5=real_2d; else arg5=[];end
    try
        CN0 = gps_est_cno(arg1, arg2, arg3, arg4, arg5); % no filters
        fprintf(1,'Failed argument check, case %d, test 2.\n',i);
        fail = 1;
    catch ex %#ok<NASGU>
        % expected
    end
end

%% TEST 3 - bad arg: 1D antenna pattern
% a 1D omnidirectional antenna with no gain
omni_1d = zeros(9,2);
omni_1d(:,1) = [0 5 10 30 60 90 120 150 180]';
try
    CN0 = gps_est_cno(phys_param, RX_link, TX_link, omni_1d, real_2d); %#ok<NASGU> % no filters
    fprintf(1,'Failed 1D pattern check, test 3a.\n');
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

try
    CN0 = gps_est_cno(phys_param, RX_link, TX_link, real_2d, omni_1d); %#ok<NASGU> % no filters
    fprintf(1,'Failed 1D pattern check, test 3b.\n');
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

%% TEST 4 - filter test 1
RX_link.Nf = 0;
RX_link.L = 0;
RX_link.freq = 1575.42e6;    % GPS L1 frequency, Hz
RX_link.Ts = 300; % K
RX_link.As = 0;
RX_link.Ae = 0;

TX_link.P_sv = 1; % reflects the truth.sv_power, xmit_power_offset, and block in gpsmeas

real_2d = zeros(10,10);
real_2d(1,2:end) = ([0 5 10 30 60 90 120 150 180]*2)'; %az
real_2d(2:end,1) = [0 5 10 30 60 90 120 150 180]; % el
for i = 2:10
    for j = 2:10
        real_2d(i,j) = real_2d(1,i) + real_2d(j,1);
    end
end

rx_el = (0:10:180)*d2r;
rx_az = (0:20:360)*d2r;
tx_el = rx_el; % totally fake, but numerically useful
tx_az = rx_az; % totally fake, but numerically useful
los_mag = ones(1,19)*20000; % static values
epoch = datenum('19-May-2000');
t = epoch + (1:length(rx_el))/86400;

tx_id = makeGpsTXID(234, 31, 1);

phys_param = struct('epoch',t,...
    'TX_az',tx_az,'TX_el',tx_el,'RX_az',rx_az,'RX_el',rx_el,...
    'range',los_mag,...
    'range_rate',los_mag,... % bogus, but not a factor in gps_gain.m
    'GPS_yaw',zeros(1,length(rx_el)),...
    'meta',[]);
phys_param.meta = struct('RX_meta',[],...
    'RX_state_source',-99,...
    'TX_ID',tx_id,'TX_state_source',[],'gen_date',now);
phys_param.meta.RX_meta.RX_ID              = -99; %user-defined reveiver system identifier
phys_param.meta.RX_meta.meas_file          = '_synthetic_no_file';
phys_param.meta.RX_meta.obs_metadata{1}    = 'Synthetic measurement data from test_gps_est_cno.m';

filter = FilterPpEl(45*d2r,'tx');

resind = find(tx_el <= 45*d2r);

[CN0, fil_phys_param] = gps_est_cno(phys_param, RX_link, TX_link, real_2d, real_2d, filter); 

% check sizes
if isempty(CN0) || any(size(CN0) ~= [2 length(resind)])
    fprintf(1,'Failed CN0 size check, test 4.\n');
end

if isempty(fil_phys_param) || any(size(fil_phys_param.epoch) ~= [1 length(resind)])
    fprintf(1,'Failed filtered phys param size check, test 4.\n');
    fail = 1;
end

% test the values
if any(abs(CN0(1,:)-t(resind)) > 1e-15)
    fprintf(1,'Failed CN0 time check, test 4.\n');
    fail = 1;
end


if any(abs(CN0(2,:)-CN0(2,1)-(rx_el(resind)*r2d)-(rx_az(resind)*r2d)-(tx_el(resind)*r2d)-(tx_az(resind)*r2d)) > GAIN_TOL)
    fprintf(1,'Failed tx_gain gain check, test 4.\n');
    fail = 1;
end

% link budget check for baseline value
C = JATConstant('c') / 1000;  % km/s Speed of light
Ad_test = 20.*log10((C/RX_link.freq)./(4*pi.*los_mag(1)));
AP_test = TX_link.P_sv + Ad_test + RX_link.Ae;
RP_test = AP_test + RX_link.As;
CN0_test = RP_test - (10*log10(RX_link.Ts)) + 228.6 + RX_link.Nf + RX_link.L;

if any(abs(CN0(2,1)-CN0_test) > GAIN_TOL)
    fprintf(1,'Failed rx_gain CN0 gain check, test 4.\n');
    fail = 1;
end

% check filtered physical param values
if any(fil_phys_param.epoch ~= phys_param.epoch(resind))
    fprintf(1,'Failed filtered phys param data check, test 4.\n');
    fail = 1;
end

%% TEST 5 - filter test 2
RX_link.Nf = 0;
RX_link.L = 0;
RX_link.freq = 1575.42e6;    % GPS L1 frequency, Hz
RX_link.Ts = 300; % K
RX_link.As = 0;
RX_link.Ae = 0;

TX_link.P_sv = 1; % reflects the truth.sv_power, xmit_power_offset, and block in gpsmeas

real_2d = zeros(10,10);
real_2d(1,2:end) = ([0 5 10 30 60 90 120 150 180]*2)'; %az
real_2d(2:end,1) = [0 5 10 30 60 90 120 150 180]; % el
for i = 2:10
    for j = 2:10
        real_2d(i,j) = real_2d(1,i) + real_2d(j,1);
    end
end

rx_el = (0:10:180)*d2r;
rx_az = (0:20:360)*d2r;
tx_el = rx_el; % totally fake, but numerically useful
tx_az = rx_az; % totally fake, but numerically useful
los_mag = ones(1,19)*20000; % static values
epoch = datenum('19-May-2000');
t = epoch + (1:length(rx_el))/86400;

tx_id = makeGpsTXID(234, 31, 1);

phys_param = struct('epoch',t,...
    'TX_az',tx_az,'TX_el',tx_el,'RX_az',rx_az,'RX_el',rx_el,...
    'range',los_mag,...
    'range_rate',los_mag,... % bogus, but not a factor in gps_gain.m
    'GPS_yaw',zeros(1,length(rx_el)),...
    'meta',[]);
phys_param.meta = struct('RX_meta',[],...
    'RX_state_source',-99,...
    'TX_ID',tx_id,'TX_state_source',[],'gen_date',now);
phys_param.meta.RX_meta.RX_ID              = -99; %user-defined reveiver system identifier
phys_param.meta.RX_meta.meas_file          = '_synthetic_no_file';
phys_param.meta.RX_meta.obs_metadata{1}    = 'Synthetic measurement data from test_gps_est_cno.m';

filter = cell(1);
filter{1} = FilterPpEl(45*d2r,'tx');
filter{2} = FilterPpEl(35*d2r,'rx');

resind = find(rx_el <= 35*d2r);

[CN0, fil_phys_param] = gps_est_cno(phys_param, RX_link, TX_link, real_2d, real_2d, filter); 

% check sizes
if isempty(CN0) || any(size(CN0) ~= [2 length(resind)])
    fprintf(1,'Failed CN0 size check, test 5.\n');
end

if isempty(fil_phys_param) || any(size(fil_phys_param.epoch) ~= [1 length(resind)])
    fprintf(1,'Failed filtered phys param size check, test 5.\n');
    fail = 1;
end

% test the values
if any(abs(CN0(1,:)-t(resind)) > 1e-15)
    fprintf(1,'Failed CN0 time check, test 5.\n');
    fail = 1;
end


if any(abs(CN0(2,:)-CN0(2,1)-(rx_el(resind)*r2d)-(rx_az(resind)*r2d)-(tx_el(resind)*r2d)-(tx_az(resind)*r2d)) > GAIN_TOL)
    fprintf(1,'Failed tx_gain gain check, test 5.\n');
    fail = 1;
end

% link budget check for baseline value
C = JATConstant('c') / 1000;  % km/s Speed of light
Ad_test = 20.*log10((C/RX_link.freq)./(4*pi.*los_mag(1)));
AP_test = TX_link.P_sv + Ad_test + RX_link.Ae;
RP_test = AP_test + RX_link.As;
CN0_test = RP_test - (10*log10(RX_link.Ts)) + 228.6 + RX_link.Nf + RX_link.L;

if any(abs(CN0(2,1)-CN0_test) > GAIN_TOL)
    fprintf(1,'Failed rx_gain CN0 gain check, test 5.\n');
    fail = 1;
end

% check filtered physical param values
if any(fil_phys_param.epoch ~= phys_param.epoch(resind))
    fprintf(1,'Failed filtered phys param data check, test 5.\n');
    fail = 1;
end
