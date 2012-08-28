function fail = gps_phys_params_test()
%
% Regression and unit test for gps_phys_params.m
% Note, the heavy lifting performed by getgpsmeas is not directly checked 
% by this test.  However comparison data is used to detect regressions.
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

fail = 0;

% Tolerance history:
%
% Initially, R2010b and R2009b on Win platform had passing at 1e-12.
% 
% Bumped TOL from 1e-12, which should have passed to 1e-11 due to
% repeatability bug on Linux x86-64 platforms discussed in Mantis issue 360:
% https://gs-fftb-collab.gsfc.nasa.gov/odtbx/mantisbt-1.2.3/view.php?id=360
%
% Bumped tolerance again for differences seen between MATLAB R2010b and
% R2011a, noted in Mantis issue 348:
% https://gs-fftb-collab.gsfc.nasa.gov/odtbx/mantisbt-1.2.3/view.php?id=348
% to 2e-11 (errors around 1.82e-11 in Test 1 .phys_param{1}.range check).
TOL = 2e-11; 

%% TEST 1 - nominal test
TESTNUM = 1;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

t1 = load('gps_phys_param_dat1.mat');
% contents:
%   t1.epoch                1x1                  8  double              
%   t1.gps_meas             1x1              12918  struct              
%   t1.qatt                 4x841            26912  double              
%   t1.t                    1x841             6728  double              
%   t1.x                    6x841            40368  double  
%   t1.phys_param: {[1x1 struct]}
% See regen_test_gps_phys_param_t1.m to recreate this test data.

RX_ID = -99; 
RX_state_source = 'gps_phys_params_test.m';
Rotation2ECI = [];
AntOrient = dcm('ax2',pi/2);
yuma_file = 'Yuma1134.txt';

% let's identify the GPS SV transmitters:
TX_ID{1} = makeGpsTXID(341, 1, 1); % PRN 1
TX_ID{2} = makeGpsTXID(342, 2, 1); % PRN 2
TX_ID{3} = makeGpsTXID(344, 4, 2); % PRN 4, different block type

phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
    t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
    yuma_file, TX_ID, []);

% basic output checks
if isempty(phys_param)
    fprintf(1,'Failed test %d: unexpected empty output.\n',TESTNUM);
    fail = 1;
end
if ~iscell(phys_param)
    fprintf(1,'Failed test %d: unexpected output type.\n',TESTNUM);
    fail = 1;
end    
if length(phys_param) ~= 1
    fprintf(1,'Failed test %d: unexpected output size.\n',TESTNUM);
    fail = 1;
end

% basic struct checks
sz = 282;
if ~isfield(phys_param{1},'epoch') || (length(phys_param{1}.epoch) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'epoch');
    fail = 1;
end
if ~isfield(phys_param{1},'TX_az') || (length(phys_param{1}.TX_az) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'TX_az');
    fail = 1;
end
if ~isfield(phys_param{1},'TX_el') || (length(phys_param{1}.TX_el) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'TX_el');
    fail = 1;
end
if ~isfield(phys_param{1},'RX_az') || (length(phys_param{1}.RX_az) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'RX_az');
    fail = 1;
end
if ~isfield(phys_param{1},'RX_el') || (length(phys_param{1}.RX_el) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'RX_el');
    fail = 1;
end
if ~isfield(phys_param{1},'range') || (length(phys_param{1}.range) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'range');
    fail = 1;
end
if ~isfield(phys_param{1},'range_rate') || (length(phys_param{1}.range_rate) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'range_rate');
    fail = 1;
end
if ~isfield(phys_param{1},'GPS_yaw') || (length(phys_param{1}.GPS_yaw) ~= sz)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'GPS_yaw');
    fail = 1;
end
if ~isfield(phys_param{1},'meta') || (length(phys_param{1}.meta) ~= 1)
    fprintf(1,'Failed test %d: unexpected output field: %s.\n',TESTNUM,'meta');
    fail = 1;
end

% some quick metadata checks
if phys_param{1}.meta.RX_meta.RX_ID ~= RX_ID
    fprintf(1,'Failed test %d: metadata failure: %s.\n',TESTNUM,'phys_param{1}.meta.RX_meta');
    fail = 1;
end
if phys_param{1}.meta.TX_ID.GPS_ID ~= TX_ID{1}.GPS_ID
    fprintf(1,'Failed test %d: metadata failure: %s.\n',TESTNUM,'phys_param{1}.meta.TX_ID.GPS_ID');
    fail = 1;
end
if phys_param{1}.meta.TX_ID.GPS_PRN ~= TX_ID{1}.GPS_PRN
    fprintf(1,'Failed test %d: metadata failure: %s.\n',TESTNUM,'phys_param{1}.meta.TX_ID.GPS_PRN');
    fail = 1;
end
if phys_param{1}.meta.TX_ID.block_type ~= TX_ID{1}.block_type
    fprintf(1,'Failed test %d: metadata failure: %s.\n',TESTNUM,'phys_param{1}.meta.TX_ID.block_type');
    fail = 1;
end

% regression check against original data
if any(abs(phys_param{1}.epoch - t1.phys_param{1}.epoch) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.epoch');
    fail = 1;
end
if any(abs(phys_param{1}.TX_az - t1.phys_param{1}.TX_az) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.TX_az');
    fail = 1;
end
if any(abs(phys_param{1}.TX_el - t1.phys_param{1}.TX_el) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.TX_el');
    fail = 1;
end
if any(abs(phys_param{1}.RX_az - t1.phys_param{1}.RX_az) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.RX_az');
    fail = 1;
end
if any(abs(phys_param{1}.RX_el - t1.phys_param{1}.RX_el) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.RX_el');
    fail = 1;
end
if any(abs(phys_param{1}.range - t1.phys_param{1}.range) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.range');
    fail = 1;
end
if any(abs(phys_param{1}.range_rate - t1.phys_param{1}.range_rate) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.range_rate');
    fail = 1;
end
if any(abs(phys_param{1}.GPS_yaw - t1.phys_param{1}.GPS_yaw) > TOL)
    fprintf(1,'Failed test %d: regression data failure: %s.\n',TESTNUM,'phys_param{1}.GPS_yaw');
    fail = 1;
end


%% TEST 2 - missing arg checks
TESTNUM = 2;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

%arguments:
% gps_meas, RX_state_source, RX_ID, epoch, tsc, xsc, qatt, Rotation2ECI, ant_bod, gps_alm_file, TX_IDs, filters)
% Note: qatt, Rotation2ECI, & ant_bod can be empty.
arglist = [...
    0           1           1       1       1   1                                   1               1;...
    1           0           1       1       1   1                                   1               1;...
    1           1           0       1       1   1                                   1               1;...
    1           1           1       0       1   1                                   1               1;...
    1           1           1       1       0   1                                   1               1;...
    1           1           1       1       1   0                                   1               1;...
    1           1           1       1       1   1                                   0               1;...
    1           1           1       1       1   1                                   1               0;...
    ];
for i = 1:size(arglist,1)
    if arglist(i,1); arg1=t1.gps_meas; else arg1=[];end
    if arglist(i,2); arg2=RX_state_source; else arg2=[];end
    if arglist(i,3); arg3=RX_ID; else arg3=[];end
    if arglist(i,4); arg4=t1.epoch; else arg4=[];end
    if arglist(i,5); arg5=t1.t; else arg5=[];end
    if arglist(i,6); arg6=t1.x; else arg6=[];end
    arg7=t1.qatt;
    arg8=Rotation2ECI;
    arg9=AntOrient;
    if arglist(i,7); arg10=yuma_file; else arg10=[];end
    if arglist(i,8); arg11=TX_ID; else arg11=[];end
    try
        [phys_param] = gps_phys_params(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, []); % no filters
        fprintf(1,'Failed test %d: argument check, case %d.\n',TESTNUM,i);
        fail = 1;
    catch ex 
        % expected
    end
end


%% TEST 3 - bad & mismatched arg checks
TESTNUM = 3;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

% bad time size
try
    phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
        zeros(3,1), zeros(6,1), t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad time size.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad state size
try
    phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, zeros(7,1), t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad state size.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% time vs state length mismatch
try
    phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x(:,1:end-1), t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: time vs state length mismatch.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad qatt size vs time & state
try
    phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt(:,1:end-1), Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad qatt size vs time & state.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad Rotation2ECI function signature
try
    phys_param = gps_phys_params(t1.gps_meas, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, @BadRotFcn, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad Rotation2ECI function signature.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad gpsmeas struct, 1
try
    phys_param = gps_phys_params(struct('x',[]), RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad gpsmeas struct, 1.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad gpsmeas struct, 2
% no meta
gmfail2 = struct('GPS_PRN',t1.gps_meas.GPS_PRN,'PRN_data',t1.gps_meas.PRN_data);
try
    phys_param = gps_phys_params(gmfail2, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad gpsmeas struct, 2.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad gpsmeas struct, 3
% no RX_ID
gmfail3 = gmfail2;
gmfail3.meta = struct('a',[]); % bogus
try
    phys_param = gps_phys_params(gmfail3, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad gpsmeas struct, 3.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% bad RX_ID cross-check
gmfail4 = t1.gps_meas;
gmfail4.RX_meta.RX_ID = 987;
try
    phys_param = gps_phys_params(gmfail4, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: bad RX_ID cross-check.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

% duplicate GPS PRNs
gmfail5 = t1.gps_meas;
gmfail5.GPS_PRN = [gmfail5.GPS_PRN gmfail5.GPS_PRN]; % dupe
gmfail5.PRN_data{2} = gmfail5.PRN_data{1};
try
    phys_param = gps_phys_params(gmfail5, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    fprintf(1,'Failed test %d: duplicate GPS PRNs.\n',TESTNUM);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end


%% TEST 4 - off-nominal cases: gps_meas data PRNs
TESTNUM = 4;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

% valid, but empty gps_meas data
gmempty = t1.gps_meas;
gmempty.GPS_PRN = [];
gmempty.PRN_data{1}.epoch = [];
gmempty.PRN_data{1}.raw_SNR = [];
gmempty.PRN_data{1}.pseudorange = [];
gmempty.PRN_data{1}.doppler = [];
gmempty.PRN_data{1}.phase = [];
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmempty, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: valid, but empty gps_meas data - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: valid, but empty gps_meas data - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

% all unmatched GPS PRNs in gps_meas data
gmunmatched = t1.gps_meas;
gmunmatched.GPS_PRN = 13;
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmunmatched, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: all unmatched GPS PRNs - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: all unmatched GPS PRNs - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

% unsupported PRN in almanac
gmunmatched.GPS_PRN = 32;
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmunmatched, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: all unmatched GPS PRNs - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: all unmatched GPS PRNs - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

% unhealthy PRN in almanac
gmunmatched.GPS_PRN = 29; % unhealthy PRN 29 in yuma_file=Yuma1134_unhealthy29.txt
unhealthy_yuma_file = 'Yuma1134_unhealthy29.txt';
% let's identify the GPS SV transmitters:
TX_ID_unh{1} = makeGpsTXID(351, 1, 1); % PRN 1
TX_ID_unh{2} = makeGpsTXID(352, 2, 1); % PRN 2
TX_ID_unh{3} = makeGpsTXID(354, 29, 1); % PRN 29
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmunmatched, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        unhealthy_yuma_file, TX_ID_unh, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: all unmatched GPS PRNs - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: all unmatched GPS PRNs - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end


%% TEST 5 - multiple PRNS
TESTNUM = 5;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

% multiple I/O: multiple phys_meas cell array structs
gmmult = t1.gps_meas;
gmmult.GPS_PRN = [1 2 4];
gmmult.PRN_data{2} = gmmult.PRN_data{1};
gmmult.PRN_data{3} = gmmult.PRN_data{1};
try
    phys_param = gps_phys_params(gmmult, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if isempty(phys_param)
        fprintf(1,'Failed test %d: multiple I/O - unexpectedly empty output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: multiple I/O - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end
if any(size(phys_param) ~= [1 3])
    fprintf(1,'Failed test %d: multiple I/O - wrong output data size.\n',TESTNUM);
    fail = 1;
end

%% TEST 6 - filtering multiple PRNS, multiple results
TESTNUM = 6;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

% multiple I/O: multiple phys_meas cell array structs
gmmult = t1.gps_meas;
gmmult.GPS_PRN = [1 2 4];
gmmult.PRN_data{2} = gmmult.PRN_data{1};
gmmult.PRN_data{3} = gmmult.PRN_data{1};

% filter on GPS SV block, only PRNs 1 & 2 are block "1" (II/IIA)
filter = FilterGpsBlock('include',TX_ID, 1);

try
    [phys_param, fgm] = gps_phys_params(gmmult, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, filter);
    if isempty(phys_param)
        fprintf(1,'Failed test %d: multiple I/O - unexpectedly empty output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: multiple I/O - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end
if any(size(phys_param) ~= [1 2])
    fprintf(1,'Failed test %d: multiple I/O - wrong output data size.\n',TESTNUM);
    fail = 1;
end
if phys_param{1}.meta.TX_ID.GPS_PRN ~= 1
    fprintf(1,'Failed test %d: multiple I/O - unexpected PRN in slot 1.\n',TESTNUM);
    fail = 1;
end
if phys_param{2}.meta.TX_ID.GPS_PRN ~= 2
    fprintf(1,'Failed test %d: multiple I/O - unexpected PRN in slot 2.\n',TESTNUM);
    fail = 1;
end
% check the filter gps meas data
if any(fgm.GPS_PRN ~= [1 2])
    fprintf(1,'Failed test %d: multiple I/O - unexpected PRN results in filtered gps data.\n',TESTNUM);
    fail = 1;
end

%% TEST 7 - off-nominal cases: gps_meas times
TESTNUM = 7;
fprintf(1,'\nExecuting gps_phys_params test %d:\n\n',TESTNUM);

% all meas times before states
gmtime = t1.gps_meas;
ebias = gmtime.PRN_data{1}.epoch(end) - convertTime('GPS','UTC',t1.epoch);
gmtime.PRN_data{1}.epoch = gmtime.PRN_data{1}.epoch - ebias - 1/86400; % and one more second so no times match
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmtime, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: all meas times before states - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: all meas times before states - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

% all meas times after states
gmtime = t1.gps_meas;
ebias = gmtime.PRN_data{1}.epoch(1) - ...
    convertTime('GPS','UTC',(t1.epoch + ((t1.t(end)-t1.t(1) + 1)/86400)) ); % and one more second so no times match
gmtime.PRN_data{1}.epoch = gmtime.PRN_data{1}.epoch - ebias;
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmtime, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if ~isempty(phys_param)
        fprintf(1,'Failed test %d: all meas times after states - unexpected output data.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: all meas times after states - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

% meas times before and after state times
gmtime = t1.gps_meas;
ebias = gmtime.PRN_data{1}.epoch(10) - gmtime.PRN_data{1}.epoch(1);
gmtime.PRN_data{1}.epoch(1:10) = gmtime.PRN_data{1}.epoch(1:10) - ebias;
ebias = gmtime.PRN_data{1}.epoch(end) - gmtime.PRN_data{1}.epoch(end-10);
gmtime.PRN_data{1}.epoch(end-10:end) = gmtime.PRN_data{1}.epoch(end-10:end) + ebias;
try
    fprintf(1,'\n\ngps_phys_params_test.m: The following warning is expected:\n\n');
    phys_param = gps_phys_params(gmtime, RX_state_source, RX_ID, t1.epoch, ...
        t1.t, t1.x, t1.qatt, Rotation2ECI, AntOrient, ...
        yuma_file, TX_ID, []);
    if isempty(phys_param)
        fprintf(1,'Failed test %d: meas times before and after state times - unexpectedly empty output data.\n',TESTNUM);
        fail = 1;
    end
    if length(phys_param) ~= 1
        fprintf(1,'Failed test %d: meas times before and after state times - incorrect output cell array length.\n',TESTNUM);
        fail = 1;
    end
    if length(phys_param{1}.epoch) ~= (length(gmtime.PRN_data{1}.epoch)-19)
        fprintf(1,'Failed test %d: meas times before and after state times - incorrect output data length.\n',TESTNUM);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: meas times before and after state times - caught unexpected exception.\n',TESTNUM);
    fail = 1;
end

end % function

%% A "bad" function for rotation that doesn't take a time argument.
function I = BadRotFcn()
    I = eye(6);
end