function fail = test_FilterEstCnoThresh()
% Unit test case for the FilterEstCnoThresh class.
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

fail = 0;

d2r = pi/180;

%% TEST 1 - constructor
TESTNUM = 1;
fprintf(1,'\nExecuting test_FilterEstCnoThresh test %d:\n\n',TESTNUM);

CASE= 'no arg';
try
    f = FilterEstCnoThresh;
    if f.cno_thresh ~= 0
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

CASE= 'with arg 1';
try
    f = FilterEstCnoThresh(99);
    if f.cno_thresh ~= 99
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

%% TEST 2 - filter - empty arg
TESTNUM = 2;
fprintf(1,'\nExecuting test_FilterEstCnoThresh test %d:\n\n',TESTNUM);

f = FilterEstCnoThresh(22);

CASE= 'filter(empty arg)';
try
    res = f.filter([]);
    if ~isempty(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'filter(empty args)';
try
    [res1, res2, res3] = f.filter([]);
    if ~isempty(res1)
        fprintf(1,'Failed test %d: %s - bad value 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(res2)
        fprintf(1,'Failed test %d: %s - bad value 2.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(res3)
        fprintf(1,'Failed test %d: %s - bad value 3.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

%% TEST 3 - filter - empty phys params data struct
TESTNUM = 3;
fprintf(1,'\nExecuting test_FilterEstCnoThresh test %d:\n\n',TESTNUM);

% physical param data
dat = makePpData; % empty
% add metadata
dat.meta.RX_meta.RX_ID = 33;
dat.meta.RX_meta.meas_file = 'test_FilterEstCnoThresh.m';
dat.meta.RX_meta.obs_metadata{1} = '1';
dat.meta.RX_meta.obs_metadata{2} = '2';
dat.meta.RX_state_source = 'test_FilterEstCnoThresh.m';
dat.meta.TX_ID = makeGpsTXID(123, 13, 1); % ID 123, PRN 13, block II
dat.meta.TX_state_source = 'test_FilterEstCnoThresh.m';
dat.meta.gen_date = now;

% GPS DATA
gdat = makeGpsData; % empty
% add metadata
gdat.RX_meta.RX_ID = 33;
gdat.RX_meta.meas_file = 'test_FilterEstCnoThresh.m';
gdat.RX_meta.obs_metadata{1} = '1';
gdat.RX_meta.obs_metadata{2} = '2';

% add sample data
gdat.GPS_PRN = 13;
gdat.PRN_data{1}.epoch = 0:200;
gdat.PRN_data{1}.raw_SNR = 0:200;
gdat.PRN_data{1}.pseudorange = 0:200;
gdat.PRN_data{1}.doppler = 0:200;
gdat.PRN_data{1}.phase = 0:200;

cdat = [1:200; 1:200]; % bogus data for filtering
cexp = [22:200; 22:200]; % expected result

f = FilterEstCnoThresh(22);

CASE= 'filter(empty struct)';
try
    [cres, res, gres] = f.filter(cdat, dat); % no gps data input
    
    % if not gps data input then output should be []
    if ~isempty(gres)
        fprintf(1,'Failed test %d: %s - gps meas results should be empty.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check C/No data
    if isempty(cres) || any(any(cres ~= cexp))
        fprintf(1,'Failed test %d: %s - improper filtering.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check PP metadata
    if ~isstruct(res) || ~isfield(res,'meta') || ~isfield(res.meta,'RX_meta') || ...
        ~isfield(res.meta.RX_meta,'RX_ID') || res.meta.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_meta.meas_file ~= dat.meta.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.meta.RX_meta.obs_metadata) ~= 2 || ...
        dat.meta.RX_meta.obs_metadata{1} ~= res.meta.RX_meta.obs_metadata{1} || ...
        dat.meta.RX_meta.obs_metadata{2} ~= res.meta.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_state_source ~= dat.meta.RX_state_source
        fprintf(1,'Failed test %d: %s - RX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_ID ~= dat.meta.TX_ID.GPS_ID
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_PRN ~= dat.meta.TX_ID.GPS_PRN
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_PRN not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.block_type ~= dat.meta.TX_ID.block_type
        fprintf(1,'Failed test %d: %s - TX_ID.block_type not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_state_source ~= dat.meta.TX_state_source
        fprintf(1,'Failed test %d: %s - TX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res.meta, 'gen_date') || isempty(res.meta.gen_date)
        fprintf(1,'Failed test %d: %s - gen_date not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % ensure empty data
    if ~isfield(res, 'epoch') || ~isempty(res.epoch)
        fprintf(1,'Failed test %d: %s - incorrect epoch.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'TX_az') || ~isempty(res.TX_az)
        fprintf(1,'Failed test %d: %s - incorrect TX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'TX_el') || ~isempty(res.TX_el)
        fprintf(1,'Failed test %d: %s - incorrect TX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'RX_az') || ~isempty(res.RX_az)
        fprintf(1,'Failed test %d: %s - incorrect RX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'RX_el') || ~isempty(res.RX_el)
        fprintf(1,'Failed test %d: %s - incorrect RX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'range') || ~isempty(res.range)
        fprintf(1,'Failed test %d: %s - incorrect range.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'range_rate') || ~isempty(res.range_rate)
        fprintf(1,'Failed test %d: %s - incorrect range_rate.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'GPS_yaw') || ~isempty(res.GPS_yaw)
        fprintf(1,'Failed test %d: %s - incorrect GPS_yaw.\n',TESTNUM,CASE);
        fail = 1;
    end

catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

expdat1 = [];

%% TEST 4 - filter - filter all data, out of range threshold
TESTNUM = 4;
fprintf(1,'\nExecuting test_FilterEstCnoThresh test %d:\n\n',TESTNUM);

% GPS DATA
gdat = makeGpsData; % empty
% add metadata
gdat.RX_meta.RX_ID = 33;
gdat.RX_meta.meas_file = 'test_FilterEstCnoThresh.m';
gdat.RX_meta.obs_metadata{1} = '1';
gdat.RX_meta.obs_metadata{2} = '2';

% add sample data, 2 PRNs
gdat.GPS_PRN = [4 5];

epoch_base_utc = datenum('1 Jan 2012');
epoch_base_gps = convertTime('GPS','UTC',epoch_base_utc);

gdat.PRN_data{1}.epoch = (0:200)+epoch_base_gps;
gdat.PRN_data{1}.raw_SNR = 0:200;
gdat.PRN_data{1}.pseudorange = 0:200;
gdat.PRN_data{1}.doppler = 0:200;
gdat.PRN_data{1}.phase = 0:200;

gdat.PRN_data{2}.epoch = (0:200)+epoch_base_gps;
gdat.PRN_data{2}.raw_SNR = 0:200;
gdat.PRN_data{2}.pseudorange = 0:200;
gdat.PRN_data{2}.doppler = 0:200;
gdat.PRN_data{2}.phase = 0:200;

% PHYS PARAMS DATA
dat = makePpData; % empty
% add metadata
dat.meta.RX_meta.RX_ID = 33;
dat.meta.RX_meta.meas_file = 'test_FilterEstCnoThresh.m';
dat.meta.RX_meta.obs_metadata{1} = '1';
dat.meta.RX_meta.obs_metadata{2} = '2';
dat.meta.RX_state_source = 'test_FilterEstCnoThresh.m';
dat.meta.TX_ID = makeGpsTXID(123, 4, 1); % ID 123, PRN 4, block II
dat.meta.TX_state_source = 'test_FilterEstCnoThresh.m';
dat.meta.gen_date = now;
% add data (different length than gps data)
dat.epoch = (1:180)+epoch_base_utc;
dat.TX_az = (1:180)*d2r; % rad
dat.TX_el = (1:180)*d2r;
dat.RX_az = (1:180)*d2r;
dat.RX_el = (1:180)*d2r;
dat.range = 1:180;
dat.range_rate = 1:180;
dat.GPS_yaw = 1:180;

cdat = [dat.epoch; 1:length(dat.epoch)]; % bogus data for filtering

f = FilterEstCnoThresh(length(dat.epoch)+1); % out of range value means filter all

cexp = zeros(2,0); % expected result

% phys param expected data
expdatdef = [];
expdatang = (expdatdef)*d2r;
expdatepoch = (expdatdef)+epoch_base_utc; % UTC

% gps expected data
expdat1 = []; % note gps times cut down to match phys param times
expdat1epoch = (expdat1)+epoch_base_gps; % GPS


CASE= 'filter(all filtered)';
try
    [cres, res, gres] = f.filter(cdat, dat, gdat);
    
    % check C/No data
    if ~isempty(cres) || any(size(cres) ~= size(cexp))
        fprintf(1,'Failed test %d: %s - improper filtering.\n',TESTNUM,CASE);
        fail = 1;
    end
        
    % check PP metadata
    if ~isstruct(res) || ~isfield(res,'meta') || ~isfield(res.meta,'RX_meta') || ...
        ~isfield(res.meta.RX_meta,'RX_ID') || res.meta.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_meta.meas_file ~= dat.meta.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.meta.RX_meta.obs_metadata) ~= 2 || ...
        dat.meta.RX_meta.obs_metadata{1} ~= res.meta.RX_meta.obs_metadata{1} || ...
        dat.meta.RX_meta.obs_metadata{2} ~= res.meta.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_state_source ~= dat.meta.RX_state_source
        fprintf(1,'Failed test %d: %s - RX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_ID ~= dat.meta.TX_ID.GPS_ID
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_PRN ~= dat.meta.TX_ID.GPS_PRN
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_PRN not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.block_type ~= dat.meta.TX_ID.block_type
        fprintf(1,'Failed test %d: %s - TX_ID.block_type not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_state_source ~= dat.meta.TX_state_source
        fprintf(1,'Failed test %d: %s - TX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res.meta, 'gen_date') || isempty(res.meta.gen_date)
        fprintf(1,'Failed test %d: %s - gen_date not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % ensure empty data
    if ~isfield(res, 'epoch') || ~isempty(res.epoch)
        fprintf(1,'Failed test %d: %s - incorrect epoch.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'TX_az') || ~isempty(res.TX_az)
        fprintf(1,'Failed test %d: %s - incorrect TX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'TX_el') || ~isempty(res.TX_el)
        fprintf(1,'Failed test %d: %s - incorrect TX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'RX_az') || ~isempty(res.RX_az)
        fprintf(1,'Failed test %d: %s - incorrect RX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'RX_el') || ~isempty(res.RX_el)
        fprintf(1,'Failed test %d: %s - incorrect RX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'range') || ~isempty(res.range)
        fprintf(1,'Failed test %d: %s - incorrect range.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'range_rate') || ~isempty(res.range_rate)
        fprintf(1,'Failed test %d: %s - incorrect range_rate.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res, 'GPS_yaw') || ~isempty(res.GPS_yaw)
        fprintf(1,'Failed test %d: %s - incorrect GPS_yaw.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check gps data
    % check metadata
    if ~isstruct(gres) || ~isfield(gres,'RX_meta') || ...
        ~isfield(gres.RX_meta,'RX_ID') || gres.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if gres.RX_meta.meas_file ~= gdat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(gres.RX_meta.obs_metadata) ~= 2 || ...
        gdat.RX_meta.obs_metadata{1} ~= gres.RX_meta.obs_metadata{1} || ...
        gdat.RX_meta.obs_metadata{2} ~= gres.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check for empty data
    if any(gres.GPS_PRN ~= -1) || length(gres.PRN_data) ~= 1 ...
            || ~isempty(gres.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - completely filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(gres.PRN_data{1}.epoch ~= expdat1) || ...
        any(gres.PRN_data{1}.raw_SNR ~= expdat1) || ...
        any(gres.PRN_data{1}.pseudorange ~= expdat1) || ...
        any(gres.PRN_data{1}.doppler ~= expdat1) || ...
        any(gres.PRN_data{1}.phase ~= expdat1)
        fprintf(1,'Failed test %d: %s - improperly filtered 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check gps data
    % check metadata
    if ~isstruct(gres) || ~isfield(gres,'RX_meta') || ...
        ~isfield(gres.RX_meta,'RX_ID') || gres.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if gres.RX_meta.meas_file ~= gdat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(gres.RX_meta.obs_metadata) ~= 2 || ...
        gdat.RX_meta.obs_metadata{1} ~= gres.RX_meta.obs_metadata{1} || ...
        gdat.RX_meta.obs_metadata{2} ~= gres.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check data
    if any(gres.GPS_PRN ~= -1) || length(gres.PRN_data) ~= 1 ...
            || ~isempty(gres.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - filtered gps data unexpectedly populated.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(gres.PRN_data{1}.epoch ~= expdat1epoch) || ...
        any(gres.PRN_data{1}.raw_SNR ~= expdat1) || ...
        any(gres.PRN_data{1}.pseudorange ~= expdat1) || ...
        any(gres.PRN_data{1}.doppler ~= expdat1) || ...
        any(gres.PRN_data{1}.phase ~= expdat1)
        fprintf(1,'Failed test %d: %s - improperly filtered 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

%% TEST 5 - filter - partial filter
TESTNUM = 5;
fprintf(1,'\nExecuting test_FilterEstCnoThresh test %d:\n\n',TESTNUM);

% use the same data as above, change the filter

% phys param expected data
expdatdef = 90:180;
expdatang = (expdatdef)*d2r;
expdatepoch = (expdatdef)+epoch_base_utc; % UTC

% gps expected data
expdat1 = [90:180]; % note gps times cut down to match phys param times
expdat1epoch = (expdat1)+epoch_base_gps; % GPS

cexp = [expdatepoch; expdatdef];

f = FilterEstCnoThresh(90);

CASE= 'filter(partial filtered)';
try
    [cres, res, gres] = f.filter(cdat, dat, gdat);
    
    % check C/No data
    if isempty(cres) || any(size(cres) ~= size(cexp))
        fprintf(1,'Failed test %d: %s - improper filtering.\n',TESTNUM,CASE);
        fail = 1;
    end
        
    % check metadata
    if ~isstruct(res) || ~isfield(res,'meta') || ~isfield(res.meta,'RX_meta') || ...
        ~isfield(res.meta.RX_meta,'RX_ID') || res.meta.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_meta.meas_file ~= dat.meta.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.meta.RX_meta.obs_metadata) ~= 2 || ...
        dat.meta.RX_meta.obs_metadata{1} ~= res.meta.RX_meta.obs_metadata{1} || ...
        dat.meta.RX_meta.obs_metadata{2} ~= res.meta.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.RX_state_source ~= dat.meta.RX_state_source
        fprintf(1,'Failed test %d: %s - RX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_ID ~= dat.meta.TX_ID.GPS_ID
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.GPS_PRN ~= dat.meta.TX_ID.GPS_PRN
        fprintf(1,'Failed test %d: %s - TX_ID.GPS_PRN not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_ID.block_type ~= dat.meta.TX_ID.block_type
        fprintf(1,'Failed test %d: %s - TX_ID.block_type not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.meta.TX_state_source ~= dat.meta.TX_state_source
        fprintf(1,'Failed test %d: %s - TX_state_source not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isfield(res.meta, 'gen_date') || isempty(res.meta.gen_date)
        fprintf(1,'Failed test %d: %s - gen_date not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check data    
    if any(res.epoch ~= expdatepoch) || ...
        any(res.TX_az ~= expdatang) || ...
        any(res.RX_az ~= expdatang) || ...
        any(res.TX_el ~= expdatang) || ...
        any(res.RX_el ~= expdatang) || ...
        any(res.range ~= expdat1) || ...
        any(res.range_rate ~= expdat1) || ...
        any(res.GPS_yaw ~= expdat1)
        fprintf(1,'Failed test %d: %s - improperly filtered 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    
    % check gps data
    % check metadata
    if ~isstruct(gres) || ~isfield(gres,'RX_meta') || ...
        ~isfield(gres.RX_meta,'RX_ID') || gres.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if gres.RX_meta.meas_file ~= gdat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(gres.RX_meta.obs_metadata) ~= 2 || ...
        gdat.RX_meta.obs_metadata{1} ~= gres.RX_meta.obs_metadata{1} || ...
        gdat.RX_meta.obs_metadata{2} ~= gres.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check data
    if any(gres.GPS_PRN ~= 4) || length(gres.PRN_data) ~= 1 ...
            || isempty(gres.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(gres.PRN_data{1}.epoch ~= expdat1epoch) || ...
        any(gres.PRN_data{1}.raw_SNR ~= expdat1) || ...
        any(gres.PRN_data{1}.pseudorange ~= expdat1) || ...
        any(gres.PRN_data{1}.doppler ~= expdat1) || ...
        any(gres.PRN_data{1}.phase ~= expdat1)
        fprintf(1,'Failed test %d: %s - improperly filtered 1.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

clear f;

