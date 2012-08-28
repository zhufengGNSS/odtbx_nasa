function fail = FilterPpExplicitTimes_test()
% Unit test case for the FilterPpExplicitTimes class.
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

d2r = pi/180;

%% TEST 1 - constructor
TESTNUM = 1;
fprintf(1,'\nExecuting FilterPpExplicitTimes test %d:\n\n',TESTNUM);

CASE= 'no arg';
try
    f = FilterPpExplicitTimes;
    if ~isempty(f.filter_times)
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
    f = FilterPpExplicitTimes(99);
    if f.filter_times ~= 99
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;


CASE= 'with arg 3';
try
    f = FilterPpExplicitTimes(99:101);
    if any(f.filter_times ~= 99:101)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;


%% TEST 2 - getInd, filter, filterMeas - empty arg
TESTNUM = 2;
fprintf(1,'\nExecuting FilterPpExplicitTimes test %d:\n\n',TESTNUM);

f = FilterPpExplicitTimes(1:100);

CASE= 'getInd(empty arg)';
try
    res = f.getInd([]); %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(empty arg)';
try
    res = f.filter([]);
    if ~isfield(res,'epoch')
        fprintf(1,'Failed test %d: %s - missing field.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(res.epoch)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end


CASE= 'filterMeas(empty arg) 1';
try
    [res] = f.filterMeas([]);
    if ~isfield(res,'epoch')
        fprintf(1,'Failed test %d: %s - missing field.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(res.epoch)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'filterMeas(empty arg) 2';
try
    [res, gres] = f.filterMeas([]);
    if ~isfield(res,'epoch')
        fprintf(1,'Failed test %d: %s - missing field.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(res.epoch)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isempty(gres)
        fprintf(1,'Failed test %d: %s - bad gres.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

clear f;


%% TEST 3 - getInd, filter & filterMeas - empty phys params data struct
TESTNUM = 3;
fprintf(1,'\nExecuting FilterPpExplicitTimes test %d:\n\n',TESTNUM);

f = FilterPpExplicitTimes(-100:-1);

% physical param data
dat = makePpData; % empty
% add metadata
dat.meta.RX_meta.RX_ID = 33;
dat.meta.RX_meta.meas_file = 'FilterPpExplicitTimes.m';
dat.meta.RX_meta.obs_metadata{1} = '1';
dat.meta.RX_meta.obs_metadata{2} = '2';
dat.meta.RX_state_source = 'FilterPpExplicitTimes.m';
dat.meta.TX_ID = makeGpsTXID(123, 13, 1); % ID 123, PRN 13, block II
dat.meta.TX_state_source = 'FilterPpExplicitTimes.m';
dat.meta.gen_date = now;

% GPS DATA
gdat = makeGpsData; % empty
% add metadata
gdat.RX_meta.RX_ID = 33;
gdat.RX_meta.meas_file = 'FilterPpExplicitTimes.m';
gdat.RX_meta.obs_metadata{1} = '1';
gdat.RX_meta.obs_metadata{2} = '2';

% add sample data
gdat.GPS_PRN = 13;
gdat.PRN_data{1}.epoch = 0:200;
gdat.PRN_data{1}.raw_SNR = 0:200;
gdat.PRN_data{1}.pseudorange = 0:200;
gdat.PRN_data{1}.doppler = 0:200;
gdat.PRN_data{1}.phase = 0:200;

expdat1 = []; % expected gps data results

CASE= 'getInd(empty struct)';
try
    res = f.getInd(dat);
    if ~isempty(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'filter(empty struct)';
try
    res = f.filter(dat);
    
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

CASE= 'filterMeas(empty struct)';
try
    [res, gres] = f.filterMeas(dat, gdat);
    
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

catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

clear f;

%% TEST 4 - getInd, filter, FilterMeas - filter all data
TESTNUM = 4;
fprintf(1,'\nExecuting FilterPpExplicitTimes test %d:\n\n',TESTNUM);

% GPS DATA
gdat = makeGpsData; % empty
% add metadata
gdat.RX_meta.RX_ID = 33;
gdat.RX_meta.meas_file = 'FilterPpExplicitTimes.m';
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
dat.meta.RX_meta.meas_file = 'FilterPpExplicitTimes.m';
dat.meta.RX_meta.obs_metadata{1} = '1';
dat.meta.RX_meta.obs_metadata{2} = '2';
dat.meta.RX_state_source = 'FilterPpExplicitTimes.m';
dat.meta.TX_ID = makeGpsTXID(123, 4, 1); % ID 123, PRN 4, block II
dat.meta.TX_state_source = 'FilterPpExplicitTimes.m';
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

% phys param expected data
expdatdef = []; %#ok<NASGU>
expdatang = []; %#ok<NASGU>

% gps expected data
expdat1 = [];

f = FilterPpExplicitTimes(dat.epoch - 200); % no times in common with dat.epoch = filter all

CASE= 'getInd(all filtered)';
try
    res = f.getInd(dat);
    if ~isempty(res)
        fprintf(1,'Failed test %d: %s - shouldn''t be not empty.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'filter(all filtered)';
try
    res = f.filter(dat);
        
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


CASE= 'filterMeas(all filtered)';
try
    [res, gres] = f.filterMeas(dat,gdat);
        
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
    
    
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

%% TEST 5 - getInd & filter - partial filter
TESTNUM = 5;
fprintf(1,'\nExecuting FilterPpExplicitTimes test %d:\n\n',TESTNUM);

% use the same data as above, change the filter

% phys param expected data
expdatdef = 90:180;
expdatang = (expdatdef)*d2r;
expdatepoch = (expdatdef)+epoch_base_utc; % UTC

% gps expected data
expdat1 = 90:180; % note gps times cut down to match phys param times
expdat1epoch = (expdat1)+epoch_base_gps; % GPS

f = FilterPpExplicitTimes(expdatepoch); % UTC times for phys meas

CASE= 'getInd(partial filtered)';
try
    res = f.getInd(dat);
    if isempty(res)
        fprintf(1,'Failed test %d: %s - output mismatch.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'filter(partial filtered)';
try
    res = f.filter(dat);
        
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
    if any(expdatepoch ~= res.epoch)
        fprintf(1,'Failed test %d: %s - bad epoch.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatang ~= res.TX_az)
        fprintf(1,'Failed test %d: %s - bad angles, TX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatang ~= res.RX_az)
        fprintf(1,'Failed test %d: %s - bad angles, RX_az.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatang ~= res.TX_el)
        fprintf(1,'Failed test %d: %s - bad angles, TX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatang ~= res.RX_el)
        fprintf(1,'Failed test %d: %s - bad angles, RX_el.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatdef ~= res.range)
        fprintf(1,'Failed test %d: %s - bad data, range.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatdef ~= res.range_rate)
        fprintf(1,'Failed test %d: %s - bad data, range_rate.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(expdatdef ~= res.GPS_yaw)
        fprintf(1,'Failed test %d: %s - bad data, GPS_yaw.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end


CASE= 'filterMeas(partial filtered)';
try
    [res, gres] = f.filterMeas(dat,gdat);
        
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

