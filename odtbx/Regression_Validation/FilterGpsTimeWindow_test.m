function fail = FilterGpsTimeWindow_test()
% Unit test case for the FilterGpsTimeWindow class.
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

%% TEST 1 - constructor
TESTNUM = 1;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

CASE= 'no arg';
try
    f = FilterGpsTimeWindow;
    if f.early ~= 0
        fprintf(1,'Failed test %d: %s - bad early value.\n',TESTNUM,CASE);
        fail = 1;
    end
    if f.late ~= 0
        fprintf(1,'Failed test %d: %s - bad late value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

CASE= 'with arg';
try
    f = FilterGpsTimeWindow(1,2);
    if f.early ~= 1
        fprintf(1,'Failed test %d: %s - bad early value.\n',TESTNUM,CASE);
        fail = 1;
    end
    if f.late ~= 2
        fprintf(1,'Failed test %d: %s - bad late value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

CASE= 'reversed arg';
try
    f = FilterGpsTimeWindow(2,1);
    if f.early ~= 1
        fprintf(1,'Failed test %d: %s - bad early value.\n',TESTNUM,CASE);
        fail = 1;
    end
    if f.late ~= 2
        fprintf(1,'Failed test %d: %s - bad late value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

CASE= 'bad arg1';
try
    f = FilterGpsTimeWindow([1:99]');
    fprintf(1,'Failed test %d: %s, expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'bad arg2';
try
    f = FilterGpsTimeWindow([1:2],[1:2]);
    fprintf(1,'Failed test %d: %s, expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'bad arg3';
try
    f = FilterGpsTimeWindow(2,3,4);
    fprintf(1,'Failed test %d: %s, expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end
clear f;

%% TEST 2 - getInd & filter - empty arg
TESTNUM = 2;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

f = FilterGpsTimeWindow(1,2);

CASE= 'getInd(empty arg)';
try
    res = f.getInd([]);
    fprintf(1,'Failed test %d: %s, missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'getInd(empty arg, bad prn)';
try
    res = f.getInd([],44);
    fprintf(1,'Failed test %d: %s, missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(empty arg)';
try
    res = f.filter([]);
    if FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

clear f;


%% TEST 3 - getInd & filter - empty gps data struct
TESTNUM = 3;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

f = FilterGpsTimeWindow(1,2);
dat = makeGpsData; % empty
% add metadata
dat.RX_meta.RX_ID = TESTNUM;
dat.RX_meta.meas_file = 'FilterGpsTimeWindow.m';
dat.RX_meta.obs_metadata{1} = '1';
dat.RX_meta.obs_metadata{2} = '2';

CASE= 'getInd(empty struct)';
try
    res = f.getInd(dat,1);
    if ~isempty(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'getInd(empty struct, bad prn)';
try
    res = f.getInd(dat,44);
    fprintf(1,'Failed test %d: %s - missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(empty struct)';
try
    res = f.filter(dat);
    if FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
    if ~isstruct(res) || ~isfield(res,'RX_meta') || ...
        ~isfield(res.RX_meta,'RX_ID') || res.RX_meta.RX_ID ~= TESTNUM
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.RX_meta.meas_file ~= dat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.RX_meta.obs_metadata) ~= 2 || ...
        dat.RX_meta.obs_metadata{1} ~= res.RX_meta.obs_metadata{1} || ...
        dat.RX_meta.obs_metadata{2} ~= res.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

clear f;

%% TEST 4 - getInd & filter - filter all of one PRN
TESTNUM = 4;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

f = FilterGpsTimeWindow( 101,105 );

dat = makeGpsData; % empty
% add metadata
dat.RX_meta.RX_ID = 33;
dat.RX_meta.meas_file = 'FilterGpsTimeWindow.m';
dat.RX_meta.obs_metadata{1} = '1';
dat.RX_meta.obs_metadata{2} = '2';

% add sample data
dat.GPS_PRN = 4;
dat.PRN_data{1}.epoch = 0:100;
dat.PRN_data{1}.raw_SNR = 0:100;
dat.PRN_data{1}.pseudorange = 0:100;
dat.PRN_data{1}.doppler = 0:100;
dat.PRN_data{1}.phase = 0:100;


CASE= 'getInd(all filtered)';
try
    res = f.getInd(dat,1);
    if ~isempty(res)
        fprintf(1,'Failed test %d: %s - shouldn''t be not empty.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'getInd(all filtered, bad prn)';
try
    res = f.getInd(dat,44);
    fprintf(1,'Failed test %d: %s - missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(all filtered)';
try
    res = f.filter(dat);
    if FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - bad value.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check metadata
    if ~isstruct(res) || ~isfield(res,'RX_meta') || ...
        ~isfield(res.RX_meta,'RX_ID') || res.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.RX_meta.meas_file ~= dat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.RX_meta.obs_metadata) ~= 2 || ...
        dat.RX_meta.obs_metadata{1} ~= res.RX_meta.obs_metadata{1} || ...
        dat.RX_meta.obs_metadata{2} ~= res.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check for empty data
    if res.GPS_PRN ~= -1 || length(res.PRN_data) ~= 1 || ~isempty(res.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - completely filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;


%% TEST 5 - getInd & filter - partial filter, one PRN
TESTNUM = 5;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

% use the same data as above, change the filter
f = FilterGpsTimeWindow(1,100);
expind = 2:101; % expected indices
expdat = 1:100;

CASE= 'getInd(partial filtered)';
try
    res = f.getInd(dat,1);
    if res ~= expind
        fprintf(1,'Failed test %d: %s - output mismatch.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'getInd(partial filtered, bad prn)';
try
    res = f.getInd(dat,44);
    fprintf(1,'Failed test %d: %s - missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(partial filtered)';
try
    res = f.filter(dat);
    if ~FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - unexpectedly empty.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check metadata
    if ~isstruct(res) || ~isfield(res,'RX_meta') || ...
        ~isfield(res.RX_meta,'RX_ID') || res.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.RX_meta.meas_file ~= dat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.RX_meta.obs_metadata) ~= 2 || ...
        dat.RX_meta.obs_metadata{1} ~= res.RX_meta.obs_metadata{1} || ...
        dat.RX_meta.obs_metadata{2} ~= res.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check for empty data
    if res.GPS_PRN ~= 4 || length(res.PRN_data) ~= 1 || isempty(res.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - completely filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(res.PRN_data{1}.epoch ~= expdat) || ...
        any(res.PRN_data{1}.raw_SNR ~= expdat) || ...
        any(res.PRN_data{1}.pseudorange ~= expdat) || ...
        any(res.PRN_data{1}.doppler ~= expdat) || ...
        any(res.PRN_data{1}.phase ~= expdat)
        fprintf(1,'Failed test %d: %s - improperly filtered.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;

%% TEST 6 - getInd & filter - filter none, one PRN
TESTNUM = 6;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

% use the same data as above, change the filter
f = FilterGpsTimeWindow(-5,1300);
expind = 1:101; % expected indices
expdat = 0:100;

CASE= 'getInd(filter none)';
try
    res = f.getInd(dat,1);
    if res ~= expind
        fprintf(1,'Failed test %d: %s - output mismatch.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'getInd(filter none, bad prn)';
try
    res = f.getInd(dat,44);
    fprintf(1,'Failed test %d: %s - missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(filter none)';
try
    res = f.filter(dat);
    if ~FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - unexpectedly empty.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check metadata
    if ~isstruct(res) || ~isfield(res,'RX_meta') || ...
        ~isfield(res.RX_meta,'RX_ID') || res.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.RX_meta.meas_file ~= dat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.RX_meta.obs_metadata) ~= 2 || ...
        dat.RX_meta.obs_metadata{1} ~= res.RX_meta.obs_metadata{1} || ...
        dat.RX_meta.obs_metadata{2} ~= res.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check for empty data
    if res.GPS_PRN ~= 4 || length(res.PRN_data) ~= 1 || isempty(res.PRN_data{1}.epoch)
        fprintf(1,'Failed test %d: %s - completely filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(res.PRN_data{1}.epoch ~= expdat) || ...
        any(res.PRN_data{1}.raw_SNR ~= expdat) || ...
        any(res.PRN_data{1}.pseudorange ~= expdat) || ...
        any(res.PRN_data{1}.doppler ~= expdat) || ...
        any(res.PRN_data{1}.phase ~= expdat)
        fprintf(1,'Failed test %d: %s - improperly filtered.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;
%% TEST 7 - filter - multiple PRNs
TESTNUM = 7;
fprintf(1,'\nExecuting FilterGpsTimeWindow test %d:\n\n',TESTNUM);

% add more data
% add sample data
dat.GPS_PRN = [4 5];
dat.PRN_data{2}.epoch = 0:200;
dat.PRN_data{2}.raw_SNR = 0:200;
dat.PRN_data{2}.pseudorange = 0:200;
dat.PRN_data{2}.doppler = 0:200;
dat.PRN_data{2}.phase = 0:200;

% change the filter
f = FilterGpsTimeWindow(50,1200);
expind1 = 51:101; % expected indices
expdat1 = 50:100;
expind2 = 51:201; % expected indices
expdat2 = 50:200;



CASE= 'getInd(filter multiple PRNs)';
try
    res = f.getInd(dat,1);
    if res ~= expind1
        fprintf(1,'Failed test %d: %s - output mismatch 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    res = f.getInd(dat,2);
    if res ~= expind2
        fprintf(1,'Failed test %d: %s - output mismatch 2.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end

CASE= 'getInd(filter multiple PRNs, bad prn)';
try
    res = f.getInd(dat,44);
    fprintf(1,'Failed test %d: %s - missed expected exception.\n',TESTNUM,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end

CASE= 'filter(filter multiple PRNs)';
try
    res = f.filter(dat);
    if ~FilterGpsData.hasData(res)
        fprintf(1,'Failed test %d: %s - unexpectedly empty.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check metadata
    if ~isstruct(res) || ~isfield(res,'RX_meta') || ...
        ~isfield(res.RX_meta,'RX_ID') || res.RX_meta.RX_ID ~= 33
        fprintf(1,'Failed test %d: %s - RX_ID not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if res.RX_meta.meas_file ~= dat.RX_meta.meas_file
        fprintf(1,'Failed test %d: %s - meas_file not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    if length(dat.RX_meta.obs_metadata) ~= 2 || ...
        dat.RX_meta.obs_metadata{1} ~= res.RX_meta.obs_metadata{1} || ...
        dat.RX_meta.obs_metadata{2} ~= res.RX_meta.obs_metadata{2}
        fprintf(1,'Failed test %d: %s - obs_metadata not preserved.\n',TESTNUM,CASE);
        fail = 1;
    end
    % check for empty data
    if any(res.GPS_PRN ~= [4 5]) || length(res.PRN_data) ~= 2 ...
            || isempty(res.PRN_data{1}.epoch) || isempty(res.PRN_data{2}.epoch)
        fprintf(1,'Failed test %d: %s - completely filtered data not properly structured.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(res.PRN_data{1}.epoch ~= expdat1) || ...
        any(res.PRN_data{1}.raw_SNR ~= expdat1) || ...
        any(res.PRN_data{1}.pseudorange ~= expdat1) || ...
        any(res.PRN_data{1}.doppler ~= expdat1) || ...
        any(res.PRN_data{1}.phase ~= expdat1)
        fprintf(1,'Failed test %d: %s - improperly filtered 1.\n',TESTNUM,CASE);
        fail = 1;
    end
    if any(res.PRN_data{2}.epoch ~= expdat2) || ...
        any(res.PRN_data{2}.raw_SNR ~= expdat2) || ...
        any(res.PRN_data{2}.pseudorange ~= expdat2) || ...
        any(res.PRN_data{2}.doppler ~= expdat2) || ...
        any(res.PRN_data{2}.phase ~= expdat2)
        fprintf(1,'Failed test %d: %s - improperly filtered 2.\n',TESTNUM,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Failed test %d: %s, unexpected exception.\n',TESTNUM,CASE);
    fail = 1;
end
clear f;