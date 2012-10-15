function [fail] = rinexo2gpsdata_test()
% Unit test for the rinexo2gps function.
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
fprintf(1,'Testing rinexo2gps...\n');

fs = filesep;

%% TEST 1: Nominal file load with struct returned
TEST = 1;
fprintf(1,'Test %d:\n',TEST);

CASE = 'Nominal struct, entire file';

try
    fprintf(1,'Test %d, (%s), reading file...\n',TEST,CASE);
    gps_data = rinexo2gpsdata(['DataFiles' fs 'RinexOTestFile.rnx'], 99);
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected exception.\n',TEST,CASE);
    fail = 1;
end

if gps_data.RX_meta.RX_ID ~= 99
    fprintf(1,'Test %d, (%s), bad RX_ID.\n',TEST,CASE);
    fail = 1;
end

if length(gps_data.GPS_PRN) ~= 22
    fprintf(1,'Test %d, (%s), bad gps_data.GPS_PRN length.\n',TEST,CASE);
    fail = 1;
end

prnind = find(gps_data.GPS_PRN == 2);

if length(gps_data.PRN_data{prnind}.epoch) ~= 13255
    fprintf(1,'Test %d, (%s), failed data length check.\n',TEST,CASE);
    fail = 1;
end

dv2 = datevec(gps_data.PRN_data{prnind}.epoch(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end

% compare against values from file
if (abs(gps_data.PRN_data{prnind}.pseudorange(end)-83632686.619) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.doppler(end)-(-1637.458)) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.phase(end)-0) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.raw_SNR(end)-35.193) > 1e-3)
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

% save off for comparison below
full_data = gps_data;

%% TEST 2: filtered file load with struct returned
TEST = 2;
fprintf(1,'Test %d:\n',TEST);

CASE = 'Filtered struct, entire file';

f = FilterGpsPrn('include',2);

try
    fprintf(1,'Test %d, (%s), reading file...\n',TEST,CASE);
    gps_data = rinexo2gpsdata(['DataFiles' fs 'RinexOTestFile.rnx'], 99, f);
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected exception.\n',TEST,CASE);
    fail = 1;
end

if gps_data.RX_meta.RX_ID ~= 99
    fprintf(1,'Test %d, (%s), bad RX_ID.\n',TEST,CASE);
    fail = 1;
end

if length(gps_data.GPS_PRN) ~= 1
    fprintf(1,'Test %d, (%s), bad gps_data.GPS_PRN length.\n',TEST,CASE);
    fail = 1;
end

prnind = find(gps_data.GPS_PRN == 2);

if length(gps_data.PRN_data{prnind}.epoch) ~= 13255
    fprintf(1,'Test %d, (%s), failed data length check.\n',TEST,CASE);
    fail = 1;
end

dv2 = datevec(gps_data.PRN_data{prnind}.epoch(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end

% compare against values from file
if (abs(gps_data.PRN_data{prnind}.pseudorange(end)-83632686.619) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.doppler(end)-(-1637.458)) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.phase(end)-0) > 1e-3) || ...
        (abs(gps_data.PRN_data{prnind}.raw_SNR(end)-35.193) > 1e-3)
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

%% TEST 3: filtered file load, write to file
TEST = 3;
fprintf(1,'Test %d:\n',TEST);

CASE = 'Filtered struct, write to file';

f = FilterGpsPrn('include',2);

try
    fprintf(1,'Test %d, (%s), reading file...\n',TEST,CASE);
    % pretty much guarantees writing to a single file
    fn = rinexo2gpsdata(['DataFiles' fs 'RinexOTestFile.rnx'], 99, f, 10000);
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected exception.\n',TEST,CASE);
    fail = 1;
end

if length(fn) ~= 1 || isempty(fn) || isempty(fn{1})
    fprintf(1,'Test %d, (%s), didn''t write a file.\n',TEST,CASE);
    fail = 1;
end

if ~exist(fn{1},'file')
    fprintf(1,'Test %d, (%s), invalid returned filename.\n',TEST,CASE);
    fail = 1;
end

gps_data = load(fn{1});

if gps_data.gps_dat.RX_meta.RX_ID ~= 99
    fprintf(1,'Test %d, (%s), bad RX_ID.\n',TEST,CASE);
    fail = 1;
end

if length(gps_data.gps_dat.GPS_PRN) ~= 1
    fprintf(1,'Test %d, (%s), bad gps_data.GPS_PRN length.\n',TEST,CASE);
    fail = 1;
end

prnind = find(gps_data.gps_dat.GPS_PRN == 2);

if length(gps_data.gps_dat.PRN_data{prnind}.epoch) ~= 13255
    fprintf(1,'Test %d, (%s), failed data length check.\n',TEST,CASE);
    fail = 1;
end

dv2 = datevec(gps_data.gps_dat.PRN_data{prnind}.epoch(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end

% compare against values from file
if (abs(gps_data.gps_dat.PRN_data{prnind}.pseudorange(end)-83632686.619) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.doppler(end)-(-1637.458)) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.phase(end)-0) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.raw_SNR(end)-35.193) > 1e-3)
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

% clean up
if ~isempty(fn)
    fprintf(1,'Test %d, (%s), removing test files.\n',TEST,CASE);
    for i = 1:length(fn)
        delete(fn{i});
    end
end

clear fn;

%% TEST 4: unfiltered file load, write to files
TEST = 4;
fprintf(1,'Test %d:\n',TEST);

CASE = 'Unfiltered struct, write to files';

try
    fprintf(1,'Test %d, (%s), reading file...\n',TEST,CASE);
    % pretty much guarantees writing to multiple files, dir specified
    mydir = [pwd fs 'rinexo2gpsdata_test_DIR'];
    if ~exist(mydir,'dir')
        fprintf(1,'Test %d, (%s), creating test directory.\n',TEST,CASE);
        mkdir( mydir);
    end
    fn = rinexo2gpsdata(['DataFiles' fs 'RinexOTestFile.rnx'], 99, [], 1, mydir);
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected exception.\n',TEST,CASE);
    fail = 1;
end

fprintf(1,'Test %d, (%s), created %d files.\n',TEST,CASE,length(fn));

if isempty(fn) || isempty(fn{1})
    fprintf(1,'Test %d, (%s), didn''t write a file.\n',TEST,CASE);
    fail = 1;
end

if ~exist(fn{1},'file')
    fprintf(1,'Test %d, (%s), invalid returned filename.\n',TEST,CASE);
    fail = 1;
end

for i = 1:length(fn)
    if isempty(strfind(fn{i},mydir))
        fprintf(1,'Test %d, (%s), wrote to incorrect directory, expected %s, got %s.\n',TEST,CASE, mydir, fn{i});
        fail = 1;
    end
end

gps_data = load(fn{end});

if gps_data.gps_dat.RX_meta.RX_ID ~= 99
    fprintf(1,'Test %d, (%s), bad RX_ID.\n',TEST,CASE);
    fail = 1;
end

prnind = find(gps_data.gps_dat.GPS_PRN == 2);

dv2 = datevec(gps_data.gps_dat.PRN_data{prnind}.epoch(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end

% compare against values from file
if (abs(gps_data.gps_dat.PRN_data{prnind}.pseudorange(end)-83632686.619) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.doppler(end)-(-1637.458)) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.phase(end)-0) > 1e-3) || ...
        (abs(gps_data.gps_dat.PRN_data{prnind}.raw_SNR(end)-35.193) > 1e-3)
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

% no cleanup for test 5

%% TEST 5: compare read types, file vs nofile
TEST = 5;
fprintf(1,'Test %d:\n',TEST);

CASE = 'compare file vs nofile';

prns = [];
numpts = 0;

% open each returned filename, get the PRNs and get the number of time
% points for PRN==2.
for i = 1:length(fn)
    if ~exist(fn{1},'file')
        fprintf(1,'Test %d, (%s), invalid returned filename: %s.\n',TEST,CASE,fn{i});
        fail = 1;
    end
    
    if gps_data.gps_dat.RX_meta.RX_ID ~= 99
        fprintf(1,'Test %d, (%s), bad RX_ID.\n',TEST,CASE);
        fail = 1;
    end
    
    gps_data = load(fn{i});
    
    prns = union(gps_data.gps_dat.GPS_PRN, prns);
    
    prnind = find(gps_data.gps_dat.GPS_PRN == 2);
    numpts = numpts + length(gps_data.gps_dat.PRN_data{prnind}.epoch); %#ok<FNDSB>
end

if length(prns) ~= 22
    fprintf(1,'Test %d, (%s), bad gps_data.GPS_PRN length.\n',TEST,CASE);
    fail = 1;
end

if numpts ~= 13255
    fprintf(1,'Test %d, (%s), failed data length check.\n',TEST,CASE);
    fail = 1;
end

% clean up
if ~isempty(fn)
    fprintf(1,'Test %d, (%s), removing test files.\n',TEST,CASE);
    for i = 1:length(fn)
        delete(fn{i});
    end
end

if exist(mydir,'dir')
    fprintf(1,'Test %d, (%s), removing test directory.\n',TEST,CASE);
    rmdir( mydir);
end


