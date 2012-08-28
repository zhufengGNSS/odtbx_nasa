function [fail] = RinexOReader_test()
% Unit test for the RinexOReader class.
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
fprintf(1,'Testing RinexOReader...\n');

fs = filesep;

%% TEST 1: Load a file
TEST = 1;
CASE = 'Bad file load';
fprintf(1,'Test %d:\n',TEST);

try
    ror = RinexOReader(['ABCDEFG',fs,'RinexOTestFile.rnx']); %#ok<NASGU>
    
    fprintf(1,'Test %d, (%s), missed expected constructor exception.\n',TEST,CASE);
    fail = 1;
catch ex %#ok<NASGU>
    % expected
end
clear ror;

CASE = 'Nominal load';

try
    ror = RinexOReader(['DataFiles',fs,'RinexOTestFile.rnx']); %#ok<NASGU>
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

%% TEST 2: check pre-read status
TEST = 2;
CASE = 'pre-read status';
fprintf(1,'Test %d:\n',TEST);

try
    if ~(ror.isReady)
        fprintf(1,'Test %d, (%s), reader not ready.\n',TEST,CASE);
        fail = 1;
    end
    if ror.isDone
        fprintf(1,'Test %d, (%s), reader unexpectedly done.\n',TEST,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

%% TEST 3: read the header
TEST = 3;
CASE = 'read header';
fprintf(1,'Test %d:\n',TEST);

try
    ror.readHeader;
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

try
    if ~(ror.isReady)
        fprintf(1,'Test %d, (%s), reader not ready.\n',TEST,CASE);
        fail = 1;
    end
    if ror.isDone
        fprintf(1,'Test %d, (%s), reader unexpectedly done.\n',TEST,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

% check header data
if ror.numobs ~= 4
    fprintf(1,'Test %d, (%s), bad numobs.\n',TEST,CASE);
    fail = 1;
end
if length(ror.header) ~= 4
    fprintf(1,'Test %d, (%s), bad header lines.\n',TEST,CASE);
    fail = 1;
end

%% TEST 4: read the body
TEST = 4;
CASE = 'read body';
fprintf(1,'Test %d:\n',TEST);

try
    fprintf(1,'Test %d, (%s), readBody() call...\n',TEST,CASE);
    ror.readBody;
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

try
    if ror.isReady
        fprintf(1,'Test %d, (%s), reader still ready.\n',TEST,CASE);
        fail = 1;
    end
    if ~ror.isDone
        fprintf(1,'Test %d, (%s), reader not done.\n',TEST,CASE);
        fail = 1;
    end
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

% Check data sizes
if any(size(ror.timegps) ~= [1 83113])
    fprintf(1,'Test %d, (%s), bad timegps size.\n',TEST,CASE);
    fail = 1;
end
if any(size(ror.C1) ~= [32 83113])
    fprintf(1,'Test %d, (%s), bad C1 size.\n',TEST,CASE);
    fail = 1;
end
if any(size(ror.D1) ~= [32 83113])
    fprintf(1,'Test %d, (%s), bad D1 size.\n',TEST,CASE);
    fail = 1;
end
if any(size(ror.L1) ~= [32 83113])
    fprintf(1,'Test %d, (%s), bad L1 size.\n',TEST,CASE);
    fail = 1;
end
if any(size(ror.S1) ~= [32 83113])
    fprintf(1,'Test %d, (%s), bad S1 size.\n',TEST,CASE);
    fail = 1;
end
if any(size(ror.bias) ~= [1 83113])
    fprintf(1,'Test %d, (%s), bad bias size.\n',TEST,CASE);
    fail = 1;
end

% test data, first data point, no measurements
dv1 = datevec(ror.timegps(1)+723186);
t1 = [1980  1  6  0  0  1.75000000]; % from file
if any((dv1-t1) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 1.\n',TEST,CASE);
    fail = 1;
end
if any(ror.C1(1) ~= 0) || any(ror.D1(1) ~= 0) || ...
        any(ror.L1(1) ~= 0) || any(ror.S1(1) ~= 0) ||...
        any(ror.bias(1) ~= 0)
    fprintf(1,'Test %d, (%s), failed data check 1.\n',TEST,CASE);
    fail = 1;
end

dv2 = datevec(ror.timegps(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end
dind = [1 3:32];
if any(ror.C1(dind, end) ~= 0) || any(ror.D1(dind, end) ~= 0) || ...
        any(ror.L1(dind, end) ~= 0) || any(ror.S1(dind, end) ~= 0) ||...
        ror.bias(end) ~= 0
    fprintf(1,'Test %d, (%s), failed zeros data check 2.\n',TEST,CASE);
    fail = 1;
end
if abs(ror.C1(2, end)-83632686.619) > 1e-3 || ...
        abs(ror.D1(2, end)-(-1637.458)) > 1e-3 || ...
        abs(ror.L1(2, end)-0 > 1e-3 || ...
        abs(ror.S1(2, end)-35.193) > 1e-3 ||...
        abs(ror.bias(end)-0)) > 1e-3
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

%% TEST 5: readlimit adjust & multiple readBody calls
TEST = 5;
CASE = 'multiple readBody calls';
fprintf(1,'Test %d:\n',TEST);

clear ror;

try
    ror = RinexOReader(['DataFiles',fs,'RinexOTestFile.rnx']); %#ok<NASGU>
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected constructor exception.\n',TEST,CASE);
    fail = 1;
end

try
    ror.readlimit = 10000;
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected set readlimit exception.\n',TEST,CASE);
    fail = 1;
end

if ~(ror.isReady)
    fprintf(1,'Test %d, (%s), reader not ready.\n',TEST,CASE);
    fail = 1;
end
if ror.isDone
    fprintf(1,'Test %d, (%s), reader unexpectedly done.\n',TEST,CASE);
    fail = 1;
end

try
    ror.readHeader;
catch ex %#ok<NASGU>
    fprintf(1,'Test %d, (%s), unexpected readHeader exception.\n',TEST,CASE);
    fail = 1;
end

ctr = 1;
numpts = 0;
while (ror.isReady && ~ror.isDone)
    try
        fprintf(1,'Test %d, (%s), readBody() call #%d...\n',TEST,CASE,ctr);
        ror.readBody;
    catch ex %#ok<NASGU>
        fprintf(1,'Test %d, (%s), unexpected readBody exception.\n',TEST,CASE);
        fail = 1;
    end
    
    numpts = numpts + length(ror.timegps);

    if ctr == 1
        fprintf(1,'Test %d, (%s), first data check.\n',TEST,CASE);
        % test data, first data point, no measurements
        dv1 = datevec(ror.timegps(1)+723186);
        t1 = [1980  1  6  0  0  1.75000000]; % from file
        if any((dv1-t1) > 1e-3)
            fprintf(1,'Test %d, (%s), failed timegps check 1.\n',TEST,CASE);
            fail = 1;
        end
        if any(ror.C1(1) ~= 0) || any(ror.D1(1) ~= 0) || ...
                any(ror.L1(1) ~= 0) || any(ror.S1(1) ~= 0) ||...
                any(ror.bias(1) ~= 0)
            fprintf(1,'Test %d, (%s), failed data check 1.\n',TEST,CASE);
            fail = 1;
        end
    end
    ctr = ctr+1;
end % while

% final time point check after last call
fprintf(1,'Test %d, (%s), last data check.\n',TEST,CASE);
dv2 = datevec(ror.timegps(end)+723186);
t2 = [2015  4 30 14 39 51.10199164 ]; % from file
if any((dv2-t2) > 1e-3)
    fprintf(1,'Test %d, (%s), failed timegps check 2.\n',TEST,CASE);
    fail = 1;
end
dind = [1 3:32];
if any(ror.C1(dind, end) ~= 0) || any(ror.D1(dind, end) ~= 0) || ...
        any(ror.L1(dind, end) ~= 0) || any(ror.S1(dind, end) ~= 0) ||...
        ror.bias(end) ~= 0
    fprintf(1,'Test %d, (%s), failed zeros data check 2.\n',TEST,CASE);
    fail = 1;
end
if abs(ror.C1(2, end)-83632686.619) > 1e-3 || ...
        abs(ror.D1(2, end)-(-1637.458)) > 1e-3 || ...
        abs(ror.L1(2, end)-0 > 1e-3 || ...
        abs(ror.S1(2, end)-35.193) > 1e-3 ||...
        abs(ror.bias(end)-0)) > 1e-3
    fprintf(1,'Test %d, (%s), failed data check 2.\n',TEST,CASE);
    fail = 1;
end

if ctr ~= 10 % 83113/10000 is 9 calls, +1 for the counter
    fprintf(1,'Test %d, (%s), improper number of calls to read data.\n',TEST,CASE);
    fail = 1;
end

% check amount of data read
if numpts ~= 83113
    fprintf(1,'Test %d, (%s), improper number of times read.\n',TEST,CASE);
    fail = 1;
end

%% Clean up
if nargout < 2
    clear ror;
end
