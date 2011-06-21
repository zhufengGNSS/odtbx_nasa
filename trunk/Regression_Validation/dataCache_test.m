function [fail] = dataCache_test()
% test for dataCache.m
% Runs through the interface and tests the data cache.
% Starts with basic input checking (somewhat thorough, but not exhaustive).
% Then tests cache creation and simple isempty call.
% Exercises data add calls with good and bogus arguments.
% After data is added exercises get calls.
% Then goes on to test delete on the built-up cache.  Delete is exercised
% in several orders to ensure deletion order doesn't affect the data.
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

% no input arg
try
    dataCache();
    fail = 1;
    disp('dataCache_test: failed input arg test');
catch
    % expected
end

% no output arg
try
    dataCache('create');
    fail = 1;
    disp('dataCache_test: failed output arg test');
catch
    % expected
end

% test bogus arg
try
    dataCache('BOGUS!');
    fail = 1;
    disp('dataCache_test: failed bogus arg test');
catch
    % expected
end

% test create1
cache = [];
try
    cache = dataCache('create');
catch exception
    fail = 1;
    disp('dataCache_test: failed create1 test');
    exception.getReport
end

if ~isstruct(cache)
    fail = 1;
    disp('dataCache_test: failed create1b arg test');
end

% test create2
cache = [];
try
    cache = dataCache('CrEaTe');
catch exception
    fail = 1;
    disp('dataCache_test: failed create2 test');
    exception.getReport
end

if ~isstruct(cache)
    fail = 1;
    disp('dataCache_test: failed create2b arg test');
end

% test isempty
e = -1;
try
    e = dataCache('isempty',cache);
catch exception
    fail = 1;
    disp('dataCache_test: failed isempty1 test');
    exception.getReport
end
if e ~= 1
    fail = 1;
    disp('dataCache_test: failed isempty1b test');
end

% test isempty2
e = -1;
try
    e = dataCache('IsEmPTY',cache);
catch exception
    fail = 1;
    disp('dataCache_test: failed isempty2 test');
    exception.getReport
end
if e ~= 1
    fail = 1;
    disp('dataCache_test: failed isempty2b test');
end

% test isempty, bogus arg
e = -1;
boguscache = struct('sdf',[1:5],'iopaeudfuio',[1:10]);
try
    e = dataCache('IsEmPTY',boguscache);
    fail = 1;
    disp('dataCache_test: failed bogus isempty test');
catch %exception
    % expected
end

% test length of zero
e = -1;
try
    e = dataCache('length',cache);
catch exception
    fail = 1;
    disp('dataCache_test: failed lengthzero test');
    exception.getReport
end
if e ~= 0
    fail = 1;
    disp('dataCache_test: failed lengthzero test');
end

% test add 1
cache = dataCache('create');
try
    cache = dataCache('add',cache,'key1',[]);
catch exception
    fail = 1;
    disp('dataCache_test: failed add1 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed add1b test');
end
if dataCache('length',cache) ~= 1
    fail = 1;
    disp('dataCache_test: failed add1c test');
end

% test add 2
cache = dataCache('create');
try
    cache = dataCache('add',cache,'key1',[1:20]);
catch exception
    fail = 1;
    disp('dataCache_test: failed add2 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed add2b test');
end
if dataCache('length',cache) ~= 1
    fail = 1;
    disp('dataCache_test: failed add2c test');
end

% test add, bogus arg
boguscache = struct('sdf',[1:5],'iopaeudfuio',[1:10]);
try
    boguscache = dataCache('add',boguscache,'key1',[1:20]);
    fail = 1;
    disp('dataCache_test: failed bogus add test1');
catch %exception
    % expected
end

% test add, bogus arg 2
try
    cache = dataCache('addsadr',cache,'key1',[1:20]);
    fail = 1;
    disp('dataCache_test: failed bogus add test2');
catch %exception
    % expected
end

% test add, bogus arg 3
try
    cache = dataCache('add',cache,'',[1:20]);
    fail = 1;
    disp('dataCache_test: failed bogus add test3');
catch %exception
    % expected
end

% test add 3
cache = dataCache('create');
try
    cache = dataCache('add',cache,'key1',[1:10]);
    cache = dataCache('add',cache,'key2',[2:20]);
    cache = dataCache('add',cache,'key3',[3:30]);
catch exception
    fail = 1;
    disp('dataCache_test: failed add3 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed add3b test');
end
if dataCache('length',cache) ~= 3
    fail = 1;
    disp('dataCache_test: failed add3c test');
end

% test get 1
check = [];
try
    check = dataCache('get',cache,'key1');
catch exception
    fail = 1;
    disp('dataCache_test: failed get1 test');
    exception.getReport
end
if check ~= [1:10]
    fail = 1;
    disp('dataCache_test: failed get1a test');
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed get1b test');
end
if dataCache('length',cache) ~= 3
    fail = 1;
    disp('dataCache_test: failed get1c test');
end

% test get 2
check = [];
try
    check = dataCache('get',cache,'key2');
catch exception
    fail = 1;
    disp('dataCache_test: failed get2 test');
    exception.getReport
end
if check ~= [2:20]
    fail = 1;
    disp('dataCache_test: failed get2a test');
end

% test get 3
check = [];
try
    check = dataCache('get',cache,'key3');
catch exception
    fail = 1;
    disp('dataCache_test: failed get3 test');
    exception.getReport
end
if check ~= [3:30]
    fail = 1;
    disp('dataCache_test: failed get3a test');
end

% test get 4
check = [];
try
    check = dataCache('gET',cache,'key3');
catch exception
    fail = 1;
    disp('dataCache_test: failed get4 test');
    exception.getReport
end
if check ~= [3:30]
    fail = 1;
    disp('dataCache_test: failed get4a test');
end

% test get, missing key
check = [5:10];
try
    check = dataCache('gET',cache,'key999');
catch exception
    fail = 1;
    disp('dataCache_test: failed get missing key test');
    exception.getReport
end
if ~isempty(check)
    fail = 1;
    disp('dataCache_test: failed get missing key test');
end


% test get, bogus cache
boguscache = struct('sdf',[1:5],'iopaeudfuio',[1:10]);
try
    check = dataCache('get',boguscache,'key1');
    fail = 1;
    disp('dataCache_test: failed bogus get test');
catch %exception
    % expected
end

% check the cache before testing delete (the above operations shouldn't
% have affected the cache's data)
check = [];
try
    check = dataCache('get',cache,'key1');
catch exception
    fail = 1;
    disp('dataCache_test: failed get1 test (reprise)');
    exception.getReport
end
if check ~= [1:10]
    fail = 1;
    disp('dataCache_test: failed get1a test (reprise)');
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed get1b test (reprise)');
end
if dataCache('length',cache) ~= 3
    fail = 1;
    disp('dataCache_test: failed get1c test (reprise)');
end

% test delete, incorrect key (should be no change to cache)
try
    cache = dataCache('delete',cache,'key999');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete1 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete1b test');
end
if dataCache('length',cache) ~= 3
    fail = 1;
    disp('dataCache_test: failed delete1c test');
end
if dataCache('get',cache,'key3') ~= [3:30]
    fail = 1;
    disp('dataCache_test: failed delete1d test');
end

% test delete (2), last item
try
    cache = dataCache('delete',cache,'key3');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete2 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete2b test');
end
if dataCache('length',cache) ~= 2
    fail = 1;
    disp('dataCache_test: failed delete2c test');
end
if ~isempty(dataCache('get',cache,'key3'))
    fail = 1;
    disp('dataCache_test: failed delete2d test');
end
if dataCache('get',cache,'key1') ~= [1:10]
    fail = 1;
    disp('dataCache_test: failed delete2e test');
end
if dataCache('get',cache,'key2') ~= [2:20]
    fail = 1;
    disp('dataCache_test: failed delete2f test');
end

% test delete (3), last item again
try
    cache = dataCache('delete',cache,'key2');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete3 test');
    exception.getReport
end
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete3b test');
end
if dataCache('length',cache) ~= 1
    fail = 1;
    disp('dataCache_test: failed delete3c test');
end
if ~isempty(dataCache('get',cache,'key3'))
    fail = 1;
    disp('dataCache_test: failed delete3d test');
end

if ~isempty(dataCache('get',cache,'key2'))
    fail = 1;
    disp('dataCache_test: failed delete3e test');
end
if isempty(dataCache('get',cache,'key1'))
    fail = 1;
    disp('dataCache_test: failed delete3f test');
end
if dataCache('get',cache,'key1') ~= [1:10]
    fail = 1;
    disp('dataCache_test: failed delete3g test');
end

% test delete (4), the last item
try
    cache = dataCache('delete',cache,'key1');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete4 test');
    exception.getReport
end
if ~dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete4b test');
end
if dataCache('length',cache) ~= 0
    fail = 1;
    disp('dataCache_test: failed delete4c test');
end
if ~isempty(dataCache('get',cache,'key3'))
    fail = 1;
    disp('dataCache_test: failed delete4d test');
end
if ~isempty(dataCache('get',cache,'key2'))
    fail = 1;
    disp('dataCache_test: failed delete4e test');
end
if ~isempty(dataCache('get',cache,'key1'))
    fail = 1;
    disp('dataCache_test: failed delete4f test');
end

% test delete (5), on an empty cache
try
    cache = dataCache('delete',cache,'key1');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete5 test');
    exception.getReport
end

% rebuild the cache again
cache = [];
try
    cache = dataCache('create');
    cache = dataCache('add',cache,'key1',1);
    cache = dataCache('add',cache,'key2',2);
    cache = dataCache('add',cache,'key3',3);
catch exception
    fail = 1;
    disp('dataCache_test: failed rebuild 5-6 test');
    exception.getReport
end


% test delete (6), the first item
try
    cache = dataCache('delete',cache,'key1');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete6 test');
    exception.getReport
end

% check that the other items are there
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete6b test');
end
if dataCache('length',cache) ~= 2
    fail = 1;
    disp('dataCache_test: failed delete6c test');
end
if isempty(dataCache('get',cache,'key3'))
    fail = 1;
    disp('dataCache_test: failed delete6d test');
end
if isempty(dataCache('get',cache,'key2'))
    fail = 1;
    disp('dataCache_test: failed delete6e test');
end
if ~isempty(dataCache('get',cache,'key1'))
    fail = 1;
    disp('dataCache_test: failed delete6f test');
end
if dataCache('get',cache,'key3') ~= 3
    fail = 1;
    disp('dataCache_test: failed delete6g test');
end
if dataCache('get',cache,'key2') ~= 2
    fail = 1;
    disp('dataCache_test: failed delete6g test');
end


% rebuild the cache again
cache = [];
try
    cache = dataCache('create');
    cache = dataCache('add',cache,'key1',1);
    cache = dataCache('add',cache,'key2',2);
    cache = dataCache('add',cache,'key3',3);
catch exception
    fail = 1;
    disp('dataCache_test: failed rebuild 6-7 test');
    exception.getReport
end


% test delete (7), the middle item
try
    cache = dataCache('delete',cache,'key2');
catch exception
    fail = 1;
    disp('dataCache_test: failed delete7 test');
    exception.getReport
end

% check that the other items are there
if dataCache('isempty',cache)
    fail = 1;
    disp('dataCache_test: failed delete7b test');
end
if dataCache('length',cache) ~= 2
    fail = 1;
    disp('dataCache_test: failed delete7c test');
end
if isempty(dataCache('get',cache,'key3'))
    fail = 1;
    disp('dataCache_test: failed delete7d test');
end
if ~isempty(dataCache('get',cache,'key2'))
    fail = 1;
    disp('dataCache_test: failed delete7e test');
end
if isempty(dataCache('get',cache,'key1'))
    fail = 1;
    disp('dataCache_test: failed delete7f test');
end
if dataCache('get',cache,'key3') ~= 3
    fail = 1;
    disp('dataCache_test: failed delete7g test');
end
if dataCache('get',cache,'key1') ~= 1
    fail = 1;
    disp('dataCache_test: failed delete7g test');
end
