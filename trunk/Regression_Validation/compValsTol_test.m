function fail = compValsTol_test()
% Unit test for compValsTol.
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
%
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

len = 100;

t = 1:len;% times to check
tol = 0.3;

% explicit times
% grab 10's, and 5, 6, 17, 19
checkvals = sort([ 10:10:len 5 6 17 19]);
%
% 5.3-5.7 are a series in between t's, the 5.3 & 5.7 match 5 & 6
% 10.5 is a standalone e that doesn't match
% 17.1 is a standalone e that matches left
% 18.9 is a standalone e that matches right
% the 22.5-22.6 are a series in between t's that don't match
e = sort([10:10:len 5.3 5.6 5.7 10.5 17.1 18.9 22.5 22.55 22.555 22.6]);

fprintf(1,'test_CompValsTol: Number of times to check: %d\n',len);
fprintf(1,'test_CompValsTol: Number of comparison times: %d\n',length(checkvals));


%% Test 1 - full data test
try
    tic
    [indA, matchA] = compValsTol(t,e,tol);
    toc
catch ex
    fprintf(1,'unexpected exception in test 1.\n');
    fail = 1;
end

if length(indA) ~= length(checkvals) || ~any(t(indA) == checkvals)
    fprintf(1,'bad indices in test 1.\n');
    fail = 1;
end

if length(matchA) ~= length(checkvals) || ~any(matchA == checkvals)
    fprintf(1,'bad values in test 1.\n');
    fail = 1;
end

%% Test 2 - full data test
try
    tic
    [indA, matchA] = compValsTol(t,e); % zero tol
    checkvals = sort([ 10:10:len]);
    toc
catch ex
    fprintf(1,'unexpected exception in test 2.\n');
    fail = 1;
end

if length(indA) ~= length(checkvals) || ~any(t(indA) == checkvals)
    fprintf(1,'bad indices in test 2.\n');
    fail = 1;
end

if length(matchA) ~= length(checkvals) || ~any(matchA == checkvals)
    fprintf(1,'bad values in test 2.\n');
    fail = 1;
end

%% Test 3 - empty arg
try
    tic
    [indA, matchA] = compValsTol(t,[],tol);
    toc
catch ex
    fprintf(1,'unexpected exception in test 3.\n');
    fail = 1;
end

if ~isempty(indA)
    fprintf(1,'unexpected indices in test 3.\n');
    fail = 1;
end

if ~isempty(matchA)
    fprintf(1,'unexpected indices in test 3.\n');
    fail = 1;
end


%% Test 4 - empty arg
try
    tic
    [indA, matchA] = compValsTol([],e);
    toc
catch ex
    fprintf(1,'unexpected exception in test 4.\n');
    fail = 1;
end

if ~isempty(indA)
    fprintf(1,'unexpected indices in test 4.\n');
    fail = 1;
end

if ~isempty(matchA)
    fprintf(1,'unexpected indices in test 4.\n');
    fail = 1;
end

