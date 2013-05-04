function failed = estval_test
%
% estval_test Regression test wrapper for estval
% See also: estval.m
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%
%   Ravi Mathur         08/28/2012      Created
%   Ravi Mathur         05/03/2013      Fully extracted regression test
%                                       from original estval.m

totaltests = 2;

% Run all the test cases
for k = 1:totaltests,
    disp(['Case ',num2str(k),'...'])
    fail(k) = run_test(k);
end

failed = any(fail);

end %function

% Runs the actual estval tests. The contents of this function have been
% extracted from estval.
function fail = run_test(testnum)

% Common test parameters
n = 19;
N = 100;
K = 10;
randn('state', 0);
e=randn(n,N,K);
fs=e+1.5;
t = linspace(0,100,N);
Pt = [];

% This test will be opening several figures, and we don't want them to
% overwrite any existing figures. So find the largest open figure handle.
% There may be figures with non-integer handles, e.g. if the 'publish'
% command was recently used. We need to ignore these.
allfigs = [get(0, 'children');0]; % Get all open figures
badidx = find(allfigs ~= floor(allfigs)); % Find non-integer handles
allfigs(badidx) = []; % Remove non-integer handles
fhs = max(allfigs); % Get largest of remaining integer handles

% Set test-dependent parameters
load estval_test;
switch testnum
    case 1
        iusemn = 0;
        g_save = g_save0;
        
    case 2
        iusemn = 1;
        g_save = g_save1;
end

% Run estval test, getting plot info as the result
h = estval(t, e, fs, Pt, fhs, iusemn);

% Check results
lh = length(h);
fail = 0;
for i = lh:-1:1
    g(i) = get(h(i));
    fail = any([fail ( any (abs(g(i).YData - g_save(i).YData) > eps) )]);
end

% Close generated figures
close (fhs+1:gcf);

end