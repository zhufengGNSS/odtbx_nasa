function failed = getIERSTimes_test()

% Test for JAT Adaptor getIERSTimes.
% This test is built upon the results from the JAT EOP.dat datafile found
% in the jat.spacetime classpath.  If the data in that file changes then
% this test should be updated accordingly.
%
%   OUTPUTS:
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%   failed      1x1         1 indicates test failure, 0 indicates success
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

% NOTE: There are two possible "valid" results:
% 1) If a FitIERS has not been created in the JVM then its cached data
% won't have been initialized and we'll get zeros.
% 2) If a FitIERS has been created in the JVM then we'll get the answers
% we expect.

failed = 1; %default to failed

% a = the earliest valid time in EOP.dat
% b = the latest valid time in EOP.dat
[a, b] = getIERSTimes();

if (a + b) == 0.0
    % FitIERS has not been initialized, so let's do that now and re-test
    o = jat.spacetime.FitIERS();
    
    % re-get a & b
    [a, b] = getIERSTimes();
end

% floating point numerical comparison tolerance:
% (these numbers should be pretty exact since there's no calculation)
tol = 1e-15;

% see EOP.dat and the FitIERS unit test (FitIERSTest.java)
a_test = 41684.0;
b_test = 54054.0;

if (abs(a - a_test) < tol)
    if (abs(b - b_test) < tol)
        failed = 0; % passed
    else
        disp(sprintf(...
            'TEST FAILURE: getIERSTimes_test: latest valid time returned %f, but expected %f.',...
            b, b_test));
    end
else
    disp(sprintf(...
        'TEST FAILURE: getIERSTimes_test: earliest valid time returned %f, but expected %f.',...
        a, a_test));
end
    
    