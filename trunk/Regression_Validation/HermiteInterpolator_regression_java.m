function [failed,results]=HermiteInterpolator_regression_java

% Regression unit test for interpolation in java HermiteInterpolator
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

% Author Stephen D Metcalfe, Emergent Space Technologies
% Date 9-29-2010

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Stephen D Metcalfe  09/29/2010   	Original
%   Stephen D Metcalfe  10/04/2010      Fixed MLINs.
%                                       Only outputs result file when asked.

WRITERESULTS = 0;       % set to non-zero to write the results file

% preallocate vaiables in case of total failure
test.neqns      = 0;
test.dimensions = 0;
test.entries    = 0;
test.start      = 0.0;
test.step       = 0.0;
test.current    = 0.0;
test.failed     = 0;
resultset(1)    = test;
results         = ones(1);
failed          = -1;

% Run java regression test 
rtn = jat.matlabInterface.unittest.HermiteInterpolator_regression.run_regression();
[m n] = size(rtn);

% check for total failure due to an exception
if n >= 8
    % preallocate dataspace 
    count = sum(fix(rtn(:,1)));
    resultset(1:count) = test;
    results(1:count)   = 0;
    
    failed  = 0;
    testnum = 1;
    for t = 1:m
        expected = fix(rtn(t,1));   % get number of tests expected
        neqns    = fix(rtn(t,2));   % get number of equations
        if WRITERESULTS ~= 0 % save the test results for detailed analysis
            % format resultset similar to HermiteInterpolator_regression.m
            test.neqns      = neqns;            % number of equations
            test.dimensions = fix(rtn(t,3));    % number of dimensions
            test.entries    = fix(rtn(t,4));    % number of entries
            test.start      = rtn(t,5);         % start time
            test.step       = rtn(t,6);         % interval size
            test.current    = rtn(t,7);         % current time
        end
        
        for r = 1:expected
            test.failed = fix(rtn(t,7+r));
            
            % check for test failure due to an exception
            if(expected ~= neqns)
                test.failed = 1;
            end
            
            % count the failures
            if test.failed ~= 0
                failed = failed +1;
            end
            
            if WRITERESULTS ~= 0 % save the test results for detailed analysis
                % resultset is preallocated above to count items
                resultset(testnum) = test; %#ok<AGROW>
            end
            results(testnum)    = test.failed;
            testnum = testnum + 1;
        end
    end
end

if failed ~= 0
    warning('JAT:HermiteInterpolator:regression','%d failures.\n',failed);
end

if WRITERESULTS ~= 0 % save the test results for detailed analysis
    RUNTIME = datestr(now);
    save HermiteInterpolator_regression_java_results.mat RUNTIME failed results resultset;
end
