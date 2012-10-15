function [failed,results] = HermiteInterpolator_test

% Regression unit test for interpolation in HermiteInterpolator
% Uses test datasets for HermiteInterpolator verification
% Test datasets are created by HermiteInterpolator_test_data.m
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
% Date 6-13-2010

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Stephen D Metcalfe  06/13/2010   	Original
%   Stephen D Metcalfe  09/29/2010      Added revision history block.
%                                       Modified to explicitly reference 
%                                       test data in the DataFiles folder.
%                                       Added WHOGEN and GENTIME to results
%                                       to identify test data used.
%                                       Added timestamp to results to 
%                                       identify when test was run.
%   Stephen D Metcalfe  10/04/2010      Fixed MLINs and preallocated structs.
%                                       Only outputs result file when asked.
%                                       Added default values that should be
%                                       in the test data.
%   Ravi Mathur         08/27/2012      Rename to conform to new
%                                       regression test format

WRITERESULTS = 0;       % set to non-zero to write the results file
GENTIME = datestr(now); % default time generated if not in load
WHOGEN  = 'Unknown';    % default name of who generated if not in load
DETAIL  = 0;            % default flag if not in load

failed = 0;

% Setup the test parameters for HermiteInterpolator test
load DataFiles/HermiteInterpolator_test_data

% preallocate test results only if saving
if WRITERESULTS ~= 0
    test = testset(1);      % get test struct from testset
    test.neqns   = 0;       % number of equations
    test.start   = 0.0;     % start time
    test.step    = 0.0;     % time interval
    test.entries = 0;       % number of entries
    test.current = 0;       % current time to test
    test.result  = zeros(6);% interpolated result
    test.error   = zeros(6);% error
    test.percent = zeros(6);% percentage error
    test.failed  = 0;       % pass/fail
    resultset(1:length(testset)) = test;
end

results = zeros(length(testset),1);
for testnum=1:length(testset)
    test = testset(testnum);

    % display test parameters
    if DETAIL ~= 0
        fprintf('%-5d start=%f step=%f entries=%d currentTime=%f - ',...
            testnum,test.start,test.step,test.entries,test.current);
    end
    
    % run the interpolator
    test.result = HermiteInterpolator(test.x,test.current,test.neqns,test.f);            
    
    % examine results for errors
    test.error   = test.expected - test.result;
    test.percent = (test.error./test.expected)*100;
    
    % failed if result out of bounds
    if abs(test.percent) > PERR
        test.failed = 1;
        failed = failed+1;
    else
        test.failed = 0;
    end
    
    results(testnum) = test.failed;
    if WRITERESULTS ~= 0    % save results only if detailed analysis is required
        resultset(testnum) = test;
    end

    % display detailed status
    if DETAIL ~= 0 
        if test.failed ~= 0 
            fprintf('Failed\n'); 
        else
            fprintf('Passed\n'); 
        end
    end
    
    % abort test if failed and ENDONFAILURE flag is set
    if test.failed ~= 0 && ENDONFAILURE ~= 0 
        break;
    end
end

if failed ~= 0
    % GENTIME is loaded from DataFiles/HermiteInterpolator_test_data
    warning('ODTBX:HermiteInterpolator:test', 'Testset %s ... %d failures.\n', GENTIME, failed); 
end

if WRITERESULTS ~= 0
    RUNTIME = datestr(now);
    % WHOGEN and GENTIME are loaded from DataFiles/HermiteInterpolator_test_data
    save HermiteInterpolator_test_results.mat WHOGEN GENTIME RUNTIME failed results resultset; 
end
