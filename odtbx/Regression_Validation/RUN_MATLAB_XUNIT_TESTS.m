%% RUN matlab.unittest.TestCase Test Cases
% NOTE: This must be run in a directory containing the unit tests
clear classes; clear all; close all; clc;

%% Create a TestRunner object
%runner = matlab.unittest.TestRunner.withNoPlugins;  % silent runner
runner = matlab.unittest.TestRunner.withTextOutput; % verbose runner
%% Create a test suite from the TestCase classes contained in the pwd
suite = matlab.unittest.TestSuite.fromFolder(pwd,'IncludingSubfolders',true);

if isempty(suite)
    error('No TestCase classes were found, are you in the correct directory?')
end

%% Run the Tests with Profiler Turned on
profile on -history

test_results = run(runner,suite);

p = profile('info');
profile off
% You can now go to the directory where the tested code is located, then go 
% to the arrow on the upper right of the Current Folder window.  Then go 
% to Reports->Coverage Report to view your code coverage.

%% To Run a Specific Test, see the code below

% Run a specific TestCase Clase
%example: test_result = run(test_tdrssmeas_basic)
%test_result = run(test_tdrssmeas_hifi)


% Run a specific test within a TestCase Clase
%example: test_result = run(test_tdrssmeas_basic,'testBadTdrsType')
