function regressionTesting(varargin)

% REGRESSIONTESTING - ODTBX Regression Test Driver
%
% This script runs a pre-defined set of test cases on the ODTBX codebase.
% The tests can be external or built-in tests.  The results (pass, fail,
% or broken) are collected and written to a hard-coded log file.  
%
% Additional information can be appended to these results, if desired.
%
% This script can optionally email the results to a designated email
% recipient and SMTP server.
%
% Usage
%
%    local machine - "regressionTesting"
%
%    auto regression:
%       "regressingTesting(logpath, email, smptserver)"
%           where:
%           logpath = (optional, recommended) the path where the log file
%                     will be placed, will write to the screen if not
%                     supplied
%           email = (optional) the email of the person to receive the 
%                    report, leave empty ('') if specifying other arguments
%           smptserver = (optional) the fully-qualified name of the SMTP
%                        server to send the email, leave empty ('') if
%                        specifying other arguments
%           infofile = (optional) the full filepath of any additional text
%                        information to append to the log file
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

% Author: Kate Gregory, Emergent Space Technologies
% Constructed for the Orbit Determination Toolbox (ODTBX) at NASA Goddard.
% Date Completed: 25 April 2007
% Modified: 20 September 2007 Derek Surka - added R2.0 tests
% Modified: March 2008  Mike Southward - corrected autoregression
% Modified: 6 April 2008 Keith Speckman - Modified for use on local and cypher
% Modified 16 Dec 2008 Allen Brown - documented new arguments
% Modified 22 Dec 2008 Kevin Berry - replaced old measurement subfunction tests
% Modified 28 Jan 2009 Allen Brown - added getIERSTimes_test
% Modified 10 Feb 2009 Allen Brown - renamed jatRK*_tests
% Modified Mar 2009 Allen Brown - added estbat & estseq tests
%                                 added infofile output capability
%                                 added testing with and without toolboxes
% Modified 16 Jun 2009 Kevin Berry - added tdrssmeas test
%                                    added lnrmeas test
% Modified 10 Nov 2009 Kevin Berry - added ddormeas test
% Modified 29 Sept 2010 Stephen D Metcalfe - added HermiteInterpolator_regression test
%                                            added HermiteInterpolator_regression_java test
% Modified 15 Oct 2010 Allen Brown - Reorganized by release (again)
% Modified Mar-Apr 2011 Allen Brown - Added Release 4.5 tests.

% Regression Test Cases are described in the headers of the individual test case files.

close all
clc

% Cell array list of test case ODTBX MATLAB function names and any 
% required function arguments.  They are run with one output argument,
% as in "failed = eval(testCases{i})".  Each is required to return 0
% if the case passed, and a 1 if the case failed.
testCases = {
    % Release 5.0
        'gasdyn_test'
        'opt_sched_test'
        'initialize_cov_test'
        'estpar_test'
        'sensmos_test'
        'fixzoom'
        'gpspseudomeas_test'
        'clk_all_test'
        'nbody_test'
        'restartrecord_test(''SRIF'')'
        'polydyn_test'
        'integ_test'
        'integev_test'
    % Release 4.5
        'compValsTol_test'
        'FilterGpsCnoThresh_test'
        'FilterGpsExplicitTimes_test'
        'FilterGpsTimeWindow_test'
        'FilterGpsPrn_test'
        'FilterGpsBlock_test'
        'FilterGpsMulti_test'
        'FilterPpEl_test'
        'FilterPpExplicitTimes_test'
        'FilterEstCnoThresh_test'
        'gpslinkbudget_test'
        'gps_est_cno_test'
        'gps_gain_test'
        'slerp_test'
        'gps_phys_params_test'
        'RinexOReader_test'
        'rinexo2gpsdata_test'
        'gpsenh_vs_gpsmeas_test'
        'fgprop_test'
        'opnavmeas_test'
        'opnavmeas_test(''CCD'')'
    % Release 4.0
        'HermiteInterpolator_regression'
        'HermiteInterpolator_regression_java'
        'xnavmeas(''RegressionTest'')'
        'testMeasPartials_regression'
        'restartrecord_test'
        'test_comp_bs_3d_56'
    % Release 3.5
        'gsCoverage(''RegressionTest'')'
        'estsrif' % run internal self-test
    % Release 3.0
        'LOSRange(''RegressionTest'')'
        'LOSRangeRate(''RegressionTest'')'
        'LOSDoppler(''RegressionTest'')'
        'lightTimeCorrection(''RegressionTest'')'
        'rrdotlt(''RegressionTest'')'
        'srpAccel(''RegressionTest'')'
        'ttdelay(''RegressionTest'')'
        'stationData(''RegressionTest'')'
        'staEarthTide(''RegressionTest'')'
        'EarthOrbitPlot(''RegressionTest'')'
        'relativityLightDelay(''RegressionTest'')'
        'tdrssmeas(''RegressionTest'')'
        'lnrmeas(''RegressionTest'')'
        'kepprop2b(''RegressionTest'')'
        'ddormeas(''RegressionTest'')'
        'kep2cart_test'
        'dataCache_test'
        'estspf_regression'
        'varpiles'
        'covmake_test'
        'gpsmeas_test_driver'
        'estbat_no_proc_noise_test'
        'estval'
        'editflag_test'
  %      'plot_ominusc'
  %      'montecarloseed_test'
        'updatevectorized_test'
        'estbat_tspan_test'
    % Release 2.0
        'gsmeas_test'
        'odtbxOptions_test'
        'rrdotang(''RegressionTest'')'
    % JAT_Adapters
        'convertTime_test'
        'createJATWorld_test'
        'ecef2LLA(''RegressionTest'')'
        'ephemDE405_test'
        'getGroundStationInfo_test'
        'getIndex_test'
        'JATConstant_test'
        'jatDCM_test'
        'jatHCW_test'
    %    'jatWorldPropagatorRK4_test'
    %    'jatWorldPropagatorRK8_test'
    %    'jatWorldPropagatorRKF78_test'
        'jatTropoModel(''RegressionTest'')'
        'matlabTimeJDMJD_test'
        'LLA2ecef(''RegressionTest'')'
        'jatStaAzEl(''RegressionTest'')'
        'jatIonoDelayModel(''RegressionTest'')'
        'jatIonoDelayModel_test' % independent validation check
        'jatChargedParticleModel(''RegressionTest'')'
        'getIERSTimes_test'
        'kepel_test'
        'jatWorldPropagatorRKF78_regression'
        'jatRK8_regression'
        'jatRK8_time_io_test'
        'jatForces_regression'
    % Release 1.0
    %     'estbat_001'
    %     'estbat_002'
    %     'estbat_003'
    %     'estseq_001'
    %     'estseq_002'
    %     'estseq_003'
    %     'estseq_004'
    %     'estseq_005'
    %     'estseq_006'
        'estbat' % run internal self-test
        'estseq' % run internal self-test
    };

brokenTests = {}; % cell array of the names and run name of all broken tests
                  % in {:,1} and stack dumps in {:,2}
failedTests = {}; % cell array of the names and run name of all failed tests
passedTests = {}; % cell array of the names and run name of all passing tests
results = 0; % the count of all broken or failed test cases for all runs

% Set up the log file for writing:
if nargin
    logPath = varargin{1};
    file = 'ODTBX_Pass_Fail_Log_File.txt'; % the hardcoded log file name
    fileName = [date, file];
    wholeFile = fullfile(logPath,fileName);
    fid1 = fopen(wholeFile,'w');
else
    fid1 = 1;
end

runNames = { '' 'No toolboxes' };

for run = 1:length(runNames)
    
    if run == 2
        removeToolboxPaths(); % remove all optional toolboxes
    end

    % Run all the tests in the default environment:
    [runresults, runpassedTests, runfailedTests, runbrokenTests] = runAllCases(testCases,runNames{run});

    % collect the results for all runs:
    results = results + runresults;
    if ~isempty(runpassedTests)
        for i = 1:length(runpassedTests)
            passedTests{end+1} = runpassedTests{i};
        end
    end
    if ~isempty(runfailedTests)
        for i = 1:length(runfailedTests)
            failedTests{end+1} = runfailedTests{i};
        end
    end
    if ~isempty(runbrokenTests)
        for i = 1:size(runbrokenTests,1)
            brokenTests{end+1,1} = runbrokenTests{i,1};
            brokenTests{end,2} = runbrokenTests{i,2};
        end
    end
end % run

% Write the log file with the Failure Names (if any)
if(exist('failedTests','var') && ~isempty(failedTests));
    fprintf(fid1,'The Following Test Cases Failed the Regression Testing:\n--------------------------------------------------------\n');
    for i = 1:length(failedTests)
        fprintf(fid1,'%s\n', failedTests{i});
    end
else
    fprintf(fid1,'\nNo Test Cases Failed the Regression Testing!\n\n');
end

% Write the log file with the Broken Case Details (if any)
if(exist('brokenTests','var') && ~isempty(brokenTests));
    fprintf(fid1,'\nThe Following Test Cases Returned Errors During the Regression Testing:\n-----------------------------------------------------------------------\n');
    for i = 1:size(brokenTests,1)
        fprintf(fid1,'\n\t%s\n', brokenTests{i,1});        fprintf(fid1,'Error in\n');        tab = '\t';
        for j = length(brokenTests{i,2}.stack):-1:1
            fprintf(fid1,[tab '==> %s at %i\n'],brokenTests{i,2}.stack(j).name,brokenTests{i,2}.stack(j).line); 
            tab = [tab '\t'];
        end
        fprintf(fid1,'Error Description:\n');
        fprintf(fid1,'%s\n\n',brokenTests{i,2}.message);
    end
else
    fprintf(fid1,'\n\nNo Test Cases Returned Errors During the Regression Testing!\n\n');
end

% Write the log file with the Passing Case Names
if(exist('passedTests','var') && ~isempty(passedTests));
    fprintf(fid1,'\nThe Following Test Cases Passed the Regression Testing:\n---------------------------------------------------------\n');
    for i = 1:length(passedTests)
        fprintf(fid1,'%s\n', passedTests{i});
    end
else
    fprintf(fid1,'No Test Cases Passed the Regression Testing!\n\n');
end

% If we have an infofile, then concatenate its contents to the log file.
if nargin >= 4 && ~isempty(varargin{4})
    finfo=fopen(varargin{4},'rt'); % text read on PC & unix
    if finfo < 0
        fprintf(fid1,'\n\nScript Error: regressionTestin.m failed to open given info file: %s',varargin{4});
    else
        % copy the contents of finfo into fid1
        fprintf(fid1,'\n\n---------------------------------------------------------\n');
        fprintf(fid1,'Additional Information (from: %s)\n',varargin{4});
        fprintf(fid1,'---------------------------------------------------------\n');
        while 1
            tline = fgetl(finfo);
            if ~ischar(tline), break, end
            fprintf(fid1,'%s\n',tline);
        end
        fclose(finfo);
    end
end

fclose all;

% If needed, Email out results if a test case failed or all test cases passed.
% Attach log message with results of passed and failed regression tests.
if nargin > 1
    email = varargin{2};
    if nargin > 2
        smtpServer = varargin{3};
    else
        smtpServer = 'mailhost.gsfc.nasa.gov';
    end
else
    email = '';
    smtpServer = '';
end
if ~isempty(email) && ~isempty(smtpServer)
    % Set up email preferences for all the function emails.
    setpref('Internet','SMTP_Server',smtpServer);
    setpref('Internet','E_mail','odtbx-builder@emergentspace.com');
    if (results > 0)
        message = sprintf('%d ODTBX regression test case(s) failed/broke. Please refer to the attached log message for details.',results);
        subject = strcat('An ODTBX Regression Test Case Failed/Broke on: ',date);
        sendmail(email,subject, message, {wholeFile});
    else
        message = {strcat('All of the ODTBX regression test cases passed. Refer to the attached log message for details.')};
        subject = strcat('All ODTBX Regression Test Cases Passed on: ',date);
        sendmail(email,subject, message, {wholeFile});
    end
end

%
% done.
%

function [results, passedTests, failedTests, brokenTests] = runAllCases(testCases, runName)
% Support function for running all the test cases and classifying them.
% This can be named via runName so that all cases can be re-run with
% various external differences (if required).
%
% Inputs:
%   testCases = nx1 cell array of funtions and arguments to test
%   runName = (optional) identifier string of the run's name
%
% Outputs:
%   results = 1x1 count of failed and broken cases
%   passedTests = nx1 cell array with the names of passing cases 
%                (from testCases)
%   failedTests = mx1 cell array with the names of failing cases 
%                (from testCases)
%   brokenTests = jx2 cell array with the names of broken cases 
%                ({:,1} from testCases, {:,2} struct from lasterror() )

if isempty(runName)
    localrunName = '';
else
    localrunName = strcat('(',runName,')');
end

results = 0;
passedTests = {};
failedTests = {};
brokenTests = {};

for i=1:length(testCases)
    
    try
        disp(sprintf('\nRunning %s.m %s...',testCases{i},localrunName));
        failed = eval(testCases{i});
        disp(sprintf('Done!\n'));

        % check that the test results conform to the interface
        if length(failed) > 1
            error('%s %s returned an incorrect result to the regression harness: length()=%d',...
                testCases{i},localrunName,length(failed));
        end

        if ~((failed == 1) || (failed == 0))
            error('%s %s returned an incorrect result to the regression harness: failed=%d',...
                testCases{i},localrunName,failed);
        end
        
        results = results + failed;
        if (failed);
            failedTests{end+1} = strcat(testCases{i},localrunName);
            disp(sprintf('Failed!\n'));
        else
            passedTests{end+1} = strcat(testCases{i},localrunName);
            disp(sprintf('Passed!\n'));
        end

    catch
        results = results + 1;
        brokenTests{end+1,1} = strcat(testCases{i},localrunName);
        brokenTests{end,2} = lasterror;
    end
    
end


function removeToolboxPaths()
% This function removes any optional toolbox paths from the MATLAB path
% and should restore MATLAB behavior to a non-toolbox environment.

[toolboxPaths] = findToolboxPaths();

for i = 1:length(toolboxPaths)
    rmpath(toolboxPaths{i});
end
