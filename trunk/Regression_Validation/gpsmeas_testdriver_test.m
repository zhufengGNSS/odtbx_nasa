%
% gpsmeas_testdriver_test
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

% Brent Wm. Barbee

% Begin function.
function [failed] = gpsmeas_testdriver_test(testcase, dump_results)
% gpsmeas_testdriver_test  Run a set of predefined test cases
%
% Can run a single test case or all test cases.  Will return whether any
% test cases failed and, if one failed, will display detailed result
% information.
%
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%   INPUTS
%    testcase       1x1     Optional.  The single test case to run.
%                           If not specified, will run all test cases.
%    dump_results   1x1     Optional.  If 1, will store the simulation
%                           results in a .mat file which can later be used
%                           as a truth file for the test case.  The created
%                           file will have the same name as the truth file
%                           with "_new" appended.
%                           If not specified, will not dump results.
%
%   OUTPUTS
%    failed         1x1     0 if all test cases passed.  
%                           1 if one of the test cases failed.

if nargin < 2
    dump_results = 0;
end

% Optional first argument, testcase, is checked later in the code

    % Load the truth data.
    truth_file{1} = 'RESULTS_leo_inertialpoint_1_gpsdef.mat';
      % LEO InertialPoint setup's results do not differ with antenna pointing
    truth_file{2} = 'RESULTS_mms_hemi_test_1_gpsdef.mat';
    truth_file{3} = 'RESULTS_geosat_highgainant_bwb_1_gpsdef.mat';
    truth_file{4} = 'RESULTS_geosat_highgainant_bwb_2_gpsdef.mat';
    truth_file{5} = 'RESULTS_mms_hemi_y2k_2_gpsdef.mat';
    truth_file{6} = 'RESULTS_leo_inertialpoint_3_gpsdef.mat';
    truth_file{7} = 'RESULTS_mms_hemi_test_3_gpsdef.mat';
    truth_file{8} = 'RESULTS_mms_hemi_y2k_3_gpsdef.mat';
    truth_file{9} = 'RESULTS_geosat_highgainant_bwb_4_gpsdef.mat';
    truth_file{10} = 'RESULTS_mms_hemi_y2k_4_gpsdef.mat';
    truth_file{11} = 'RESULTS_mms_hemi_test_56_gpsdef.mat';
    truth_file{12} = 'RESULTS_yaw_implimentation.mat';
    truth_file{13} = 'RESULTS_2D_user_and_GPS_antenna.mat';
    truth_file{14} = 'RESULTS_2D_user_1D_GPS_antenna.mat';
    truth_file{15} = 'RESULTS_leo_inertialpoint_1_gpsdef_omni.mat';

    % Set accuracy tolerance for pass/fail.
%     TOL(1) = 4.5e-3;
%     TOL(2) = 1.2e-2;
%     TOL(3) = 2.4e-2;
%     TOL(4) = 3.0e-2;
%     TOL(5) = 1.5e-2;
%     TOL(6) = 4.5e-3;
%     TOL(7) = 1.2e-2;
%     TOL(8) = 1.5e-2;
%     TOL(9) = 2.4e-2;
%     TOL(10) = 1.5e-2;
%     TOL(11) = 1.2e-2;
%     TOL(12) = 1.2e-2;
    TOL(1) = 4.5e-4;
    TOL(2) = 1.2e-4;
    TOL(3) = 2.4e-4;
    TOL(4) = 3.0e-4;
    TOL(5) = 1.5e-4;
    TOL(6) = 4.5e-4;
    TOL(7) = 1.2e-4;
    TOL(8) = 1.5e-4;
    TOL(9) = 2.4e-4;
    TOL(10) = 1.5e-4;
    TOL(11) = 1.2e-4;
    TOL(12) = 1.2e-4;
    TOL(13) = 1.2e-4;
    TOL(14) = 1.2e-4;
    TOL(15) = 4.5e-4;
    
    % Specify the appropriate Yuma file.
    yuma_file{1} = 'Yuma1134.txt';
    yuma_file{2} = 'Yuma1134.txt';
    yuma_file{3} = 'Yuma24-6_nominal.txt';
    yuma_file{4} = 'Yuma24-6_nominal.txt';
    yuma_file{5} = 'Yuma1042.txt';
    yuma_file{6} = 'Yuma1134.txt';
    yuma_file{7} = 'Yuma1134.txt';
    yuma_file{8} = 'Yuma1042.txt';
    yuma_file{9} = 'Yuma24-6_nominal.txt';
    yuma_file{10} = 'Yuma1042.txt';
    yuma_file{11} = 'Yuma1134.txt';
    yuma_file{12} = 'Yuma1134.txt';
    yuma_file{13} = 'Yuma1134.txt';
    yuma_file{14} = 'Yuma1134.txt';
    yuma_file{15} = 'Yuma1134.txt';

    % Specify the appropriate antenna pointing.
    ant_point(1).val = [1, -1];
    ant_point(2).val = [1, -1];
    ant_point(3).val = [-1];
    ant_point(4).val = [2, -2];
    ant_point(5).val = [2, -2];
    ant_point(6).val = [3, -3];
    ant_point(7).val = [3, -3];
    ant_point(8).val = [3, -3];
    ant_point(9).val = [4];
    ant_point(10).val = [4, -4];
    ant_point(11).val = [5, 6, -6];
    ant_point(12).val = [1];
    ant_point(13).val = [1];
    ant_point(14).val = [1];
    ant_point(15).val = [1]; % omni receive antenna

    % Specify the appropriate receive antenna pattern settings.
    ant_pat{1} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{2} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{3} = '{truth.rec_pattern1}';
    ant_pat{4} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{5} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{6} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{7} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{8} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{9} = '{truth.rec_pattern1}';
    ant_pat{10} = '{truth.rec_pattern1, truth.rec_pattern1}';
    ant_pat{11} = '{truth.rec_pattern1, truth.rec_pattern1,  truth.rec_pattern1}';
    ant_pat{12} = '{truth.rec_pattern1}';
    ant_pat{13} = '{truth.rec_pattern}'; % 3D antenna
    ant_pat{14} = '{truth.rec_pattern}'; % 3D antenna
    ant_pat{15} = '{truth.rec_pattern1}'; % omni receive antenna

if nargin < 1
    start_testcase=1;
    end_testcase=length(TOL);
    fprintf(1, '\n*** Executing %d test cases ...\n', end_testcase);
else
    start_testcase=testcase;
    end_testcase=testcase;
end

    failedTestCase = 0;
    % Loop over all test cases.
    for a=start_testcase:end_testcase

        if (~failedTestCase) 
            % Print a start message to the screen for the current test case.
            fprintf(1, '\n*** Executing gpsmeas test case %d ...\n', a);

            % Execute the current test case.
            failedTestCase = gpsmeas_test(truth_file{a}, TOL(a), ...
                yuma_file{a}, ant_point(a).val, ant_pat{a}, dump_results);
        end

    end
    [failed] = failedTestCase;
    
% End of function.
end
