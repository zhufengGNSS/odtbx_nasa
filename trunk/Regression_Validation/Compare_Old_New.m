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

function [agreement, measdiff, trajdiff] = Compare_Old_New(plotResults, regressionTest, testcaseName, ValidationDirectory)

if( nargin < 4 )
    ValidationDirectory = pwd;
end

%Compare New and Old Measurements for Regression Testing
measFileNameOld = strcat('measurements_baseline_',testcaseName,'.meas');
OldMeasFile = fullfile(ValidationDirectory, measFileNameOld);

fid1 = fopen(OldMeasFile);
OldMeasMat = textscan(fid1,'%f%f%f','headerLines',2);
fclose(fid1);

measFileNameNew = 'ODTBX.meas';
NewMeasFile = fullfile(ValidationDirectory, measFileNameNew);

fid2 = fopen(NewMeasFile);
NewMeasMat = textscan(fid2,'%f%f%f','headerLines',2);
fclose(fid2);

%Compare New and Old Trajectories for Regression Testing
stateFileNameOld = strcat('state_baseline_',testcaseName,'.trj');
OldStateFile = fullfile(ValidationDirectory, stateFileNameOld);

fid3 = fopen(OldStateFile);
OldStateMat = textscan(fid3,'%f%f%f%f%f%f%f','headerLines',2);
fclose(fid3);

stateFileNameNew = 'ODTBX.trj';
NewStateFile = fullfile(ValidationDirectory, stateFileNameNew);

fid4 = fopen(NewStateFile);
NewStateMat = textscan(fid4,'%f%f%f%f%f%f%f','headerLines',2);
fclose(fid4);

%Compare Measurements
measTime = OldMeasMat{1};
meas_TimeDiff = OldMeasMat{1} - NewMeasMat{1};
RangeDiff = OldMeasMat{2} - NewMeasMat{2};
RangeRateDiff = OldMeasMat{3} - NewMeasMat{3};
measdiff1 = max(abs(RangeDiff));
measdiff2 = max(abs(RangeRateDiff));
measTimeDiff = max(abs(meas_TimeDiff));

%Compare Trajectories
trajTime = OldStateMat{1};
traj_TimeDiff = OldStateMat{1} - NewStateMat{1};
PosDiff = sqrt((OldStateMat{2}+OldStateMat{3}+OldStateMat{4}).^2) - sqrt((NewStateMat{2}+NewStateMat{3}+NewStateMat{4}).^2);
VelDiff = sqrt((OldStateMat{5}+OldStateMat{6}+OldStateMat{7}).^2) - sqrt((NewStateMat{5}+NewStateMat{6}+NewStateMat{7}).^2);
trajdiff1 = max(abs(PosDiff));
trajdiff2 = max(abs(VelDiff));
trajTimeDiff = max(abs(traj_TimeDiff));
agreement = 0;
measdiff = 0;
trajdiff = 0;
currentDir = pwd;

% Evaluate if the Code has changed by comparing against regression test
if(regressionTest == 1);
    if ((measdiff1 >= 1e-5) | (measdiff2 >= 1e-6) | (measTimeDiff >= 1e-10));
        file = '_measdiff.txt';
        fileName = strcat(date, testcaseName, file);
        wholefile = fullfile(currentDir, fileName);
        fid1 = fopen(wholefile, 'w');

        fprintf(fid1,'Time (s)\t     Range Difference(m)     RangeRate Difference(m/s)\n--------        ------------------           ------------------\n');
        for i=1:1:length(measTime);
            fprintf(fid1,'%3.4f %30.13f %23.13f\n', measTime(i), RangeDiff(i), RangeRateDiff(i));
        end
        fclose(fid1);
        agreement = 1;
        measdiff = 1;
    end

    if ((trajdiff1 >= 1e-5) | (trajdiff2 >= 1e-8) | (trajTimeDiff >= 1e-10));
        file = '_trajdiff.txt';
        fileName = strcat(date, testcaseName, file);
        wholefile = fullfile(currentDir, fileName);                
        fid = fopen(wholefile, 'w');

        fprintf(fid,'Time (s)\t     Position Difference(m)     Velocity Difference(m/s)\n--------         -------------------         -----------------\n');
        for i=1:1:length(trajTime);
            fprintf(fid,'%3.4f %30.13f %23.13f\n', trajTime(i), PosDiff(i), VelDiff(i));
        end
        fclose(fid);
        agreement = 1;
        trajdiff = 1;
    end
end

%Plot Measurement Differences over simulation time
if(plotResults==1)
    figure
    subplot(1,2,1)
    plot(OldMeasMat{1},RangeDiff,'r*')
    title('Measurement difference in Range over time in meters')
    xlabel('Time(s)')
    ylabel('Range(m)')

    subplot(1,2,2)
    plot(OldMeasMat{1},RangeRateDiff,'r*')
    title('Measurment difference in Range-Rate over time in meters/second')
    xlabel('Time(s)')
    ylabel('Range-Rate(m/s)')

%Plot Trajectory Differences over simulation time
    figure
    subplot(1,2,1)
    plot(OldStateMat{1}*60,PosDiff,'r*')
    title('Trajectory difference in Position over time in meters')
    xlabel('Time(s)')
    ylabel('Position(m)')

    subplot(1,2,2)
    plot(OldStateMat{1}*60,VelDiff,'r*')
    title('Trajectory difference in Velocity over time in meters/second')
    xlabel('Time(s)')
    ylabel('Velocity(m/s)')
end


