% Covariance Comparison Function
% This function was created to compare the baseline covariance against the
% current covariance values for the ODTBX regression testing.
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

% Created by Kathryn Gregory, Emergent Space Technologies
% Date Created: 25 April 2007

function [agreement] = compareCovariance(P, testcaseName, validationDirectory);

if( nargin < 3)
    validationDirectory = pwd;
end

% Initialize agreement
agreement = 0;

% Load the Initial and Final Baseline Covariance Matrices
fileName = strcat('covariance_baseline_',testcaseName,'.mat');
baselineCovariance=load(fullfile(validationDirectory, fileName));

if(~isnan(P(1,end)))
    [x,y] = size(P);
    theEnd = y;
else
    theEnd = min(find(isnan(P(1,:))))-1;
end

% ODTBX Covariance Matrices at first and last times
ODTBX_Covariance_Time_Zero = (unscrunch(P(:,1)))*10^6;
ODTBX_Covariance_Time_Last = (unscrunch(P(:,theEnd)))*10^6;

% Covariance difference for first and last times
Covariance_Difference_Time_Zero = ODTBX_Covariance_Time_Zero-baselineCovariance.ODTBX_Covariance_Time_Zero;
Covariance_Difference_Time_Last = ODTBX_Covariance_Time_Last-baselineCovariance.ODTBX_Covariance_Time_Last;
timeZeroDiff = max(abs(Covariance_Difference_Time_Zero));
timeLastDiff = max(abs(Covariance_Difference_Time_Last));

% Check if ODEAS Covariance and current run ODTBX Covariance still agree and
% Write out difference file if they don't. 
if((timeZeroDiff >= 1e-8) | (timeLastDiff >= 1e-8))
    file = '_covdiff.txt';
    fileName = strcat(date,testcaseName,file);
    wholefile = fullfile(validationDirectory, fileName);                
    fid3 = fopen(wholefile, 'w');
    fprintf(fid3,'This file contains the covariance matrices for ODTBX Current and ODTBX Baseline:\n 1. ODTBX Covariance at First time step \n 2. ODTBX Baseline Covariance at First time step \n 3. ODTBX Covariance at Last time step \n 4. ODTBX Baseline Covariance at Last time step \n');
    fclose(fid3);
    dlmwrite(fileName, ODTBX_Covariance_Time_Zero,'-append','roffset', 1, ...
        'delimiter', '\t', 'precision', 6);
    dlmwrite(fileName, baselineCovariance.ODTBX_Covariance_Time_Zero,'-append','roffset', 1, ...
        'delimiter', '\t', 'precision', 6);
    dlmwrite(fileName, ODTBX_Covariance_Time_Last,'-append','roffset', 1, ...
        'delimiter', '\t', 'precision', 6);
    dlmwrite(fileName, baselineCovariance.ODTBX_Covariance_Time_Last,'-append','roffset', 1, ...
        'delimiter', '\t', 'precision', 6);
    agreement = 1; % There is a disagreement so the flag is set to indicate the disagreement
end