function failed = stationData_test(type)
%
% stationData_test Regression test for stationData
% See also: stationData.m
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
%   Ravi Mathur         08/28/2012      Extracted from stationData.m

if strcmpi(type,'ValidationTest')

	failed = stationData_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = stationData_regression_test();
    
else
    disp('stationData_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stationData_validation_test - Validation for stationData.m

function failed = stationData_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/10/2008	 	Original

disp(' ')
disp('Performing Validation Test....')
disp(' ')
disp('This test passes if there is no error.')
disp(' ')

Sta = stationData([12 27 43 66])

Sta = stationData([23 26 61],Sta)

Sta = stationData

failed = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stationData_regression_test - Regression for stationData.m

function failed = stationData_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/10/2008	 	Original

disp(' ')
disp('Performing Regression Test....')
disp(' ')

tol = 1e-10;

Sta(1) = stationData;
Sta(2) = stationData([12 13 14 15 16 17 23 24 25 26 27 28 33 34 42 43 45 46 53 54 61 63 65 66]);
Sta(3) = stationData([12 13 14 15 16 17 23 24 25 26 27 28 33 34 42 43 45 46 53 54 61 63 65 66],...
                     Sta(2))

StaResults = [];
for k = 1:3
	StaResults = [StaResults;Sta(k).DSS_IDs;Sta(k).siteIDs;Sta(k).lat;Sta(k).lon;...
	              Sta(k).height;Sta(k).staPos;Sta(k).posUnc;Sta(k).Press;Sta(k).Temp;...
		      Sta(k).RelHum;Sta(k).TempLapse];
end

%PreviousStaResults = StaResults;
%save stationData_RegressionData6_08 PreviousStaResults

failed = 0;

if exist('stationData_RegressionData6_08.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load the stationData outputs
	load stationData_RegressionData6_08

	Diff = PreviousStaResults - StaResults;

	if (any(any( abs(Diff) > tol )) |  any(any(isnan(Diff))))
		failed = 1;
	end

else
	failed = 1;
end


end