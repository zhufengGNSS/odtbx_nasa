function failed = staEarthTide_test(type)
%
% staEarthTide_test Regression test for staEarthTide
% See also: staEarthTide.m
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
%   Ravi Mathur         08/28/2012      Extracted from staEarthTide.m

if strcmpi(type,'ValidationTest')

	failed = staEarthTide_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = staEarthTide_regression_test();
else
    disp('staEarthTide_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = staEarthTide_validation_test();

% STAEARTHTIDE_VALIDATION_TEST - Validation Test for staEarthTide

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/11/2008	 	Original
%   Ravi Mathur         08/28/2012      Extracted from staEarthTide.m

% The portion of this algorithm corresponding to equations 4.18 and 4.22 in Will Campbell's
% thesis were verified against the thesis.

n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 1 % ECEF Input
Sta(Case).staPos = JATConstant('rEarth') * [1/sqrt(2) 1/sqrt(2) 0.01]';
Sta(Case).Epoch = datenum('August 1, 2003');
[AdStaPos(:,1),AdLat(1),AdLon(1),AdHeight(1)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2 % LLA Input
Sta(Case).lat = 25*pi/180;
Sta(Case).lon = 318*pi/180;
Sta(Case).height = 115;
Sta(Case).StationInfo = 'LatLonHeight';
Sta(Case).Epoch = datenum('Jan 1, 2000');
[AdStaPos(:,2),AdLat(2),AdLon(2),AdHeight(2)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 3 % Multiple ECEF Inputs
Sta(Case).staPos = [Sta(1).staPos,LLA2ecef(Sta(2).lat*180/pi,Sta(2).lon*180/pi,...
                    Sta(2).height/1e3)*1e3];
Sta(Case).Epoch = [Sta(1).Epoch,Sta(2).Epoch];
Sta(Case).StationInfo = 'ECEF';
[AdStaPos(:,3:4),AdLat(3:4),AdLon(3:4),AdHeight(3:4)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 4 % Multiple LLH Inputs
[lat,lon,height] = ecef2LLA(Sta(1).staPos/1e3);
Sta(Case).lat = [lat*pi/180,Sta(2).lat];
Sta(Case).lon = [lon*pi/180,Sta(2).lon];
Sta(Case).height = [height*1e3,Sta(2).height]; 
Sta(Case).StationInfo = 'LatLonHeight';
Sta(Case).Epoch = [Sta(1).Epoch,Sta(2).Epoch];
[AdStaPos(:,5:6),AdLat(5:6),AdLon(5:6),AdHeight(5:6)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 5 % Multiple ECEF Inputs - single epoch
Sta(Case).staPos = Sta(3).staPos;
Sta(Case).Epoch = Sta(1).Epoch;
Sta(Case).StationInfo = 'ECEF';
[AdStaPos(:,7:8),AdLat(7:8),AdLon(7:8),AdHeight(7:8)] = staEarthTide(Sta(Case))

ResultsMatrix = [AdStaPos;AdLat;AdLon;AdHeight];

failed = 0;
tol = 1e-4;

if exist('staEarthTide_ValidationData8_08.mat')

	%%%%%%%%%%%%%%%%%%%%%%%%%
	% The following lines were commented out after being used to generate validation data

% 	currentdir = pwd;
% 	cd /matlab/dsn
% 	format long
% 
% 	global drift_flag we PM_data
% 	pm_reader; % load PM_data
% 	we = JATConstant('wEarth');
% 	drift_flag = [0 1 0 0];
%
% 	Epoch = matlabTime2JD(Sta(1).Epoch);
% 	r1 = stn_drift(Sta(1).staPos'/1e3,[],[],[],Epoch)'*1e3
% 
% 	Epoch = matlabTime2JD(Sta(2).Epoch);
% 	r2 = stn_drift(Sta(3).staPos(:,2)'/1e3,[],[],[],Epoch)'*1e3
% 
% 	Epoch = matlabTime2JD(Sta(1).Epoch);
% 	r3 = stn_drift(Sta(3).staPos(:,2)'/1e3,[],[],[],Epoch)'*1e3
% 
%  	cd(currentdir)
% 
% 	save staEarthTide_ValidationData6_09 r1 r2 r3

	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	% Load validation values for AzEx and ElEx
	load staEarthTide_ValidationData8_08

	[lat,lon,height] = ecef2LLA(r1/1e3);
	lat1 = lat*pi/180;
	lon1 = lon*pi/180;
	height1 = height*1e3;

	[lat,lon,height] = ecef2LLA(r2/1e3);
	lat2 = lat*pi/180;
	lon2 = lon*pi/180;
	height2 = height*1e3;

	[lat,lon,height] = ecef2LLA(r3/1e3);
	lat3 = lat*pi/180;
	lon3 = lon*pi/180;
	height3 = height*1e3;

	ExpectedMatrix = [r1,r2,r1,r2,r1,r2,r1,r3
	                  lat1,lat2,lat1,lat2,lat1,lat2,lat1,lat3
	                  lon1,lon2,lon1,lon2,lon1,lon2,lon1,lon3
	                  height1,height2,height1,height2,height1,height2,height1,height3];

	Diff = [ExpectedMatrix - ResultsMatrix]

	if any( any( abs(Diff) > tol )) | any(isnan(Diff))
		failed = 1;
	end
else
	failed = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = staEarthTide_regression_test();

% STAEARTHTIDE_REGRESSION_TEST - Regression Test for staEarthTide

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         08/04/2008	 	Original
%   Ravi Mathur         08/28/2012      Extracted from staEarthTide.m

n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 1 % ECEF Input
Sta(Case).staPos = JATConstant('rEarth') * [1/sqrt(2) 1/sqrt(2) 0.01]';
Sta(Case).Epoch = datenum('August 1, 2003');
[AdStaPos(:,1),AdLat(1),AdLon(1),AdHeight(1)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2 % LLA Input
Sta(Case).lat = 25*pi/180;
Sta(Case).lon = 318*pi/180;
Sta(Case).height = 115;
Sta(Case).StationInfo = 'LatLonHeight';
Sta(Case).Epoch = datenum('Jan 1, 2000');
[AdStaPos(:,2),AdLat(2),AdLon(2),AdHeight(2)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 3 % Multiple ECEF Inputs
Sta(Case).staPos = [Sta(1).staPos,LLA2ecef(Sta(2).lat*180/pi,Sta(2).lon*180/pi,...
                    Sta(2).height/1e3)*1e3];
Sta(Case).Epoch = [Sta(1).Epoch,Sta(2).Epoch];
Sta(Case).StationInfo = 'ECEF';
[AdStaPos(:,3:4),AdLat(3:4),AdLon(3:4),AdHeight(3:4)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 4 % Multiple LLH Inputs
[lat,lon,height] = ecef2LLA(Sta(1).staPos/1e3);
Sta(Case).lat = [lat*pi/180,Sta(2).lat];
Sta(Case).lon = [lon*pi/180,Sta(2).lon];
Sta(Case).height = [height*1e3,Sta(2).height]; 
Sta(Case).StationInfo = 'LatLonHeight';
Sta(Case).Epoch = [Sta(1).Epoch,Sta(2).Epoch];
[AdStaPos(:,5:6),AdLat(5:6),AdLon(5:6),AdHeight(5:6)] = staEarthTide(Sta(Case));

%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 5 % Multiple ECEF Inputs - single epoch
Sta(Case).staPos = Sta(3).staPos;
Sta(Case).Epoch = Sta(1).Epoch;
Sta(Case).StationInfo = 'ECEF';
[AdStaPos(:,7:8),AdLat(7:8),AdLon(7:8),AdHeight(7:8)] = staEarthTide(Sta(Case))

ResultsMatrix = [AdStaPos;AdLat;AdLon;AdHeight];

failed = 0;
tol = 1e-4;

% PreviousResultsMatrix = ResultsMatrix;
% save staEarthTide_RegressionData6_09 PreviousResultsMatrix

if exist('staEarthTide_RegressionData6_09.mat')

	disp(' ')
	disp('Performing Regression Tests...')
	disp(' ')

	% Load validation values for AzEx and ElEx
	load staEarthTide_RegressionData6_09

	Diff = [PreviousResultsMatrix - ResultsMatrix]

	if any( any( abs(Diff) > tol )) | any(isnan(Diff))
		failed = 1;
	end
else
	failed = 1;
end

end