function [AdStaPos,AdLat,AdLon,AdHeight] = staEarthTide(Sta)

% STAEARTHTIDE Ground station position adjustment for solid earth tides
%
%  Like the oceans, tidal forces affect the land, though to a lesser extent.
%  Nevertheless, the displacement of tracking stations due to solid Earth tides
%  is sufficiently large to be included when determining a station's ECEF
%  position.  Tidal displacements may be as large as 50cm. (Will Campbell)
%
%  staEarthTide adjusts the station position to account for these tidal
%  effects.  The adjusted position is output in ECEF and lat, lon height.
%
%   [AdStaPos,AdLat,AdLon,AdHeight] = staEarthTide(Sta)
%
%  INPUTS 
%      VARIABLE            SIZE         DESCRIPTION (Optional/Default)
%      Sta.Epoch         (1xN)        UTC time in Matlab datenum format
%      Sta.StationInfo   string       Method of specifying ground station
%                                          position
%
%      CASE: Sta.StationInfo = 'ECEF'  (Default)
%      -------------------------------------------------------------------------
%      Sta.staPos        (3xM)        Station ECEF coordinates (m)
%
%      CASE: Sta.StationInfo = 'LatLonHeight'
%      -------------------------------------------------------------------------
%      Sta.lat           (1xM)        Station geodetic latitude (rad)
%      Sta.lon           (1xM)        Station geodetic longitude (rad)
%      Sta.height        (1xM)        Station geodetic height (m)
%
%  OUTPUTS 
%      AdStaPos          (3xM)        Adjusted Station ECEF coordinates (m)
%      AdLat             (1xM)        Adjusted Station geodetic latitude (rad)
%      AdLon             (1xM)        Adjusted Station geodetic longitude (rad)
%      AdHeight          (1xM)        Adjusted Station geodetic height (m)
%
% NOTES
%
%   This The algorithm in this function was extracted from the stn_drift.m
%   script developed by William Campbell (8.1.02).  The theory is explained in
%   his thesis.
%
%  VALIDATION TEST
%
%   To perform a validation test, replace the Sta structure with 
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
%  REGRESSION TEST
%
%   To perform a regression test, replace the Sta structure with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
%
%   keywords: ground station, earth tide
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman      06/11/2008	 	Original
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy and
%                                       updated the regression test

% Initialize LLH outputs
AdLat = 0;
AdLon = 0;
AdHeight = 0;

% Determine whether this is an actual call to the program or a test

if strcmpi(Sta,'ValidationTest')

	AdStaPos = staEarthTide_validation_test();

elseif strcmpi(Sta,'RegressionTest')

	AdStaPos = staEarthTide_regression_test();

else
	if isfield(Sta,'StationInfo') & strcmpi(Sta.StationInfo,'LatLonHeight')

		Sta.staPos = LLA2ecef(Sta.lat*180/pi,Sta.lon*180/pi,Sta.height/1e3)*1e3;

	end

	AdStaPos = getAdjustedPos(Sta.staPos,Sta.Epoch);

	[AdLat,AdLon,AdHeight] = ecef2LLA(AdStaPos/1e3);

	% Convert height to m
	AdLat = AdLat * pi/180;
	AdLon = AdLon * pi/180;
	AdHeight = AdHeight * 1e3;

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function AdStaPos = getAdjustedPos(staPos,EpochVec)

Epoch = EpochVec(1);

for m = 1:size(staPos,2)

	% Convert station position to km

	% Update the Epoch if multiple Epochs are provided
	if length(EpochVec) > 1
		Epoch = EpochVec(m);
	end

	% Get earth, moon, and sun mu
 	mu_e = JATConstant('muEarth')/(1e3)^3;
 	mu_m = JATConstant('muMoon')/(1e3)^3;
 	mu_s = JATConstant('muSun')/(1e3)^3;
    
 	% Get earth, moon, and sun vecotors in IRCF frame
 	earth = ephemDE405('Earth',Epoch,'UTC'); % km
 	moon = ephemDE405('Moon',Epoch,'UTC');   % km
 	sun = ephemDE405('Sun',Epoch,'UTC');     % km
    
	r0 = staPos(:,m)/1e3;
 	% Calculate relative moon-earth and sun-earth vectors
 	RmECI = (moon - earth); % km
 	RsECI = (sun - earth);  % km

 	% Convert to ECEF
 	Rm = jatDCM('eci2ecef',Epoch) * RmECI;
 	Rs = jatDCM('eci2ecef',Epoch) * RsECI;

	% Calcualte unit vectors
 	r_hat = unit(r0);
 	Rm_hat = unit(Rm);
 	Rs_hat = unit(Rs);
    
	% Calucate vector magnitudes
 	r4 = norm(r0)^4;
 	Rm3 = norm(Rm)^3;
 	Rs3 = norm(Rs)^3;
    
	% Love numbers in the solid Earth tide problem (Equations 4.18 in Will Campbell's thesis)
 	h2 = .6090;
 	l2 = .0852;

	% Solid Earth tide equation (Equation 4.22 in Will Campbell's thesis)
 	AdStaPos(:,m) = [mu_m/mu_e*r4/Rm3 * (3*l2*dot(Rm_hat,r_hat)*Rm_hat + ...
	                 r_hat*(3*(h2/2-l2) * dot(Rm_hat,r_hat)^2 - h2/2)) ...
                       + mu_s/mu_e*r4/Rs3 * (3*l2*dot(Rs_hat,r_hat)*Rs_hat + ...
			 r_hat*(3*(h2/2-l2) * dot(Rs_hat,r_hat)^2 - h2/2)) ...
	     	         + r0] * 1e3; % Convert back to m
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = staEarthTide_validation_test();

% STAEARTHTIDE_VALIDATION_TEST - Validation Test for staEarthTide

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/11/2008	 	Original

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
