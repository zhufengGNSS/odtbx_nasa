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
%  VALIDATION/REGRESSION TEST
%
%   These tests have been moved to staEarthTide_test.m to conform to the
%   new regression testing format.
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
%   Ravi Mathur         08/28/2012      Extracted regression test

% Initialize LLH outputs
AdLat = 0;
AdLon = 0;
AdHeight = 0;


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