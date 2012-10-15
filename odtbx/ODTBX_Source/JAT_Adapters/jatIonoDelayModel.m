function [IonoDelay] = jatIonoDelayModel(Ephem)

% JATIONODELAYMODEL  Ionospheric delay distance between satellites and ground stations
%
%  [IonoDelay] = jatIonoDelayModel(Ephem)
%   
%   jatIonoDelayModel returns the Ionospheric delay distance between satellites
%   and ground stations.  The Ephem structure may contain arrays of N satellites
%   and M ground stations.  Only the Ephem fields necessary for the chosen
%   StationInfo case are necessary.  A default frequency of 1.6 Ghz is used
%   unless otherwise specified.  The array lengths of satellite position
%   and epoch must be the same.  If signal frequency is specified then it
%   should also be the same length as sat position and epoch or it should
%   be a scalar (assumed to be the same value for all combinations).
%
%   The ground station location may be specified in Latitude, Longitude, and 
%   Height coordinates, or in ECEF coordinates.  The spacecraft location may be
%   specified in ECI J2000 coordinates or in ECEF coordinates.
%
%  INPUTS 
%      VARIABLE            SIZE         DESCRIPTION (Optional/Default)
%      Ephem.SatPos        (3xN)        ECEF or ECI satellite coordinates (m)
%      Ephem.SatCoords     string       (Optional) Satellite coordinate system
%                                          'ECI' (Default) or 'ECEF'
%      Ephem.Epoch         (1xN)        UTC epoch in Matlab datenum time
%                                          format
%      Ephem.SignalFreq (1xN or 1x1)    (Optional) Signal Frequency (Hz)
%                                          Default = 1.6e9
%      Ephem.StationInfo   string       Method of specifying ground station
%                                          position
%
%      CASE: Ephem.StationInfo = 'ECEF'
%      -------------------------------------------------------------------------
%      Ephem.StaPos        (3xM)        Station ECEF coordinates (m)
%
%      CASE: Ephem.StationInfo = 'LatLonHeight'
%      -------------------------------------------------------------------------
%      Ephem.Lat           (1xM)        Station geodetic latitude (deg)
%      Ephem.Lon           (1xM)        Station geodetic longitude (deg)
%      Ephem.Height        (1xM)        Station geodetic height (km)
%
%  OUTPUTS 
%      IonoDelay           (NxM)        Ionospheric delay distance between
%					   satellites and ground stations (m)
%
%  REFERENCE
%
%   Burkhart, Paul D., "Adaptive Orbit Determination for Interplanetary
%   Spacecraft," Ph.D. dissertation, The University of Texas at Austin,
%   May 1995, pp. 22-27.
%
% VALIDATION/REGRESSION TEST
%
%   These tests were moved to jatIonoDelayModel_test.m to conform to the
%   new regression testing framework.
%
%   keywords: JAT Adapter, delay, time, Ionosphere
%   See also: jatTropoModel, jatGPSIonoModel
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
%   Keith Speckman      05/09/2008	 	Original
%   Allen Brown         05/01/2009      updated with new JAT comparison
%                                       data, removed old regression test
%                                       in favor of existing validation
%                                       test (and external validation test)
%   Kevin Berry         06/22/2009      Corrected the help to specify the
%                                       input time scale of UTC  
%   Ravi Mathur         08/29/2012      Extracted regression test

% Determine size of outputs (M,N)
N = size(Ephem.SatPos,2);

switch Ephem.StationInfo

	% Station Coordinates given in ECEF Coordinates
	case 'ECEF'
		M = size(Ephem.StaPos,2);
	case 'LatLonHeight'
		M = size(Ephem.Lat,2);
end

% Determine size of SignalFrequency
if isfield(Ephem,'SignalFreq');
	FreqSize = size(Ephem.SignalFreq,2);
else
	FreqSize = 0;
end

% Initialize IonoDelay for quicker execution time
IonoDelay = zeros(N,M);

for n = 1:N

	% Create a java object containing multiple time formats sync'd to the input time (converted
	% to MJD)
	T=jat.spacetime.Time(matlabTime2MJD(Ephem.Epoch(n)));

	if isfield(Ephem,'SatCoords') && strcmpi(Ephem.SatCoords,'ECEF')

		% Determine ECEF satPos
		satPosECEF = Ephem.SatPos(:,n);

	else

		% Convert satPos to ECEF
		satPosECEF = jatDCM('eci2ecef',Ephem.Epoch(n))*Ephem.SatPos(:,n);

	end

	for m = 1:M

		% Convert Lat Lon Height to ECEF Vector if necessary
		if strcmpi(Ephem.StationInfo,'LatLonHeight') && n == 1
			Ephem.StaPos(:,m)= LLA2ecef(Ephem.Lat(m),Ephem.Lon(m),Ephem.Height(m))*1e3;
		end

		% Update Signal Frequency if desired and create ionosphericDelay java object
		if FreqSize == 1
			ionosphericDelayModel = ...
			   jat.ground_tracking.IonosphericDelayModel(Ephem.StaPos(:,m),...
			                                             Ephem.SignalFreq);
		elseif FreqSize > 1
			ionosphericDelayModel = ...
			   jat.ground_tracking.IonosphericDelayModel(Ephem.StaPos(:,m),...
			                                             Ephem.SignalFreq(n));
		else
			ionosphericDelayModel = ...
			        jat.ground_tracking.IonosphericDelayModel(Ephem.StaPos(:,m));
		end

		% Calculate Delay
		slant_ionoDelay = ionosphericDelayModel.computeDelay(T,satPosECEF);
        % slant angle is computed by JAT but never used here.
		%SlantAngle(n,m) = slant_ionoDelay(1);
		IonoDelay(n,m) = slant_ionoDelay(2);
	end
end

end