function [azimuth,elevation] = jatStaAzEl(Ephem)

% JATSTAAZEL Returns Azimuth and Elevation of satellites relative to ground stations
%
%   [azimuth,elevation] = jatStaAzEl(Ephem)
%
%   This function outputs the azimuth and elevation of a spacecraft relative to
%   a ground station.  Solutions are calculated in the Tropocentric (NED)
%   reference frame.  Azimuth is defined as the angle on the surface tangent
%   plane from the North reference point in the clockwise direction (due east is
%   90 degrees azimuth).  Elevation is the angle above the surface tangent
%   plane.
%
%   The ground station location may be specified in Latitude, Longitude, and 
%   Height coordinates, or in ECEF coordinates.  The spacecraft location may be
%   specified in ECI J2000 coordinates or in ECEF coordinates.
%
%   Since Azimuth is undefined at an elevation of pi/2 radians (90 degrees), The
%   azimuth is set to zero if the elevation is within 1e-15 radians.
%   Co-locating the spacecraft and ground station results in zero for
%   azimuth and elevation.
%   
% INPUTS 
%    VARIABLE           SIZE      DESCRIPTION (Optional/Default)
%    Ephem.satPos       (3xN)     ECEF or ECI satellite coordinates,
%                                 in meters
%    Ephem.SatCoords    string    (Optional) Satellite coordinate system
%                                   'ECEF' (Default) or 'ECI'
%    Ephem.Epoch        (1xN)     UTC epoch in Matlab datenum time format
%                                   (Optional - Required only for ECI SatCoords)
%    Ephem.StationInfo  string    Method of specifying ground station position
%
%   CASE: Ephem.StationInfo = 'ECEF'
%   ----------------------------------------------------------------------------
%    Ephem.staPos       (3xM)     Station ECEF coordinates (meters)
%
%   CASE: Ephem.StationInfo = 'LatLonHeight'
%   ----------------------------------------------------------------------------
%    Ephem.lat          (1xM)     Station geodetic latitude (deg) (note 
%                                 outputs are in radians)
%    Ephem.lon          (1xM)     Station geodetic longitude (deg) (note 
%                                 outputs are in radians)
%    Ephem.height       (1xM)     Station geodetic height (meters)
%
% OUTPUTS 
%    azimuth            (NxM)     Satelite azimuth relative to ground station 
%                                   (0 to 2*pi rad, note that Ephem.lat &
%                                   lon use degrees)
%    elevation          (NxM)     Satelite elevation relative to ground station
%                                   (-pi/2 to +pi/2 rad, note that 
%                                   Ephem.lat & lon use degrees)
%
% VALIDATION/REGRESSION TEST
%
%   These tests were moved to jatStaAzEl_test.m to conform to the new
%   regression testing framework.
%
%  keywords: JAT Adapter, azimuth, elevation
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
%   Author      		Date        	  Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman      05/15/2008      Original
%   Keith Speckman      06/27/2008      Incorporated Jat change from SEZ to
%                                       NED
%   Allen Brown         02/27/2009      Fixed documentation &
%                                       units problems
%   Allen Brown         04/09/2009      Updated data for coincident
%                                       satellite and station
%   Kevin Berry         06/22/2009      Corrected the help to specify the
%                                       input time scale of UTC  
%   Ravi Mathur         08/29/2012      Extracted regression test

% Determine size of outputs (M,N)
N = size(Ephem.satPos,2);

% Set epsilon for elevations near 90 degrees
epsilon = 1e-15;

switch Ephem.StationInfo

	% Station Coordinates given in ECEF Coordinates
	case 'ECEF'
		M = size(Ephem.staPos,2);
	case 'LatLonHeight'
		M = size(Ephem.lat,2);
end

% Initialize azimuth and elevation for quicker execution time
azimuth = zeros(N,M);
elevation = zeros(N,M);

for n = 1:N

	% Determine ECEF satPos

	satPosECEF = Ephem.satPos(:,n);

	if isfield(Ephem,'SatCoords')
		if strcmpi(Ephem.SatCoords,'ECI')

			% Convert satPos to ECEF (m)
			satPosECEF = jatDCM('eci2ecef',Ephem.Epoch(n))*Ephem.satPos(:,n);

		end
	end

	for m = 1:M

		% Convert Lat Lon Height to ECEF Vector if necessary
		if strcmpi(Ephem.StationInfo,'LatLonHeight') & n == 1
            % LLA2ecef is deg & km in, km out, so convert to meters
			Ephem.staPos(:,m) = LLA2ecef(Ephem.lat(m),Ephem.lon(m),Ephem.height(m))*1000;
		end

		% Create Java station object (input in m)
		station = jat.groundstations.GroundStation('tmp',Ephem.staPos(:,m));

		% Calculate Az and El
		azel = station.azEl(satPosECEF); % input in m
		azimuth(n,m) = azel(1);
		elevation(n,m) = azel(2);

		% Set azimuth to zero if elevation is near 90 degrees
		if abs(abs(elevation(n,m)) - pi/2) < epsilon
			azimuth(n,m) = 0;
		end
	end
end

end