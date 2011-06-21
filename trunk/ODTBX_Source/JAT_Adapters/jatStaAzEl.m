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
% VALIDATION TEST
%
%   To perform a validation test, replace the Ephem structure with 
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
% REGRESSION TEST
%
%   To perform a regression test, replace the Ephem structure with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
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

% Determine whether this is an actual call to the program or a test

if strcmpi(Ephem,'ValidationTest')

	azimuth = jatAzEl_validation_test();
	elevation = 0;

elseif strcmpi(Ephem,'RegressionTest')

	azimuth = jatAzEl_regression_test();
	elevation = 0;

else

	[azimuth,elevation] = getAzEl(Ephem);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function [azimuth,elevation] = getAzEl(Ephem);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jalAzEl_validation_test - jatAzEl Validation Tests

function failed = jatAzEl_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original
%   Keith Speckman         06/27/2008       Changed examples from SEZ to NED
%   Allen Brown            04/09/2009       Updated data for coincident
%                                           satellite and station

disp(' ')
disp('Performing Test....')
disp(' ')

fprintf('%17s%17s%17s%17s\n','Station Pos','Satellite Pos','Expected Az','Az Results')
fprintf('%17s%17s%17s%17s\n','/ Rearth   ','/ Rearth   ',' El (deg)','  El (deg)')
fprintf('%s\n\n',char(ones(1,17*4)*'-'));

tol = 1e-7;

Re = JATConstant('rEarth');

Ephem.Epoch = datenum('May 19, 2008');

Case = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 90;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);


Case = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 270;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = 90;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);
disp('Azimuth undefined')
disp(' ')

Case = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[0 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

% A coincident satellite and station should have az & el of zero.
AzEx = 0;
ElEx = 0;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);
disp('Azimuth undefined')
disp('Elevation undefined')
disp(' ')

Case = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[0 2 1]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = atan(2/1)*180/pi;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 270;
ElEx = -45;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.lat = 0;
Ephem.lon = 0;
Ephem.height = -0.0007;
Ephem.satPos = Ephem.staPos + Re*[1 0 1]';
Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = jatDCM('ecef2eci',Ephem.Epoch) * [Ephem.staPos + Re*[1 0 0]'];
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

AzEx = 0;
ElEx = 90;

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);


Case = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos(:,1) = Re*[1 0 0]';
Ephem.satPos(:,1) = Ephem.staPos(:,1) + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

Ephem.staPos(:,2) = Re*[0 1 0]';
Ephem.satPos(:,2) = Ephem.staPos(:,2) + Re*[1 1 0]';

Ephem.satPos(:,3) = Ephem.staPos(:,2) + Re*[0 1 2]';

AzEx(1,1) = 90;
ElEx(1,1) = 45;
AzEx(1,2) = 270;
ElEx(1,2) = 0;
AzEx(2,1) = 90;
ElEx(2,1) = 0;
AzEx(2,2) = 270;
ElEx(2,2) = 45;
AzEx(3,1) = 45;
ElEx(3,1) = -atan(1/sqrt(8))*180/pi;
AzEx(3,2) = 0;
ElEx(3,2) = atan(1/2)*180/pi;

[azimuthtemp,elevationtemp] = getAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
disp('Matrix output')
disp(' ')

casenum = Case-1;
for n = 1:3
	for m = 1:2
		casenum = casenum + 1;
		azimuth(casenum,1) = azimuthtemp(n,m);
		elevation(casenum,1) = elevationtemp(n,m);
		fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1,m)/Re,Ephem.satPos(1,n)/Re,...
		        AzEx(n,m),azimuth(casenum)*180/pi);
		fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2,m)/Re,Ephem.satPos(2,n)/Re,...
		        ElEx(n,m),elevation(casenum)*180/pi);
		fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3,m)/Re,Ephem.satPos(3,n)/Re);
	end
end

failed = 0;

if exist('jatStaAzEl_ValidationData4_09.m') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	clear AzEx ElEx

	% Load validation values for AzEx and ElEx
	jatStaAzEl_ValidationData4_09

	Diff = [AzEx - azimuth;ElEx - elevation]

	if any( abs(Diff) > tol ) | any(isnan(Diff([1:end])))
		failed = 1;
	end
else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jalAzEl_regression_test - jatAzEl Regression Tests

function failed = jatAzEl_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original
%   Keith Speckman         06/27/2008       Changed examples from SEZ to NED
%   Allen Brown            04/09/2009       Updated data for coincident
%                                           satellite and station

disp(' ')
disp('Performing Demonstration....')

failed = 0;
tol = 1e-7;

Re = JATConstant('rEarth');

Ephem.Epoch = datenum('May 19, 2008');

Case = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[0 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[0 2 1]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.lat = 0;
Ephem.lon = 0;
Ephem.height = -0.0007;
Ephem.satPos = Ephem.staPos + Re*[1 0 1]';
Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = jatDCM('ecef2eci',Ephem.Epoch) * [Ephem.staPos + Re*[1 0 0]'];
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

[azimuth(Case,1),elevation(Case,1)] = getAzEl(Ephem);

Case = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos(:,1) = Re*[1 0 0]';
Ephem.satPos(:,1) = Ephem.staPos(:,1) + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

Ephem.staPos(:,2) = Re*[0 1 0]';
Ephem.satPos(:,2) = Ephem.staPos(:,2) + Re*[1 1 0]';

Ephem.satPos(:,3) = Ephem.staPos(:,2) + Re*[0 1 2]';

[azimuthtemp,elevationtemp] = getAzEl(Ephem);

casenum = Case-1;
for n = 1:3
	for m = 1:2
		casenum = casenum + 1;
		azimuth(casenum,1) = azimuthtemp(n,m);
		elevation(casenum,1) = elevationtemp(n,m);
	end
end
azimuth
elevation

format long g

current_values = [azimuth,elevation];

%save jatStaAzEl_RegressionData6_08 azimuth elevation

failed = 0;

if exist('jatStaAzEl_RegressionData4_09.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	clear azimuth elevation

	% Load regression values for AzEx and ElEx
	load jatStaAzEl_RegressionData4_09

	Diff = [azimuth,elevation] - current_values

	if any(any(abs(Diff)>tol)) | any(isnan([Diff(:,1);Diff([1:end],2)]))
		failed = 1;
	end
else
	failed = 1;
end


end

