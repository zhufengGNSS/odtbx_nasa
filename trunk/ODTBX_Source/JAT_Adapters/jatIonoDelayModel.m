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
%  VALIDATION TEST
%
%   The validation test for this function is in a separate file:
%   jatIonoDelayModel_test.m.  This compares the results against an
%   independent formulation and dataset.  There is no validation
%   of the algorithm in this file, only a regression test.
%
%  REGRESSION TEST
%
%   This regression test exercises the various input options for this
%   function and compares against data generated from the JAT
%   IonosphericDelayModel.java's main() routine.
%
%   To perform a regression test, replace the Ephem structure with 
%   'RegressionTest' as the input argument.  If the data file is not in the
%   path this will perform as an example.  The return argument represents
%   regression test failure (1) or success (0).
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

% Determine whether this is an actual call to the program or a test

if strcmpi(Ephem,'RegressionTest')

	IonoDelay = jatIonoDelayModel_regression_test();
else

	IonoDelay = getIonoDelay(Ephem);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function IonoDelay = getIonoDelay(Ephem)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = jatIonoDelayModel_regression_test()
% jatIonoDelayModel_regression_test - jatIonoDelayModel Regression Tests
% 
% This function checks the jatIonoDelayModel function against the test case
% contained in IonosphericDelayModel.java.  The java file contains just a
% station lat, lon, height example.  This script converts the lat, lon, height
% to an ECEF position in order to validate both input paths.

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/16/2008	 	Original
%   Allen Brown        05/01/2009       Updated data file from JAT and
%                                       fixed unit conversion bugs.

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-10; % changed to 1e-10 for 5_09 data, where differences from JAT were ~2e-12

% Test case taken from IonosphericDelayModel.java - main()
deg = [-90:0.01:90-0.01];
lat = 0;  % java: rad
lon = 0; % java: rad
height = 10; % java: m
Ephem.Lat = lat*180/pi; % jatIonoDelayModel: deg
Ephem.Lon = lon*180/pi; % jatIonoDelayModel: deg
Ephem.Height = height/1000; % jatIonoDelayModel: km
Ephem.Epoch = mJD2MatlabTime(54103.0 + 1.7361e-006)*ones(size(deg)); 
Ephem.StationInfo = 'LatLonHeight';

x = jat.constants.WGS84.R_Earth*cos(deg*pi/180);
y = jat.constants.WGS84.R_Earth*sin(deg*pi/180);
z = zeros(size(x));
Ephem.SatPos = [x;y;z];
Ephem.SatCoords = 'ECEF';

% Initialize Ephem.StaPos to speed things up
Ephem.StaPos = Ephem.SatPos * 0;

% Convert test case to ECEF to test that call method (rad & m)
station = jat.groundstations.GroundStation('tmp',lat,lon,height);
Ephem.StaPos = station.getECEFPosition();

% Test ECEF input method
Ephem.StationInfo = 'ECEF';
ionoDelay_ECEF = jatIonoDelayModel(Ephem);

% Test Lat, Lon, Height input method
Ephem.StationInfo = 'LatLonHeight';
ionoDelay_LLH = jatIonoDelayModel(Ephem);

disp(' ')
disp(['Station Lat = ',num2str(Ephem.Lat(1)),'deg'])
disp(['Station Lon = ',num2str(Ephem.Lon(1)),'deg'])
disp(['Station Height = ',num2str(Ephem.Height(1)),'km'])

disp(' ')
disp('******************************************** Satellite ECEF ************************************************')
disp('                X                         Y                         Z                      Delay')
disp([x',y',z',ionoDelay_ECEF])

disp(' ')
disp('********************************************* Satellite  LLH ************************************************')
disp('               Lat                       Lon                      Height                   Delay')
disp([x',y',z',ionoDelay_LLH])


if exist('IonoDelayModelJatOutputs5_09.txt','file') == 2
	disp(' ')
	disp('Performing Regression...')
	disp(' ')


	% Load the JAT IonosphericDelayModel outputs.  This file is iono.txt
	% from JAT's IonosphericDelayModel.main() output.
	load IonoDelayModelJatOutputs5_09.txt

	format long g

	% Compare the above results against the JAT outputs
	maxerrors_ECEF = max(abs([x',y',z',ionoDelay_ECEF]-IonoDelayModelJatOutputs5_09));
	maxerrors_LLH = max(abs([x',y',z',ionoDelay_LLH]-IonoDelayModelJatOutputs5_09));
    % we're interested in the delays:
    Diff = [ maxerrors_ECEF(:,4); maxerrors_LLH(:,4)];

    if any( abs(Diff) > tol ) ||  any(isnan(Diff))
        failed = 1;
        disp('FAILED: Comparison with JAT IonosphericDelayModel.');
    else
        disp('Comparison with JAT IonosphericDelayModel passed.');
    end
    
	% Next test: Test different sized inputs and satellite ECI input path

    % tweak the tolerance to 1mm since we're adjusting station location
    tol = 0.001; %m
    
	Ephem.Epoch = Ephem.Epoch(1:4);
	for k = 1:4
		Ephem.SatPos(:,k) = jatDCM('ecef2eci',Ephem.Epoch(k))* Ephem.SatPos(:,k);
	end
	Ephem.SatPos = Ephem.SatPos(:,1:4);
	Ephem.SatCoords = 'ECI';

	staPos_km = LLA2ecef(0,0,10/1000); % deg, deg, km
    Ephem.StaPos = staPos_km*1000; % Ephem.StaPos must be in m
 	Ephem.StaPos(:,2) = Ephem.StaPos(:,1) + 1; % delta by 1 m
 	Ephem.StationInfo = 'ECEF';
 
 	Ephem.SignalFreq = 1.6001e9; % tweak just to test the input
 
 	ionoDelayMat = jatIonoDelayModel(Ephem);
    
    if size(ionoDelayMat) ~= [4 2]
        failed = 1;
        disp('Unexpected results size with multiple stations and satellite positions.');
    end
    
    disp(' ')
    disp('******************************************** Satellite ECI ************************************************')
    disp('Station 1')
    disp('                X                         Y                         Z                      Delay')
    disp([Ephem.SatPos(1,:)',Ephem.SatPos(2,:)',Ephem.SatPos(3,:)',ionoDelayMat(:,1)]);
    disp('Station 2')
    disp('                X                         Y                         Z                      Delay')
    disp([Ephem.SatPos(1,:)',Ephem.SatPos(2,:)',Ephem.SatPos(3,:)',ionoDelayMat(:,2)]);
    
    dlys = IonoDelayModelJatOutputs5_09(1:4,4);
    Diff = [ionoDelayMat(:,1); ionoDelayMat(:,2)] - [dlys; dlys];
    if any( abs(Diff) > tol ) ||  any(isnan(Diff))
        failed = 1;
        disp('FAILED: Comparison of multiple stations and multiple satellites.');
    else
        disp('Comparison of multiple stations and multiple satellites passed.');
    end
 
% 	% Test multiple sat frequencies
% 
% 	Ephem.SatPos = Ephem.SatPos(:,[1 1 1]);
% 	Ephem.StaPos = Ephem.StaPos(:,1);
% 	Ephem.SignalFreq = [1.6,2,2.4]*1e9;
% 
% 
% 	ionoDelayVec = jatIonoDelayModel(Ephem)


else
	failed = 1;
end

end
