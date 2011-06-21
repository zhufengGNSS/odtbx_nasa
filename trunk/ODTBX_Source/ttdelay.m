function TTDelay = ttdelay(Ephem,randoms)
% TTDELAY Compute the station time tag delay error
%
%  TTDelay = ttdelay(Ephem,randoms)
%
%   The station clock time tag error is considered stochastic variables with 
%   mean value of zero.  The standard deviation is based on the current 
%   estimates of deep space navigators.  
%
%   DSN station timekeeping is based on atomic clocks utilizing a hydrogen 
%   maser oscillator [1].  These clocks are maintained to within 1 microsecond 
%   of Coordinated Universal Time.  The range increase due to a time tag error 
%   is described as 
%
%      delta_rho = delta_T * rho_dot
%
%   where delta_T is the time offset.  The time tag error has an expected 
%   value of zero and standard deviation of 1 microsecond.  The delay due to 
%   time tag uncertainty is on the order of 1-10 cm depending on spacecraft 
%   velocity.  
% 
% INPUTS
%      VARIABLE            SIZE         DESCRIPTION (Optional/Default)
%      Ephem.SatPos        (3xN)        ECI satellite coordinates (m)
%      Ephem.SatVel        (3xN)        ECI satellite velocity (m/s)
%      Ephem.Epoch         (1xN)        UTC time in Matlab datenum format
%      randoms              (1x1)       Random timetag error standard deviation
%                                          (us)
%      Ephem.StationInfo   string       Method of specifying ground station
%                                          position
%
%      CASE: Ephem.StationInfo = 'ECEF'
%      -------------------------------------------------------------------------
%      Ephem.staPos        (3xM)        Station ECEF coordinates (m)
%
%      CASE: Ephem.StationInfo = 'LatLonHeight'
%      -------------------------------------------------------------------------
%      Ephem.lat           (1xM)        Station geodetic latitude (rad)
%      Ephem.lon           (1xM)        Station geodetic longitude (rad)
%      Ephem.height        (1xM)        Station geodetic height (m)
%
% OUTPUT:
%      TTDelay              (NxM)       Time Tage Delay (meters)
%
% NOTES:
%   This code is based on Will Campbell's work.  However, his original code
%   looks at the relative velocity magnitude rather than the range rate.  The
%   latter is how the error is described in his thesis.
%   S. Hur-Diaz
%
% VALIDATION TEST
%
%   To perform a validation test, replace the Ephem input with
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
% REGRESSION TEST
%
%   To perform a regression test, replace the Ephem input with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
%
%   keywords: delay, time, time tag
%   See also: jatTropoModel, jatGPSIonoModel, jatIonoDelayModel, lightTimeCorrection
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
%   Keith Speckman      06/10/2008      Original
%   Brent Wm. Barbee    04/21/2009      Removed unnecessary semicolon after
%                                       the function declaration.
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy and
%                                       updated the regression test

% Determine whether this is an actual call to the program or a test

if strcmpi(Ephem,'ValidationTest')

	TTDelay = ttdelay_validation_test();

elseif strcmpi(Ephem,'RegressionTest')

	TTDelay = ttdelay_regression_test();

else

	TTDelay = getTTDelay(Ephem,randoms);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function

function TTDelay = getTTDelay(Ephem,randoms)

% Determine size of outputs (M,N)
N = size(Ephem.satPos,2);

switch Ephem.StationInfo

	% Station Coordinates given in ECEF Coordinates
	case 'ECEF'
		M = size(Ephem.staPos,2);
	case 'LatLonHeight'
		M = size(Ephem.lat,2);
end

% Initialize TTDelay for quicker execution time
TTDelay = zeros(N,M);

e = 1;

for n = 1:N

	if ~isempty(Ephem.Epoch > 1)
		e = n;
	end

	for m = 1:M

		% Convert Lat Lon Height to ECEF Vector if necessary
		if strcmpi(Ephem.StationInfo,'LatLonHeight') && n == 1
			Ephem.staPos(:,m)= LLA2ecef(Ephem.lat(m),Ephem.lon(m),Ephem.height(m))*1e3;
		end

		% Compute the velocity of the station in the ECEF frame (m/sec)
		v_stn_ECEF = cross([0 0 JATConstant('wEarth')]',Ephem.staPos(:,m));

		% Convert the station velocity from ECEF to ECI J2000 (m)
		v_stn_ECI = jatDCM('ecef2eci',Ephem.Epoch(e))*v_stn_ECEF;

		% Determine Station Position in ECI (m)
		r_stn_ECI = jatDCM('ecef2eci',Ephem.Epoch(e))*Ephem.staPos(:,m);

		% Station to satellite relative velocity vecotr in ECI (m/sec)
		RelVelECI = Ephem.satVel(:,n) - v_stn_ECI;

		% Station to satellite Relative position vector in ECI (m)
		RelPosECI = Ephem.satPos(:,n) - r_stn_ECI;

		% Compute range rate by dotting the relative velocity vector with the range unit
		% vector and taking the magnitude (m/sec)
		rho_dot = RelVelECI'*unit(RelPosECI);

		% Convert random error to seconds and calclate the time tag error in 
		% range measurement (m)
		TTDelay(n,m) = randoms*1e-6 * rho_dot;

	end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ttdelay_validation_test - Validation for ttdelay.m

function failed = ttdelay_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/5/2008	 	Original

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 1; % Mars orbit
randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.staPos=[JATConstant('rEarth') 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.Epoch = datenum('June 10, 2008');
[posircf,velircf] = ephemDE405('mars',Ephem.Epoch,'UTC');
Ephem.satPos = (posircf - ephemDE405('earth',Ephem.Epoch,'UTC'))*1e3;
Ephem.satVel = velircf * 1e3;

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = getTTDelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2; % Test predictable output

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.Epoch = datenum('June 10, 2008');
Ephem.StationInfo = 'ECEF';
Ephem.staPos=jatDCM('eci2ecef',Ephem.Epoch) * [JATConstant('rEarth') 0 0]';
Ephem.satPos= [JATConstant('rEarth') 0 0]' * 10;
Ephem.satVel= unit(Ephem.satPos);

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = getTTDelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 3; % Test multiple satellites

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
[lat,lon,height]=ecef2LLA([JATConstant('rEarth') 0 0]');
Ephem.lat = [lat,lat-pi/3];
Ephem.lon = [lon,lon+pi/3];
Ephem.height = [height,height+100];
Ephem.StationInfo = 'LatLonHeight';
Ephem.Epoch = datenum('June 10, 2008');
Ephem.Epoch(2) = datenum('June 10, 2008');
Ephem.Epoch(3) = datenum('June 10, 2008');
[posircf,velircf] = ephemDE405('mars',Ephem.Epoch(1),'UTC');
Ephem.satPos = (posircf - ephemDE405('earth',Ephem.Epoch(1),'UTC'))*1e3;
Ephem.satVel = velircf * 1e3;
[posircf,velircf] = ephemDE405('neptune',Ephem.Epoch(2),'UTC');
Ephem.satPos(:,2) = (posircf - ephemDE405('earth',Ephem.Epoch(2),'UTC'))*1e3;
Ephem.satVel(:,2) = velircf * 1e3;
[posircf,velircf] = ephemDE405('venus',Ephem.Epoch(3),'UTC');
Ephem.satPos(:,3) = (posircf - ephemDE405('earth',Ephem.Epoch(3),'UTC'))*1e3;
Ephem.satVel(:,3) = velircf * 1e3;

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelayTemp = getTTDelay(Ephem,randoms) %#ok<NOPRT>
TTDelay(Case:Case+5) = reshape(TTDelayTemp,1,6) %#ok<NOPRT>

if exist('ttdelay_ValidationData6_09.mat','file') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	% Load the ttdelay outputs
	load ttdelay_ValidationData6_09

	Diff = Case2TTDelay - TTDelay(2);

	if any( abs(Diff) > tol ) ||  any(isnan(TTDelay))
		failed = 1;
	end

else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ttdelay_regression_test - Regression for ttdelay.m

function failed = ttdelay_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/5/2008	 	Original

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 1; % Mars orbit
randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.staPos=[JATConstant('rEarth') 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.Epoch = datenum('June 10, 2008');
[posircf,velircf] = ephemDE405('mars',Ephem.Epoch,'UTC');
Ephem.satPos = (posircf - ephemDE405('earth',Ephem.Epoch,'UTC'))*1e3;
Ephem.satVel = velircf * 1e3;

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = getTTDelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2; % Test predictable output

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.Epoch = datenum('June 10, 2008');
Ephem.StationInfo = 'ECEF';
Ephem.staPos=jatDCM('eci2ecef',Ephem.Epoch) * [JATConstant('rEarth') 0 0]';
Ephem.satPos= [JATConstant('rEarth') 0 0]' * 10;
Ephem.satVel= unit(Ephem.satPos);

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = getTTDelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 3; % Test multiple satellites

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
[lat,lon,height]=ecef2LLA([JATConstant('rEarth') 0 0]');
Ephem.lat = [lat,lat-pi/3];
Ephem.lon = [lon,lon+pi/3];
Ephem.height = [height,height+100];
Ephem.StationInfo = 'LatLonHeight';
Ephem.Epoch = datenum('June 10, 2008');
Ephem.Epoch(2) = datenum('June 10, 2008');
Ephem.Epoch(3) = datenum('June 10, 2008');
[posircf,velircf] = ephemDE405('mars',Ephem.Epoch(1),'UTC');
Ephem.satPos = (posircf - ephemDE405('earth',Ephem.Epoch(1),'UTC'))*1e3;
Ephem.satVel = velircf * 1e3;
[posircf,velircf] = ephemDE405('neptune',Ephem.Epoch(2),'UTC');
Ephem.satPos(:,2) = (posircf - ephemDE405('earth',Ephem.Epoch(2),'UTC'))*1e3;
Ephem.satVel(:,2) = velircf * 1e3;
[posircf,velircf] = ephemDE405('venus',Ephem.Epoch(3),'UTC');
Ephem.satPos(:,3) = (posircf - ephemDE405('earth',Ephem.Epoch(3),'UTC'))*1e3;
Ephem.satVel(:,3) = velircf * 1e3;

%TTDelay = getTTDelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelayTemp = getTTDelay(Ephem,randoms) %#ok<NOPRT>
TTDelay(Case:Case+5) = reshape(TTDelayTemp,1,6) %#ok<NOPRT>

% PreviousTTDelay = TTDelay;
% save ttdelay_RegressionData2_10.mat PreviousTTDelay

if exist('ttdelay_RegressionData2_10.mat','file') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load the ttdelay outputs
	load ttdelay_RegressionData2_10.mat

	Diff = PreviousTTDelay - TTDelay;

	if any( abs(Diff) > tol ) ||  any(isnan(Diff))
		failed = 1;
	end

else
	failed = 1;
end


end
