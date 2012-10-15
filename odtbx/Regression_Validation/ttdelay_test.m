function failed = ttdelay_test(type)
%
% ttdelay_test Regression test for ttdelay
% See also: ttdelay.m
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
%   Ravi Mathur         08/28/2012      Extracted from ttdelay.m

if strcmpi(type,'ValidationTest')

	failed = ttdelay_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = ttdelay_regression_test();
    
else
    disp('ttdelay_test: Unsupported test type ')
    failed = true;
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

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = ttdelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2; % Test predictable output

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.Epoch = datenum('June 10, 2008');
Ephem.StationInfo = 'ECEF';
Ephem.staPos=jatDCM('eci2ecef',Ephem.Epoch) * [JATConstant('rEarth') 0 0]';
Ephem.satPos= [JATConstant('rEarth') 0 0]' * 10;
Ephem.satVel= unit(Ephem.satPos);

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = ttdelay(Ephem,randoms);

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

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelayTemp = ttdelay(Ephem,randoms) %#ok<NOPRT>
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

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = ttdelay(Ephem,randoms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case = 2; % Test predictable output

randoms=1; % time tag std (microsec) (according to Will Campbell's thesis)
Ephem.Epoch = datenum('June 10, 2008');
Ephem.StationInfo = 'ECEF';
Ephem.staPos=jatDCM('eci2ecef',Ephem.Epoch) * [JATConstant('rEarth') 0 0]';
Ephem.satPos= [JATConstant('rEarth') 0 0]' * 10;
Ephem.satVel= unit(Ephem.satPos);

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelay(Case) = ttdelay(Ephem,randoms);

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

%TTDelay = ttdelay(Ephem.Epoch,r_stn_ECEF,scPV,tt_std_micro)
TTDelayTemp = ttdelay(Ephem,randoms) %#ok<NOPRT>
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
