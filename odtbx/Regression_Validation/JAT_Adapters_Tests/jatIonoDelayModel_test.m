function failed = jatIonoDelayModel_test(type, ploton)
% Regression and Validation test for jatIonoDelayModel. The validation test
% compares with an independent formulation of jatIonoDelayModel.
%
% Inputs:
%   type               'RegressionTest' or 'ValidationTest'
%   ploton  (optional) true/false whether the comparison plot should be
%                      displayed (not usually provided when running as
%                      part of a test suite)
%
% Outputs:
%   failed    1=test failure, 0=test success
%
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
%
%   Ravi Mathur         08/29/2012      Extracted regression test from
%                                       jatIonoDelayModel.m, and factored
%                                       Reg/Val tests into subfunctions.

if strcmpi(type,'ValidationTest')
    if nargin == 2 % Run validation with specified plotting state
        failed = jatIonoDelayModel_validation_test(ploton);
        
    else % Run validation without plotting 
        failed = jatIonoDelayModel_validation_test(false);
    end

elseif strcmpi(type,'RegressionTest')

	failed = jatIonoDelayModel_regression_test();

else
    disp('jatIonoDelayModel_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jatIonoDelayModel_validation_test - jatIonoDelayModel Validation Test
%
% Equations and data for validation test generated from Campbell's thesis
% https://fftbwall.gsfc.nasa.gov/odtbx/trac.cgi/attachment/wiki/IncrementSix/WillCampbell.zip
% equations and data: S. Hur-Diaz
% regression test adaptation: A. Brown
%
%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Ravi Mathur         08/29/2012      Created subroutine

function failed = jatIonoDelayModel_validation_test(ploton)

failed = 0; % default to not failed

%
% Comparison data from Campbell thesis:
%
% Ground station data
lam = 0; % longitue of ground station in deg
phi = 0; % latitude of ground station in deg
f = 7650e6; % frequency (Hz)

% Satellite data
el = 90; % elevation of satellite in deg
az = 0;  % azimuth of satellite in deg

d2r=pi/180;

rE=6378; % Earth radius in km
h = 350; % mean ionospheric height in km
T=32; % hours
thd=14; % hours
N = 3.192e18; % Hz^2-m
D = 3.192e19; % Hz^2-m

%
% equations coded from Campbell thesis:
%
zeta=asin(rE*cos(el*d2r)/(rE+h)); % rad

phip=asin(sin(phi*d2r)*sin(el*d2r+zeta)+cos(phi*d2r)*cos(az*d2r)*cos(el*d2r+zeta)); % rad
lamp = lam+asin((sin(az*d2r)*cos(el*d2r+zeta))/cos(phip))/d2r; % deg

st=0:1:24;
drho=zeros(length(st),1);
for i=1:length(st)
    UT=st(i);
    tp = mod(lamp/15 + UT,24); % hours
    chi = min(pi/2, 2*pi*abs(tp-thd)/T);

    drho(i)=(N+D*cos(chi))/f^2/cos(zeta);
end

if ploton
    figure
    plot(st,drho,'xb-'),grid
    title('Iono Delay vs. Station local time')
    xlabel('hrs'),ylabel('meters')
end

%
% Setup case to run the above using jatIonoDelayModel.m
%
ephem.Lat = 0;
ephem.Lon = 0;
ephem.Height = 0;
ephem.StationInfo = 'LatLonHeight';
ephem.SatCoords = 'ECEF';
ephem.SatPos = [1 0 0]'*7000e3;
ephem.SignalFreq = 7650e6;

drhoiono=zeros(length(st),1);
for i=1:length(st)
    
    ephem.Epoch = datenum(2000,9,19,st(i),0,0);
    drhoiono(i) = jatIonoDelayModel(ephem);
    
end

if ploton
    hold on
    plot(st,drhoiono,'or--'),grid
    title('Iono Delay from jatIonoDelayModel vs. UTC')
    xlabel('hrs'),ylabel('meters')
    hold off
end

delta=drho-drhoiono;
maxerr=max(abs(delta));

% Test criteria: less than 1 mm difference in calculated delay:
if maxerr >= 0.001
    failed = 1;
    disp(sprintf('jatIonoDelayModel_test failed with a maximum error of %f, which is greater than 1mm.',maxerr));
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
%   Ravi Mathur        08/29/2012      Extracted from jatIonoDelayModel.m

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

