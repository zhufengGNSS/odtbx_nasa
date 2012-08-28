function failed = relativityLightDelay_test(type)
%
% relativityLightDelay_test Regression test for relativityLightDelay
% See also: relativityLightDelay.m
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
%   Ravi Mathur         08/28/2012      Extracted from relativityLightDelay.m

if strcmpi(type,'ValidationTest')

	failed = relativityLightDelay_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = relativityLightDelay_regression_test();
else

    disp('relativityLightDelay_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation Test

function failed = relativityLightDelay_validation_test()

rEarth = JATConstant('rEarth');
rotmat = [cos(5*pi/180) -sin(5*pi/180) 0
          sin(5*pi/180) cos(5*pi/180)  0
	  0             0              1 ];
rotmat2 = [cos(2*pi/180) -sin(2*pi/180) 0
          sin(2*pi/180) cos(2*pi/180)  0
	  0             0              1 ];
rotmat3 = [cos(1*pi/180) -sin(1*pi/180) 0
          sin(1*pi/180) cos(1*pi/180)  0
	  0             0              1 ];

%%%%%%%%%%%%%% Case 1 - Check multiple stations and spacecraft
Ephem.Epoch = [datenum('Aug 28, 2008'),datenum('Jan 4, 2005')]; %UTC

Ephem.satPos = JATConstant('rEarth')*[20*[1 1 0.3]',500*[1.2 0.8 -.2]'];

% ECEF Stn Pos (m)
Ephem.staPos(:,1) = [-2350443.812; -4651980.837; 3665630.988 ];
Ephem.staPos(:,2) = JATConstant('rEarth')*[1 0 0]';
Ephem.staPos(:,3) = JATConstant('rEarth')*[0 1 0]';
Ephem.SatCoords = 'ECI';

Ephem.StationInfo = 'ECEF';

RelLightDelay{1} = relativityLightDelay(Ephem);
disp(RelLightDelay{1})


%%%%%%%%%%%%%% Case 2 - Check the same case with different input coords
Ephem.Epoch = [datenum('Aug 28, 2008'),datenum('Jan 4, 2005')]; %UTC

Ephem.satPos = JATConstant('rEarth')*[20*[1 1 0.3]',500*[1.2 0.8 -.2]'];
Ephem.satPos = [jatDCM('eci2ecef',Ephem.Epoch(1))*Ephem.satPos(:,1), ...
                jatDCM('eci2ecef',Ephem.Epoch(2))*Ephem.satPos(:,2)];

% ECEF Stn Pos (m)
[Ephem.lat(1,1),Ephem.lon(1,1),Ephem.height(1,1)] = ...
                                  ecef2LLA([-2350443.812; -4651980.837; 3665630.988 ]/1e3);
[Ephem.lat(1,2),Ephem.lon(1,2),Ephem.height(1,2)] = ...
                                  ecef2LLA(JATConstant('rEarth')*[1 0 0]'/1e3);
[Ephem.lat(1,3),Ephem.lon(1,3),Ephem.height(1,3)] = ...
                                  ecef2LLA(JATConstant('rEarth')*[0 1 0]'/1e3);

Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

RelLightDelay{2} = relativityLightDelay(Ephem);
disp(RelLightDelay{2})

%%%%%%%%%%%%%% Case 3 - Check multiple stations and spacecraft

Ephem.Epoch = [datenum('Aug 28, 2008'),...
               datenum('Jan 4, 2005'),...
	       datenum('June 25, 2002'),...
	       datenum('March 3, 2007'),...
	       datenum('December 20, 2004'),...
	       datenum('December 20, 2004'),...
	       datenum('December 20, 2004')];

% ECEF Stn Pos (m)
Ephem.staPos(:,1) = [-2350443.812; -4651980.837; 3665630.988 ];
Ephem.staPos(:,2) = rEarth*[1 0 0]';
Ephem.staPos(:,3) = rEarth*[0 1 0]';

% LEO
Ephem.satPos(:,1) = (rEarth+400e3)*[1 1 0.3]';

% GEO
Ephem.satPos(:,2) = (rEarth+35786e3)*[1 1 0.3]';

% Lunar
Ephem.satPos(:,3) = rotmat*[ephemDE405('Moon',Ephem.Epoch(3),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(3),'UTC')] *1e3* 1.05;

% Mars
Ephem.satPos(:,4) = rotmat*[ephemDE405('Mars',Ephem.Epoch(4),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(4),'UTC')] *1e3* 1.05;

% Jupiter
Ephem.satPos(:,5) = rotmat*[ephemDE405('Jupiter',Ephem.Epoch(5),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(5),'UTC')] *1e3* 1.05;
Ephem.satPos(:,6) = rotmat2*[ephemDE405('Jupiter',Ephem.Epoch(6),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(6),'UTC')] *1e3* 1.05;
Ephem.satPos(:,7) = rotmat3*[ephemDE405('Jupiter',Ephem.Epoch(7),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(7),'UTC')] *1e3* 1.05;

Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

RelLightDelay{3} = relativityLightDelay(Ephem);

disp(' ')
disp('LEO, GEO, Lunar, Mars, and Jupiter for three groundstations (m)')
disp(' ')
disp(RelLightDelay{3})

%%%%%%%%%%%%% Case 4 - Check against original dsn equations
Ephem4.Epoch = datenum('Aug 28, 2008');

Ephem4.satPos = JATConstant('rEarth')*20*[1 1 0.3]';

% ECEF Stn Pos (m)
Ephem4.staPos = [-2350443.812
              -4651980.837
               3665630.988 ];

Ephem4.StationInfo = 'ECEF';
Ephem4.SatCoords = 'ECI';

RelLightDelay{4} = relativityLightDelay(Ephem4);

%%%%%%%%%%%%%% Case 3 - Run original equations
%addpath /matlab/dsn
%t3 = matlabTime2JD(Ephem4.Epoch);
%sc_t3 = (Ephem4.satPos'+ ephemDE405('EARTH',Ephem4.Epoch)' * 1e3) / JATConstant('au');
%receive_stn = Ephem4.staPos'/JATConstant('au');
%RelLightDelaydsn{1} = dsn_relativity(t3,sc_t3,receive_stn)*JATConstant('c')/...
%                      JATConstant('sec2Days');
%
%for k = 1:3
%	for m = 1:7
%		t3 = matlabTime2JD(Ephem.Epoch(m));
%		sc_t3 = (Ephem.satPos(m)'+ ephemDE405('EARTH',Ephem.Epoch(m),'UTC')' * 1e3) / ...
%		         JATConstant('au');
%		receive_stn = Ephem.staPos(k)'/JATConstant('au');
%		RelLightDelaydsn{2}(m,k) = ...
%		     dsn_relativity(t3,sc_t3,receive_stn)*JATConstant('c')/JATConstant('sec2Days');
%	end
%end
%
%save relativityLightDelay_ValidationData11_08 RelLightDelaydsn
%
%tol = 1.1e-1;
%failed = 0;
%if exist('relativityLightDelay_ValidationData11_08.mat') == 2
%	disp(' ')
%	disp('Performing Validation...')
%	disp(' ')
%
%	% Load validation values
%	load relativityLightDelay_ValidationData11_08
%
%	Diff1 = RelLightDelaydsn{1} - RelLightDelay{4}
%	Diff2 = RelLightDelaydsn{2} - RelLightDelay{3}
%
%	if any( abs(Diff1) > tol ) | any(isnan(Diff1)) ...
%	 | any( abs(Diff2) > tol ) | any(isnan(Diff2))
%		failed = 1;
%	end
%else
%	failed = 1;
%end

% Matching Plot from Landis Markley

dayRange = [-300:10:300];
% Add 2 Venus synodic periods to the central epoch to avoid MJD out of
% bounds error from the being outside the leap-second table when calling
% the jatDCM function below.
Ephem.Epoch = dayRange*(datenum('19-May-2000')-datenum('18-May-2000'))+ ...
              datenum('25-January-1970')+583.9*2;  %UTC
numcases = size(Ephem.Epoch,2);
for k = 1:numcases
	Ephem.satPos(:,k) = (ephemDE405('Venus',Ephem.Epoch(k),'UTC') - ephemDE405('Earth',Ephem.Epoch(k),'UTC'))*1e3;
end
Ephem.SatCoords = 'ECI';
Ephem.StationInfo = 'ECEF';
midDateIndex = floor(numcases/2)+1;
Ephem.staPos = (jatDCM('eci2ecef',Ephem.Epoch(midDateIndex))*...
               Ephem.satPos(:,midDateIndex))/norm(Ephem.satPos(:,midDateIndex))*JATConstant('rEarth');

% Calculate 2-way delay time in microsecs
RelLightDelayPlot = relativityLightDelay(Ephem)/JATConstant('c')*1e6*2;
FiniteIndex = find(isfinite(RelLightDelayPlot));

figure(1)
	subplot(121)
		plot(dayRange(FiniteIndex),RelLightDelayPlot(FiniteIndex),'*')
		axis([dayRange(1),dayRange(end),0 200])
		title('ODTBX Relativity Light Delay Results')
		ylabel('Delay (\musec)')
		xlabel('Days from Jan 25, 1970 plus 2 synodic periods')
	handle = subplot(122)
		RGB = imread('Rel Plot.jpg');
		image(RGB)
		set(handle,'XTickMode','manual')
		set(handle,'XTickLabelMode','manual')
		set(handle,'XTickLabel',[])
		set(handle,'XTick',[])
		set(handle,'YTickMode','manual')
		set(handle,'YTickLabelMode','manual')
		set(handle,'YTickLabel',[])
		set(handle,'YTick',[])
		
failed = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Regression Test

function failed = relativityLightDelay_regression_test()

rEarth = JATConstant('rEarth');
rotmat = [cos(5*pi/180) -sin(5*pi/180) 0
          sin(5*pi/180) cos(5*pi/180)  0
	  0             0              1 ];
rotmat2 = [cos(2*pi/180) -sin(2*pi/180) 0
          sin(2*pi/180) cos(2*pi/180)  0
	  0             0              1 ];
rotmat3 = [cos(1*pi/180) -sin(1*pi/180) 0
          sin(1*pi/180) cos(1*pi/180)  0
	  0             0              1 ];

%%%%%%%%%%%%%% Case 1 - Check multiple stations and spacecraft
Ephem.Epoch = [datenum('Aug 28, 2008'),datenum('Jan 4, 2005')]; %UTC

Ephem.satPos = JATConstant('rEarth')*[20*[1 1 0.3]',500*[1.2 0.8 -.2]'];

% ECEF Stn Pos (m)
Ephem.staPos(:,1) = [-2350443.812; -4651980.837; 3665630.988 ];
Ephem.staPos(:,2) = JATConstant('rEarth')*[1 0 0]';
Ephem.staPos(:,3) = JATConstant('rEarth')*[0 1 0]';
Ephem.SatCoords = 'ECI';

Ephem.StationInfo = 'ECEF';

RelLightDelay{1} = relativityLightDelay(Ephem);
disp(RelLightDelay{1})


%%%%%%%%%%%%%% Case 2 - Check the same case with different input coords
Ephem.Epoch = [datenum('Aug 28, 2008'),datenum('Jan 4, 2005')]; %UTC

Ephem.satPos = JATConstant('rEarth')*[20*[1 1 0.3]',500*[1.2 0.8 -.2]'];
Ephem.satPos = [jatDCM('eci2ecef',Ephem.Epoch(1))*Ephem.satPos(:,1), ...
                jatDCM('eci2ecef',Ephem.Epoch(2))*Ephem.satPos(:,2)];

% ECEF Stn Pos (m)
[Ephem.lat(1,1),Ephem.lon(1,1),Ephem.height(1,1)] = ...
                                  ecef2LLA([-2350443.812; -4651980.837; 3665630.988 ]/1e3);
[Ephem.lat(1,2),Ephem.lon(1,2),Ephem.height(1,2)] = ...
                                  ecef2LLA(JATConstant('rEarth')*[1 0 0]'/1e3);
[Ephem.lat(1,3),Ephem.lon(1,3),Ephem.height(1,3)] = ...
                                  ecef2LLA(JATConstant('rEarth')*[0 1 0]'/1e3);

Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

RelLightDelay{2} = relativityLightDelay(Ephem);
disp(RelLightDelay{2})

%%%%%%%%%%%%%% Case 3 - Check multiple stations and spacecraft

Ephem.Epoch = [datenum('Aug 28, 2008'),...
               datenum('Jan 4, 2005'),...
	       datenum('June 25, 2002'),...
	       datenum('March 3, 2007'),...
	       datenum('December 20, 2004'),...
	       datenum('December 20, 2004'),...
	       datenum('December 20, 2004')];

% ECEF Stn Pos (m)
Ephem.staPos(:,1) = [-2350443.812; -4651980.837; 3665630.988 ];
Ephem.staPos(:,2) = rEarth*[1 0 0]';
Ephem.staPos(:,3) = rEarth*[0 1 0]';

% LEO
Ephem.satPos(:,1) = (rEarth+400e3)*[1 1 0.3]';

% GEO
Ephem.satPos(:,2) = (rEarth+35786e3)*[1 1 0.3]';

% Lunar
Ephem.satPos(:,3) = rotmat*[ephemDE405('Moon',Ephem.Epoch(3),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(3),'UTC')] *1e3* 1.05;

% Mars
Ephem.satPos(:,4) = rotmat*[ephemDE405('Mars',Ephem.Epoch(4),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(4),'UTC')] *1e3* 1.05;

% Jupiter
Ephem.satPos(:,5) = rotmat*[ephemDE405('Jupiter',Ephem.Epoch(5),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(5),'UTC')] *1e3* 1.05;
Ephem.satPos(:,6) = rotmat2*[ephemDE405('Jupiter',Ephem.Epoch(6),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(6),'UTC')] *1e3* 1.05;
Ephem.satPos(:,7) = rotmat3*[ephemDE405('Jupiter',Ephem.Epoch(7),'UTC') - ...
                            ephemDE405('Earth',Ephem.Epoch(7),'UTC')] *1e3* 1.05;

Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

RelLightDelay{3} = relativityLightDelay(Ephem);

disp(' ')
disp('LEO, GEO, Lunar, Mars, and Jupiter for three groundstations (m)')
disp(' ')
disp(RelLightDelay{3})

% RelLightDelayRegress = RelLightDelay;
% save relativityLightDelay_RegressionData06_09 RelLightDelayRegress

tol = 1e-11;
failed = 0;
if exist('relativityLightDelay_RegressionData06_09.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load regression values
	load relativityLightDelay_RegressionData06_09

	Diff1 = RelLightDelayRegress{1} - RelLightDelay{1}
	Diff2 = RelLightDelayRegress{2} - RelLightDelay{2}
	Diff3 = RelLightDelayRegress{3} - RelLightDelay{3}

	if any( abs(Diff1) > tol ) | any(isnan(Diff1)) | ...
	   any( abs(Diff2) > tol ) | any(isnan(Diff2)) | ...
	   any( abs(Diff3) > tol ) | any(isnan(Diff3)) 
		failed = 1;
	end
else
	failed = 1;
end

end