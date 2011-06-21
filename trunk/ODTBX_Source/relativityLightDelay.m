function [RelLightDelay] = relativityLightDelay(Ephem)

% RELATIVITYLIGHTDELAY  Returns light time delay error from relativistic effects
%
%  [RelLightDelay] = relativityLightDelay(Ephem)
%   
%   relativityLightDelay returns the relativity effects delay distance between
%   satellites and groundstations.  The equations account for "the reduction 
%   in the coordinate velocity of light below c and the bending of the light
%   path. The bending increases the path length but also increases the 
%   coordinate velocity of light because the curved light path is further away
%   from the gravitating body than the straight-line path. The net effect of
%   the bending is to decrease the light time by the increase in the path
%   length divided by c."  The resulting delay should be treated as a reduction
%   in the range measurement.
%
%   The Parameterized Post-Newtonian (PPN) parameter (gamma in the Moyer
%   document) was set to unity, which represents the value in general
%   relativity.  This was done to match the equations in the heritage code. 
%   The actual value of gamma is slightly above 0.997.
%
%   This code was modified from code contained in the dsn utilities developed
%   by William Campbell.  The equations were verified against the Moyer
%   document, which is referenced below.
%
% INPUTS 
%    VARIABLE           SIZE      DESCRIPTION (Optional/Default)
%    Ephem.satPos       (3xN)     ECEF or ECI satellite coordinates (m)
%    Ephem.SatCoords    string    (Optional) Satellite coordinate system
%                                   'ECEF' (Default) or 'ECI'
%    Ephem.Epoch        (1xN)     UTC time in Matlab datenum time format
%                                   (Optional - Required only for ECI SatCoords)
%    Ephem.StationInfo  string    Method of specifying ground station position
%
%   CASE: Ephem.StationInfo = 'ECEF'
%   ----------------------------------------------------------------------------
%    Ephem.staPos       (3xM)     Station ECEF coordinates (m)
%
%   CASE: Ephem.StationInfo = 'LatLonHeight'
%   ----------------------------------------------------------------------------
%    Ephem.lat          (1xM)     Station geodetic latitude (deg)
%    Ephem.lon          (1xM)     Station geodetic longitude (deg)
%    Ephem.height       (1xM)     Station geodetic height (km)
%
%  OUTPUTS 
%    RelLightDelay      (NxM)    Relativity delay distance between satellites
%                                  and ground stations (m)
%
% REFERENCE
%
%    Deep Space Communications and Navigation Series
%      Issued by the Deep Space Communications and Navigation Systems 
%      Center of Excellence 
%      Jet Propulsion Laboratory 
%      California Institute of Technology 
%      Joseph H. Yuen, Editor-in-Chief
%      (Moyer Document - October 2000)
%
%  VALIDATION TEST
%
%   To perform a validation test, replace the Ephem structure with 
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
%  REGRESSION TEST
%
%   To perform a regression test, replace the Ephem structure with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
%
%   keywords: delay, time, relativity, relativistic
%   See also: lightTimeCorrection
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

% REVISION HISTORY
%  Author               Date            Comment
%              	    (MM/DD/YYYY)
%  Keith Speckman    08/15/2008         Original
%  Sun Hur-Diaz      12/10/2008         Added a check for singularity
%                                       in the planetary effects
%  Kevin Berry       06/25/2009         Fixed time scale discrepancy and
%                                       updated the regression test

% Determine whether this is an actual call to the program or a test

if strcmpi(Ephem,'ValidationTest')

	RelLightDelay = relativityLightDelay_validation_test();

elseif strcmpi(Ephem,'RegressionTest')

	RelLightDelay = relativityLightDelay_regression_test();
else

        RelLightDelay = getDelay(Ephem);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function RelLightDelay = getDelay(Ephem);

% Determines upleg and downleg times of a signal sent to a spacecraft.
% Inputs:
% t3: time of reception (Julian days)            		1x1
% scX: spacecraft state in barycentric (AU & AU/day)       1x6

%%%% Define Constants

% Speed of light (AU/day)
c = JATConstant('c')/JATConstant('au')/JATConstant('sec2Days');

% Astronomical Unit (m)
au = JATConstant('au');

% Gravitational Constants for 9 planets, the moon, and the sun respectively (AU^3/day^2)
mu = [JATConstant('muMercury'); 
      JATConstant('muVenus');
      JATConstant('muEarth');
      JATConstant('muMars');
      JATConstant('muJupiter');
      JATConstant('muSaturn');
      JATConstant('muUranus');
      JATConstant('muNeptune');
      JATConstant('muPluto');
      JATConstant('muMoon');
      JATConstant('muSun')]/JATConstant('au')^3*86400^2;

% Astronomical bodies of Influence      
BodyNames = {'MERCURY'
             'VENUS'
             'EARTH'
             'MARS'
             'JUPITER'
             'SATURN'
             'URANUS'
             'NEPTUNE'
             'PLUTO'
             'MOON'
             'SUN'};

% Determine size of outputs (M,N)
N = size(Ephem.satPos,2);

switch Ephem.StationInfo

	% Station Coordinates given in ECEF Coordinates
	case 'ECEF'
		M = size(Ephem.staPos,2);
	case 'LatLonHeight'
		M = size(Ephem.lat,2);
end


% Initialize azimuth and elevation for quicker execution time
RelLightDelay = zeros(N,M);

for n = 1:N

	if isfield(Ephem,'SatCoords') & strcmpi(Ephem.SatCoords,'ECEF')

		% Determine ECI satPos
		satPosECI = jatDCM('ecef2eci',Ephem.Epoch(n))*Ephem.satPos(:,n);

	else

		% Determine ECI satPos
		satPosECI = Ephem.satPos(:,n);

	end

	% Earth Position (m)
	EarthPos = ephemDE405('EARTH',Ephem.Epoch(n),'UTC')*1e3;

	% Convert to Barycentric coordinates to match original code (m)
	satPosBary = satPosECI + EarthPos;

	for m = 1:M

		% Convert Lat Lon Height to ECEF Vector if necessary (m)
		if strcmpi(Ephem.StationInfo,'LatLonHeight') & n == 1
			Ephem.staPos(:,m)= LLA2ecef(Ephem.lat(m),Ephem.lon(m),Ephem.height(m))*1e3;
		end

		% ECI Position (m)
		stn_ECI = jatDCM('ecef2eci',Ephem.Epoch(n)) * Ephem.staPos(:,m);

		% Barycentric Position (m)
		stn_Bary = stn_ECI + EarthPos;

		% 1-2 SOLUTION (au)
		dt12_old = 0;
		r1_vec = stn_Bary/au;
		r2_vec = satPosBary/au;

		r12 = norm(r1_vec - r2_vec);

		% term1 is from heritage code and represents other range corrections
		term1 = 0;

		% Solar Relativistic Effect
		rsun = ephemDE405('Sun',Ephem.Epoch(n),'UTC')*1000/au;
		r1 = norm(r1_vec-rsun);
		r2 = norm(r2_vec-rsun);
		term2 = 2/c^3*mu(11)*...
				log((r1 + r2 + r12 + 2*mu(11)/c/c)/(r1 + r2 - r12 + 2*mu(11)/c/c));

		term3 = 0;


		% Planetary Relativistic Effects
		for j = 1:10
			rb_t2a = ephemDE405(BodyNames{j},Ephem.Epoch(n),'UTC')*1e3 /au;
			rb_t1 = ephemDE405(BodyNames{j},Ephem.Epoch(n)-dt12_old,'UTC')*1e3 /au;
			r1b = norm(r1_vec - rb_t1(1:3));
			r2b = norm(r2_vec - rb_t2a(1:3));
            if abs(r1b + r2b - r12)>eps
                term3 = term3 + 2/c/c/c*mu(j)*log((r1b + r2b + r12)/(r1b + r2b - r12));
            end
		end

		% The following comment statement is from the heritage code.  It is unknown whether
		% an imaginary component may result, but the "real" function was left in just in 
		% case.

		% A glitch was encountered, so real component taken
		RelLightDelay(n,m) = real(term1 + term2 + term3) * c * au;

	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation Test

function failed = relativityLightDelay_validation_test();

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

RelLightDelay{1} = getDelay(Ephem);
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

RelLightDelay{2} = getDelay(Ephem);
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

RelLightDelay{3} = getDelay(Ephem);

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

RelLightDelay{4} = getDelay(Ephem4);

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
RelLightDelayPlot = getDelay(Ephem)/JATConstant('c')*1e6*2;
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

function failed = relativityLightDelay_regression_test();

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

RelLightDelay{1} = getDelay(Ephem);
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

RelLightDelay{2} = getDelay(Ephem);
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

RelLightDelay{3} = getDelay(Ephem);

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



