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
%  VALIDATION/REGRESSION TEST
%
%   These tests have been moved to EarthOrbitPlot_test.m to conform to
%   the new regression testing format.
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
%  Ravi Mathur       08/28/2012         Extracted regression test

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