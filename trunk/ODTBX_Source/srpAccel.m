function  accel = srpAccel(Epoch, DCi2b, SC)
   
% SRPACCEL Compute solar radiation pressure acceleration
%
%  accel = srpAccel(Epoch, DCi2b, SC)
%
%  This function computes Solar radiation pressure acceleration on a spacecraft
%  using a single or multiple flat plate spacecaft model and ECI to body
%  rotation matrix.  The body frame normal and the specular and diffuse
%  reflectivity for each plate in the model must be specified.  Shadowing from
%  planetary bodies (currently just Earth and/or the moon) can be specified.
%
%  The penumbra shadow increases as the spacecraft gets further from the planet.
%  The penumbra effect is aproximated by producing a 50% acceleration in this
%  region.  The penumbra model is turned off when the spacecraft is outside of
%  the maximum distance of the umbra.  A spacecraft that resides in the umbra
%  of one body and the penumbra of a second body will receive the full
%  shadow of the umbra.
%
%  Note that the sum of the specular and diffuse coefficients must not
%  exceed 1.0.  For a real surface, the sum will be slightly less than 1.0.
%
%  The surface is assumed to be a single flat plate.  If the angle between
%  norm_ECI and the spacecraft to Sun vector is greater than 90.0 degrees,
%  the force will be set to zero.
%
% INPUTS 
%    VARIABLE           SIZE      DESCRIPTION (Optional/Default)
%    Epoch              (1x1)     UTC time in Matlab datenum time format
%    DCi2b              (3x3)     ECI J2000 to Body direction cosine matrix
%    SC.satPos          (3x1)     ECI satellite coordinates (km)
%    SC.mass            (1x1)     Spacecraft mass (kg)
%    SC.area            (1xN)     Area of each plate (m)
%    SC.NormVec_B       (3xN)     Body frame normal vector for each plate
%    SC.c_specular      (1xN)     Specular reflectivity coefficient (0.0-1.0)
%                                   summed with c_diffuse <= 1.0
%    SC.c_diffuse       (1xN)     Diffuse reflectivity coefficient (0.0-1.0)
%                                   summed with c_specular <= 1.0
%    SC.ShadowBodies  (1xM Cell)  (optional) Shadowing planetary bodies
%                                    Empty or 'None' = no shadowing
%                                    'Earth' and/or 'Moon'
%                                   
% OUTPUT 
%    accel              (3x1)     Solar radiation pressure acceleration (m/s^2)
%
% References:
% 1. Enhanced Radiative Force Modeling of the Tracking and Data Relay Satellites
%    S.B. Luthcke, J.A. Marshall, et al
%    Journal of the Astronautical Sciences, Vol. 45, No. 3, July-September 1997
%       (19970701), pp. 349-370
%    The vector n (surface normal vector) is the surface normal of the surface
%    specified by the panel card.  If the angle between n and the spacecraft to
%    Sun vector (SC2Sun) is more than 90 degrees, there is no Solar radiation
%    pressure force; there could still be an "albedo/infared acceleration."  The
%    vector s (source incidence vector) is the projection of SC2Sun in the plane
%    of the panel.  The Solar radiation pressure force components would be
%    opposite the n and s vectors above per the equation.  I believe the s 
%    vector would be defined as follows:
%       n cross (SC2Sun cross n)
% 2. Solar Sailing, Colin R. McInnes, 1999, 2.3 (Physics of Radiation Pressure),
%    2.6 (Solar Sail Force Models)
% 3. Wertz, 1984, pg. 64, pg. 130
%
% EXAMPLE
%
%   >>Epoch = datenum('June 8, 2008'); %UTC
%   >>DCi2b = eye(3);
%   >>SC.satPos = [-19134000   -25513000     1275600]';
%   >>SC.mass = 1000;
%   >>SC.area = 5;
%   >>SC.NormVec_B = [1 0 0]';
%   >>SC.c_specular = 0.5;
%   >>SC.c_diffuse = 0.5;
%   >>SC.ShadowBodies = {'earth','moon'};
%
%   >>accel = srpAccel(Epoch, DCi2b, SC)
%   accel =
%      -3.424e-09
%     -1.8179e-09
%     -5.4831e-10
%
% VALIDATION TEST
%
%   To perform a validation test, replace the Epoch input with
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
% REGRESSION TEST
%
%   To perform a regression test, replace the Epoch input with 
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
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman      06/8/2008       Original
%   Brent Wm. Barbee    04/21/2009      Corrected the spelling of a couple
%                                       of words in the 'help' section
%                                       comments. Also transposed SC.satPos
%                                       in the EXAMPLE so that it runs
%                                       without error.
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy 

% Determine whether this is an actual call to the program or a test

if strcmpi(Epoch,'ValidationTest')

	accel = srpAccel_validation_test();

elseif strcmpi(Epoch,'RegressionTest')

	accel = srpAccel_regression_test();

else

	accel = getsrpAccel(Epoch,DCi2b,SC);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the acceleration as the sum of the individual plate forces divided by the total mass

function accel = getsrpAccel(Epoch,DCi2b,SC)

for k = 1:length(SC.area)
	if ~isfield(SC,'ShadowBodies')
		SC.ShadowBodies = {'NONE'};
	end

	SRPForce_ECI(:,k) = getsrpForce(Epoch,DCi2b,SC.NormVec_B(:,k),SC.area(k),SC.satPos,...
	                                SC.c_specular(k),SC.c_diffuse(k),SC.ShadowBodies);
end

if length(SC.area) > 1
	TotalForce = sum(SRPForce_ECI')';  % N (kg*m/s^2)
else
	TotalForce = SRPForce_ECI;
end

accel = TotalForce / SC.mass;      % m/s^2

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the force for an individual plate

function force = getsrpForce(Epoch,DCi2b,NormVec_B,area,satPos,c_specular,c_diffuse,ShadowBodies)

% testing (t_out)
% t_out=t

if c_specular + c_diffuse > 1.0
	warning('SC.c_specular + SC.c_diffuse must be <= 1')
end

aunit = JATConstant('au')/1e3; % Astronomical Unit (km)
c_light = JATConstant('c'); % Speed of light (m/s)
flux_Sun = 1367;   % W/m^2  Space Mission Analysis and Design - Wertz, Larson 1999

% compute solar e/raphemeris

SunPos = ephemDE405('GEOCENTRIC_SUN',Epoch,'UTC');  % Sun eci position (km)
rsun = SunPos - satPos;

% Determine if SpaceCraft is in shadow (c_shad = 0 means shadow)

c_shad = getShadow(Epoch,satPos,ShadowBodies,aunit);

% compute vector from the satellite (SC) to the sun unit vector (rs2su)

rs2su = unit(rsun - satPos);

% compute vector from Sun to satellite unit vector (rSun2SCu)

rSun2SCu = norm(satPos - rsun);

% compute plate normal vector in ECI frame norm_ECI

norm_ECI = inv(DCi2b) * unit(NormVec_B);

% The surface is assumed to be a single flat plate.  If the angle (ang1_deg) between norm_ECI and
% the spacecraft to Sun vector (rs2su) is greater than 90.0 degrees, the force will be set to zero.
% determine c_ang (=1 for ang1 <= pi/2 rad, =0 for ang1 > pi/2 deg)

ang1_cos = dot(norm_ECI,rs2su);
ang1 = acos(ang1_cos);

if (ang1<=pi/2)
     	c_ang=1;
else
     	c_ang=0;
end

% compute force vector

% fsrp_1: force term (scalar) in N (kg*m/s^2)

fsrp_1 = c_shad * c_ang * (flux_Sun*area/(c_light)) * ang1_cos * ((aunit/rSun2SCu)^2);

% Code prior to addition of scale factors:
% fsrp_2: force term (vector) in direction of norm_ECI
% fsrp_3: force term (vector) in direction of v2u
fsrp_2 = 2*((c_diffuse/3)+(c_specular*ang1_cos))*norm_ECI;
fsrp_3 = (1-c_specular)*rs2su;

% Calculate total force in N (kg*m/s^2)
force = -fsrp_1 * (fsrp_2 + fsrp_3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine whether any of the specified celestial bodies is shadowing the spacecraft

function [c_shad] = getShadow(Epoch,satPos,ShadowBodies,aunit)

c_shad = 1;

if ~strcmpi(ShadowBodies(1),'NONE')

	for k = 1:length(ShadowBodies)

		Body = char(ShadowBodies(k));

		if strcmpi(Body,'Earth') | strcmpi(Body,'Moon')


			% Determine sun vector relative to body (km)

			rsun_b = ephemDE405('SUN',Epoch,'UTC') - ephemDE405(Body,Epoch,'UTC');

			% Determine satellite position relative to body (km)

			satPos_b = satPos + ephemDE405('EARTH',Epoch,'UTC') - ephemDE405(Body,Epoch,'UTC');

			% planet centric distance of the sun (km)

			rmsun = norm(rsun_b);

			% planet centric unit vector of the sun

			usun = unit(rsun_b);

			% planet centric distance of the spacecraft (km)

			rscm = norm(satPos_b);

			% compute unit position vector of spacecraft

			usat = unit(satPos_b);

			% determine planetary body equitorial radius

			req = JATConstant('meanRadius',Body)/1e3;

			% determine shadow conditions

			a = usat(2) * usun(3) - usat(3) * usun(2);
			b = usat(3) * usun(1) - usat(1) * usun(3);
			c = usat(1) * usun(2) - usat(2) * usun(1);

			d = norm([a  b c]);

			e = dot(usat, usun);

			u = asin(0.00460983743 / (rmsun / aunit));
			p = asin(0.0046951089 / (rmsun / aunit));

			if (e > 0)
				q = -d;
			else
     				q = d;
			end

			if strcmpi(upper(Body),'EARTH')

				% increase the Earth's radius by 90 km
				% to account for the atmosphere

				ratm = req + 90;
			else

				% no atmosphere information is available
				% for other celestial bodies

				ratm = req;

			end

			PenCutoff = ratm/tan(u);
	
			b = asin(ratm / rscm);

			v = b - u;
			w = b + p;

			x = sin(v);
			y = sin(w);
	
			% determine shadow conditions
			% determine c_shad (=0 for in shadow, =1 for in Sunlight)

			if (q <= y & q > x) & c_shad == 1 & rscm < PenCutoff
	     			% penumbra
	     			c_shad = 0.5;
			elseif (q <= x & q >= 0)
	     			% umbra
	     			c_shad = 0;
			else
     				% sunlight
			end

		else

			warning(['Shadow model not available for ',Body])

		end

	end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% srpAccel_validation_test - Validation for srpAccel.m

function failed = srpAccel_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/5/2008	 	Original

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-7;
format compact

%%%%%%%%%
Case = 1; % Earth orbit - no shadow - plate normal in sun direction
disp('Case = 1;  Earth orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 2; % Earth orbit - no shadow - 2 plates normal in sun direction
disp('Case = 2;  Earth orbit - no shadow - 2 plates normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = [1.0,1.0];
SC(Case).c_diffuse = [0.5,0.5];
SC(Case).NormVec_B = [1 0 0;1 0 0]';
SC(Case).area = [5,5];
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 3; % Saturn orbit - no shadow - plate normal in sun direction
disp('Case = 3;  Saturn orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = -ephemDE405('earth',Epoch(Case),'UTC') + ephemDE405('saturn',Epoch(Case),'UTC')
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

RsunUnit = unit(-SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC'));
CrossProd = cross(SC(Case).NormVec_B,RsunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,RsunUnit)) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 4; % Earth orbit - no shadow - plate normal edge-on to sun direction
disp('Case = 4;  Earth orbit - no shadow - plate normal edge-on to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/2,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 5; % Earth orbit - no shadow - plate normal 45 deg to sun direction
disp('Case = 5;  Earth orbit - no shadow - plate normal 45 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/4,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 6; % Earth orbit - shadow - plate normal to sun direction
disp('Case = 6;  Earth orbit - shadow - plate normal to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
CrossProd = cross([0 1 0]',rSunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
DCecef2sun = inv(dcm('axs',RotAngle+pi,RotVector));

SC(Case).satPos = DCecef2sun*[0 2 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 7; % Earth orbit - no shadow - plate normal 135 deg to sun direction
disp('Case = 7;  Earth orbit - no shadow - plate normal 135 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+3*pi/4,RotVector+3*pi/4);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%% Validation Tests %%%%

if exist('srpAccel_ValidationData6_08.mat') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')


	%%%%%%%%%
	Case = 1;

	JatSRPObject = jat.forces.SolarRadiationPressure(SC(Case).mass,SC(Case).area,...
                                                 1+SC(Case).c_specular+2/3*SC(Case).c_diffuse);
	Rsun = ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC');
	JatsatPos = jat.matvec.data.VectorN(SC(Case).satPos*1e3);
	JatRsun = jat.matvec.data.VectorN(Rsun*1e3);
	PredictedAccel(:,Case) = JatSRPObject.accelSRP(JatsatPos,JatRsun).x;

	%%%%%%%%%
	Case = 2;
	PredictedAccel(:,Case) = 2*PredictedAccel(:,1);

	%%%%%%%%%
	Case = 3;

	JatSRPObject = jat.forces.SolarRadiationPressure(SC(Case).mass,SC(Case).area,...
                                               	1+SC(Case).c_specular+2/3*SC(Case).c_diffuse);

	Rsun = -SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC');
	JatsatPos = jat.matvec.data.VectorN(SC(Case).satPos*1e3);
	JatRsun = jat.matvec.data.VectorN(Rsun*1e3);
	PredictedAccel(:,Case) = JatSRPObject.accelSRP(JatsatPos,JatRsun).x;

	%%%%%%%%%
	Case = 4;
	PredictedAccel(:,Case) = [0 0 0]';

	%%%%%%%%%
	Case = 5;
	PredictedAccel(:,Case) = inv(DCi2b{Case}) * (-SC(Case).NormVec_B) * ...
	                       norm(PredictedAccel(:,1))*cos(pi/4);

	%%%%%%%%%
	Case = 6;
	PredictedAccel(:,Case) = [0 0 0]';

	%%%%%%%%%
	Case = 7;
	PredictedAccel(:,Case) = [0 0 0]';


	%save srpAccel_ValidationData6_08 PredictedAccel
	clear PredictedAccel % get data from data file rather than this instance

	% Load validation values PredictedAccel
	load srpAccel_ValidationData6_08

	disp(' ')
	disp('Test Results (m/s^2)')
	disp(' ')
	disp('        Test             Predicted')
	disp('----------------------------------------')

	for k = 1:7
		fprintf('\nCase %d\n\n',k)
		disp([Accel(:,k),PredictedAccel(:,k)])
	end

	disp(' ')


	Diff = [PredictedAccel - Accel]

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end


	%%%%%%%%%
	Case = 8; % Test Earth shadowing
	disp(' ')
	disp('Case 8; Testing Earth shadowing ...')
	disp(' ')

	Re = JATConstant('rEarth')/1e3;
	Xgrid = [-3:0.1:3] * Re;
	Ygrid = [5:-10:-300] * Re;
	ShadowBodies = {'EARTH'};
	Epoch(Case) = datenum('May 16, 2008'); %UTC
	aunit = JATConstant('au')/1e3; % Astronomical Unit (km)

	rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
	CrossProd = cross([0 1 0]',rSunUnit);
	RotVector = unit(CrossProd);
	RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
	DCecef2sun = inv(dcm('axs',RotAngle,RotVector));

	c_shad_totalx = zeros(length(Xgrid),length(Ygrid));
	c_shad_totalz = c_shad_totalx;

	for y = 1:length(Ygrid)
		for x = 1:length(Xgrid)
			satPosECEF = [Xgrid(x);Ygrid(y);0];
			satPos = DCecef2sun * satPosECEF;
			c_shad_totalx(x,y) = getShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
		for z = 1:length(Xgrid)
			satPosECEF = [0;Ygrid(y);Xgrid(z)];
			satPos = DCecef2sun * satPosECEF;
			c_shad_totalz(z,y) = getShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
	end

	figure(1)
		subplot(211)
			pcolor(-Ygrid/2/Re,Xgrid/2/Re,c_shad_totalx)
			title(['Case 8 - Earth Shadow Plot Varied in 2 directions.  ',...
	       		'Red = Sun, Green = Penumbra, Blue = Umbra'])
			ylabel('X-axis (Earth Diameters)')
		subplot(212)
			pcolor(-Ygrid/2/Re,Xgrid/2/Re,c_shad_totalz)
			ylabel('Z-axis (Earth Diameters)')
			xlabel('Y-axis (Earth Diameters)')

	%%%%%%%%%
	Case = 9; % Test Neptune shadowing
	% this should produce a warning since no neptune shadowing is currently available
	disp(' ')
	disp('Case 9; Testing Neptune shadowing ...')
	disp(' ')
    
    fprintf(1, '\n');
    fprintf(1, 'Note that this should produce a warning since Neptune\n');
    fprintf(1, 'is not currently available.\n');
    fprintf(1, '\n');

	Epoch(Case) = datenum('May 16, 2008'); %UTC
	Rn = JATConstant('meanRadius','neptune')/1e3;
	NeptunePos = ephemDE405('neptune',Epoch(Case),'UTC') - ephemDE405('EARTH',Epoch(Case),'UTC');


%	Xgrid = [-2:0.1:2] * Rn;
%	Ygrid = [5:-500:-8000] * Rn;
	Xgrid = [0.1] * Rn;
	Ygrid = [-500] * Rn;
	ShadowBodies = {'neptune'};
	aunit = JATConstant('au')/1e3; % Astronomical Unit (km)

	rSunUnit = unit(ephemDE405('SUN',Epoch(Case),'UTC') - ephemDE405('neptune',Epoch(Case),'UTC'));
	CrossProd = cross([0 1 0]',rSunUnit);
	RotVector = unit(CrossProd);
	RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
	DCecef2sun = inv(dcm('axs',RotAngle,RotVector));

	c_shad_totalx = zeros(length(Xgrid),length(Ygrid));
	c_shad_totalz = c_shad_totalx;

	for y = 1:length(Ygrid)
		for x = 1:length(Xgrid)
			satPosECEF = [Xgrid(x);Ygrid(y);0];
			satPos = DCecef2sun * satPosECEF + NeptunePos;
			c_shad_totalx(x,y) = getShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
%		for z = 1:length(Xgrid)
%			satPosECEF = [0;Ygrid(y);Xgrid(z)];
%			satPos = DCecef2sun * satPosECEF + NeptunePos;
%			c_shad_totalz(z,y) = getShadow(Epoch(Case),satPos,ShadowBodies,aunit);
%		end
	end

%	figure(2)
%		subplot(211)
%			pcolor(-Ygrid/2/Rn,Xgrid/2/Rn,c_shad_totalx)
%			title(['Case 9 - Neptune Shadow Plot Varied in 2 directions.  ',...
%			       'Red = Sun, Green = Penumbra, Blue = Umbra'])
%			ylabel('X-axis (Neptune Diameters)')
%		subplot(212)
%			pcolor(-Ygrid/2/Rn,Xgrid/2/Rn,c_shad_totalz)
%			ylabel('Z-axis (Neptune Diameters)')
%			xlabel('Y-axis (Neptune Diameters)')


else
	failed = 1;
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% srpAccel_regression_test - Validation for srpAccel.m

function failed = srpAccel_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/8/2008	 	Original

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-12;
format compact

%%%%%%%%%
Case = 1; % Earth orbit - no shadow - plate normal in sun direction
disp('Case = 1;  Earth orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

%%%%%%%%%
Case = 2; % Earth orbit - no shadow - 2 plates normal in sun direction
disp('Case = 2;  Earth orbit - no shadow - 2 plates normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = [1.0,1.0];
SC(Case).c_diffuse = [0.5,0.5];
SC(Case).NormVec_B = [1 0 0;1 0 0]';
SC(Case).area = [5,5];
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 3; % Saturn orbit - no shadow - plate normal in sun direction
disp('Case = 3;  Saturn orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = -ephemDE405('earth',Epoch(Case),'UTC') + ephemDE405('saturn',Epoch(Case),'UTC')
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

RsunUnit = unit(-SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC'));
CrossProd = cross(SC(Case).NormVec_B,RsunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,RsunUnit)) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 4; % Earth orbit - no shadow - plate normal edge-on to sun direction
disp('Case = 4;  Earth orbit - no shadow - plate normal edge-on to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/2,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 5; % Earth orbit - no shadow - plate normal 45 deg to sun direction
disp('Case = 5;  Earth orbit - no shadow - plate normal 45 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/4,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 6; % Earth orbit - shadow - plate normal to sun direction
disp('Case = 6;  Earth orbit - shadow - plate normal to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
CrossProd = cross([0 1 0]',rSunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
DCecef2sun = inv(dcm('axs',RotAngle+pi,RotVector));

SC(Case).satPos = DCecef2sun*[0 2 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 7; % Earth orbit - no shadow - plate normal 135 deg to sun direction
disp('Case = 7;  Earth orbit - no shadow - plate normal 135 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+3*pi/4,RotVector+3*pi/4);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

disp(' ')

%PreviousAccel = Accel;
%save srpAccel_RegressionData6_08 PreviousAccel

%%%% Validation Tests %%%%

if exist('srpAccel_RegressionData6_08.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')


	% Load validation values PreviousAccel
	load srpAccel_RegressionData6_08

	disp(' ')
	disp('Test Results (m/s^2)')
	disp(' ')
	disp('        Test             Previous')
	disp('----------------------------------------')

	for k = 1:7
		fprintf('\nCase %d\n\n',k)
		disp([Accel(:,k),PreviousAccel(:,k)])
	end

	disp(' ')


	Diff = [PreviousAccel - Accel]

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayInOut(Case,Epoch,DCi2b,SC,Accel)

disp(' ')
disp(['Case = ',num2str(Case)])
disp(['Epoch = ',num2str(Epoch),' UTC'])
DCi2b
fprintf('SC.satPos = \n')
disp(SC.satPos)
disp(['SC.c_specular = ',num2str(SC.c_specular)])
disp(['SC.c_diffuse = ',num2str(SC.c_diffuse)])
fprintf('SC.NormVec_B = \n')
disp(SC.NormVec_B)
disp(['SC.area = ',num2str(SC.area)])
disp(['SC.mass = ',num2str(SC.mass)])
fprintf('SC.ShadowBodies = \n')
disp(SC.ShadowBodies)
fprintf('Accel = \n')
disp(Accel)
disp(' ')

end
