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
% VALIDATION/REGRESSION TEST
%
%   These tests were moved to srpAccel_test.m to conform to the new
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
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman      06/8/2008       Original
%   Brent Wm. Barbee    04/21/2009      Corrected the spelling of a couple
%                                       of words in the 'help' section
%                                       comments. Also transposed SC.satPos
%                                       in the EXAMPLE so that it runs
%                                       without error.
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy 
%   Ravi Mathur         08/28/2012      Extracted regression test


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the acceleration as the sum of the individual plate forces divided by the total mass

for k = 1:length(SC.area)
	if ~isfield(SC,'ShadowBodies')
		SC.ShadowBodies = {'NONE'};
	end

	SRPForce_ECI(:,k) = getSRPForce(Epoch,DCi2b,SC.NormVec_B(:,k),SC.area(k),SC.satPos,...
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

function force = getSRPForce(Epoch,DCi2b,NormVec_B,area,satPos,c_specular,c_diffuse,ShadowBodies)

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

c_shad = getBodyShadow(Epoch,satPos,ShadowBodies,aunit);

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