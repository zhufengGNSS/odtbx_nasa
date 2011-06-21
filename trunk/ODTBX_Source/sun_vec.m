%SUN_VEC Computes unit vector from origin of ECI frame to sun
%
%  [sun_equ] = sun_vec(start_day)
%
%  Inputs:    start_day    [1,n] days since 00:00:00 01/06/1980
%
%  Outputs:   sun_equ      [3,n] Earth-Sun vector in an ECI frame
%
%  Given start_day measured in days since 01/06/80, computes true 
%  longitude of the Sun.  Earth-Sun vector is then computed in the
%  ecliptic plane, and is rotated into the ECI frame (equatorial 
%  plane).  # ephemeris days assumed equal to # Julian days
%  The Earth-Sun vector rotates approximately 1.02 deg/day
%
%  Reference: J.R. Wertz, p 141
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

function [sun_equ] = sun_vec(start_day)

%  Julian days since Jan 0,1900
%  Reference for this calculation is JD 2,415,020 which 
%  corresponds to 12:00:00 Jan 0,1900 ET (or 12:00:00 Dec 31,1899)
jd = 29224.5 + start_day;

%  Mean longitude of sun, measured in the ecliptic from mean 
%  equinox of date:
L = (279.696678 + 0.9856473354.*jd + 2.267e-13.*jd.^2);

%  Mean anomaly of sun in radians
Ms_r = (pi/180)*(358.475845 + 0.985600267.*jd - (1.12e-13).*jd.^2 - (7e-20).*jd.^3);

%  Correction between mean longitude and true longitude
dL = 1.918.*sin(Ms_r) + 0.02.*sin(2.*Ms_r);

%  True longitude of sun, in radians
L_sun = rem((pi/180)*(L+dL),2*pi);

%  Compute sun unit vector in ECI frame, where the Earth's 
%  equatorial plane is inclined inc_E radians to the ecliptic
%  R defines a rotation about the x-axis
inc_E = (pi/180)*(-23.45);
R = [1,0,0; 0,cos(inc_E),sin(inc_E); 0,-sin(inc_E),cos(inc_E)];% [3,3]
sun_ecl = [cos(L_sun);sin(L_sun);zeros(1,size(start_day,2))];  % [3,n]
%  Since R is constant through time, can do a simple matrix multiply
sun_equ = R*sun_ecl;   % [3,n]


