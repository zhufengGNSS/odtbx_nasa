function [pos,vel,vel_tot] = alm2xyz(time,alm,frame)
% ALM2XYZ  Compute R and V from GPS alamanc data
%
%    [pos,vel] = alm2xyz(time,alm,frame) 
%    Function computes position and velocity of satellite(s) from 
%    GPS almanac data, corresponding to input time vector.
%
%    INPUTS:
%       time(1,n)    vector, in linear GPS seconds, since start of 
%                    GPS time
%       alm(mm,13)   matrix, almanac data in modified format
%			               mm = #PRNs, 13 variables each PRN
%       frame		   scaler, output coord frame
%	                     0-ECI, 1-ECEF
%    OUTPUTS:
%       Size of output matrixes will be [3,n,mm]
%         If mm=1, size of outputs will be [3,n]
%       pos	      x,y,z position of each SVN at each 'n' time step
%       vel	      x,y,z velocity of each SVN at each 'n' time step
%       vel_tot   x,y,z total velocity of each SVN at each 'n' time step
%                 If frame = 0 (ECI), vel = vel_tot.
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
%   Mike Moreau         1997               Original
%   Mike Moreau         01/02/2002         Updated
%   Mike Moreau         6/2004             Changed GPS vectors to 3D arrays
%   Kevin Berry         09/21/2007         Modified for ODTBX
%   Keith Speckman      03/15/2008         Minor comment change
%   Brent Wm. Barbee    07/27/2009         Added check & warning for
%                                          propagation time > 1 day
%   Brent Wm. Barbee    07/30/2009         Updated the ndays calculation as
%                                          suggested by Sun.
%   Brent Wm. Barbee    08/03/2009         Modified warning message.
%   Brent Wm. Barbee    08/14/2009         Changed warning message
%                                          threshold from 1 day to 7 days
%                                          and added MSGID to the warning.

%*********************************************************
% Almanac file has the following format:
% Each row corresponds to one SVN
% Each SVN has 13 pieces of data associated with it:
% column   Description                       variable name
% 1        ID                                SVN
% 2        SV health                         health
% 3        eccentricity                      ecc
% 4        time of applicability [s]         toe
% 5        orbital inclination (rad)         inc     
% 6        rate of right ascension (rad/s)   OMEGA_dot
% 7        SQRT(a) (m^1/2)                   root_a
% 8        right ascension at toe (rads)     OMEGA_not
% 9        argument of perigee (rad)         arg_peri
% 10       mean anomaly (rad)                M_not
% 11       Af0 (s)                           clock bias
% 12       Af1 (s/s)                         clock drift
% 13       GPS week                          week
%********************************************************

%% These values need to obtained and called in multiple subfunctions somehow
MU = 3.986005e5;              % Earth's grav. parameter [km^3/s^2] from GPS ICD 1991
EARTH_RATE = 7.2921151467e-5;  % WGS-84 Earth's rotation rate [rad/s] from GPS ICD 1991

%%
mm = size(alm,1);                   % number of PRNs

% t_epoch is equal to zero at the start of the week
% Can be used to directly compute the delta time between the current
% time step and the almanac epoch, even across week boundaries
t_epoch = time-alm(1,13)*7*86400;

% Correct OMEGA provided in almanac for ECI                                    
if frame == 0   % ECI
    disp('ECI to ECEF conversion isn''t implimented yet')
%    alm(:,8) = mod(alm(:,8) + sidereal(gps2utc([alm(:,13),zeros(size(alm(:,13)))])),2*pi)
   % should be able to replace sidereal part with:
%    testt=datenum(gps2utc([alm(1,13),zeros(size(alm(1,13)))]));
%    jatEarthRef.GAST(matlabTime2MJD(testt),convertTime('TT', 'UTC', testt))
   % or:
%    jatEarthRef.GMST(matlabTime2MJD(testt))
   % the problem is that the input is UT1 but I only have UTC
end

% Read almanac parameters (orbital elements) from almanac
root_a = alm(:,7);
a = (root_a).^2;              % semimajor axis (kilometers)
toe = alm(:,4);               % time of reference ephemeris
M_not = alm(:,10);            % mean anomaly at toe
ecc = alm(:,3);               % eccentricity
arg_peri = alm(:,9);          % argument of periapsis
inc = alm(:,5);               % inclination
OMEGA_not = alm(:,8);         % long of asc node at epoch
OMEGA_dot = alm(:,6);         % rate of long of asc node
n_not = sqrt(MU./a.^3);       % mean motion (rad/s)

% Initialize variables
pos = zeros(3,size(t_epoch,2),mm);
vel = zeros(3,size(t_epoch,2),mm);
vel_tot = vel;

ndays = max([max(abs(t_epoch(1)-toe)) max(abs(t_epoch(end)-toe))])/86400;

if(ndays > 7)
   
    warning('ODTBX:alm2xyz:ephEpoch', 'The epoch for the ephemeris you are using is off by %0.3f days; you may want to use a different one, if available.', ndays);
    
end


for time_step = 1:size(t_epoch,2)      % Start loop through each time step
	% Compute t_k, elapsed time relative to epoch
	% The variable 't' is the current time converted to "seconds of current week"
	% The reference time 'toe' is in seconds of current week
	t = t_epoch(time_step);
	t_k = t - toe;
    
	M_k = M_not + (n_not.*t_k);   % mean anomaly at t
	M_k = rem(M_k,2*pi);          % Make 0<=M_k<=360
	
	% Iterate to find eccentric anomaly.
	% Use M_k as initial guess for E.
	% All PRNs iterated the same ammount, some may converge earlier
	% than others, but there is no logic to loof for that.
	E = M_k;                      % First guess for E
	E_old = E+1;                  % Set E-old for at least one loop
	while sum(abs(E - E_old) >= 0.0000000005) > 0
      E_old = E;
      E = M_k + (ecc.*sin(E_old));
	end
	
	% true anomaly as a function of eccentric anomaly
	v = atan2(((sqrt(1 - ecc.^2)).*sin(E)),(cos(E) - ecc));
	
	% For [x,y,z] in the ECEF, 'frame' makes 'rotate' the Earth's rotation rate.  
	% For [x,y,x] in ECI, 'frame' makes 'rotate' equal zero.
	rotate = EARTH_RATE*frame;
	
	% adjusted longitude of ascending node
	% OMEGA_not and OMEGA_dot from almanac
	OMEGA = OMEGA_not + ((OMEGA_dot-rotate).*t_k) - (rotate.*toe);
	
	phi = v + arg_peri;		      % argument of latitude
	R = a.*(1 - ecc.*cos(E));		% radius to satellite
	
	% satellite position
	x = R.*((cos(OMEGA).*cos(phi)) - (sin(OMEGA).*sin(phi).*cos(inc)));
	y = R.*((sin(OMEGA).*cos(phi)) + (cos(OMEGA).*sin(phi).*cos(inc)));
	z = R.*sin(inc).*sin(phi);
	
	p = a.*(1-ecc.^2);			   % parameter
	h = sqrt(MU*p);				   % specific angular momentum
	
	% satellite velocity
	vx = (h./R).*((x.*ecc.*sin(v))./p - (cos(OMEGA).*sin(phi) + sin(OMEGA).*cos(phi).*cos(inc)));
	vy = (h./R).*((y.*ecc.*sin(v))./p - (sin(OMEGA).*sin(phi) - cos(OMEGA).*cos(phi).*cos(inc)));
	vz = (h./R).*((z.*ecc.*sin(v))./p + sin(inc).*cos(phi));
	
	% Save position and velocity info for each satellite at this 
	% time step (3xmm).  Each column corresponds to one PRN. 
	sat_pos = [x';y';z'];
	sat_vel_int = [vx';vy';vz'];
	
	% Must subtract (EARTH_RATE X R) term from ECF velocity
	if frame % (if ECF)
      % Compute EARTH_RATE cross R in the ECF frame
      o_cross_r = cross([0 0 EARTH_RATE]'*ones(1,mm),sat_pos);
      % convert from ECI to ECEF velocity by subtracting the omega cross r term
      sat_vel = sat_vel_int - o_cross_r;
	else
      sat_vel = sat_vel_int;
	end
     
	% Reshape such that each column contains all info for all PRNs at 
	% one time step.  First 3 rows=PRN1, second 3 rows=PRN2, ...        start 3,mm
	% pos(i,:)=[x1,x2,x3,x4,...,znn
	%           y1,y2,y3,y4,...,ynn
	%           z1,z2,z3,z4,...,znn
	%           x1,x2,x3,x4,...,xnn]
	% vel ... same format
	% final matrices are ((3*mm),size(t_epoch,2))
	%pos(:,time_step)=[reshape(sat_pos,mm*3,1)];
	%vel(:,time_step)=[reshape(sat_vel,mm*3,1)];
	
	% Data saved in 3D arrays
	pos(:,time_step,:) = sat_pos;
	vel(:,time_step,:) = sat_vel;
    vel_tot(:,time_step,:) = sat_vel_int;
end % (time_step)
