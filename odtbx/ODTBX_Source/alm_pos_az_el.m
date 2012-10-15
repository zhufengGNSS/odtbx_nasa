% Calculates GPS SV antenna azimuth and elevation for a given receiver state.
%
% This function calculates satellite position, velocity, azimuth angle, and elevation
% angle for the user-defined time duration for the specified PRNs. The azimuths and
% elevations are then compared to a mask table to check satellite visibility. All
% DOPs for all visible satellites are calculated at each time. The yaw angle and rate
% of each PRN is also modeled and used to calculate azimuth.
%
% Inputs
%   epoch = (scalar datenum) time associated with start of simulation
%   time = (1xn) measurement times (secs from epoch)
%   rx_pos = (6xn) receiver spacecraft state [position;velocity] (km)
%   tx_pos = (3,n,m) GPS SV position (km)
%   tx_vel = (3,n,m) transmitter velocity (km/s)
%       where n=number of time steps
%       and m=number of SVs, at each time step
%
% Outputs
%   az = (m,n) azimuth of the satellite antenna in degrees from 0 to 359
%   el = (m,n) off-boresite angle in degrees
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

%*************************************************************************
% Copyright (c) Regents of the University of Colorado All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are
% met:
% 
% - Redistributions of source code must retain the above copyright notice, 
%   this list of conditions and the following disclaimer.
% - Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
% - Neither the name of the University of Colorado nor the names of its 
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%*************************************************************************

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

%% Original Comments and Information:
%***************************************************************************************
% Title: alm_pos_az_el
% Author: Lisa Reeh
% Date written: 1 Feb. 2002
% Date modified: 27 May 2003
% Copyright 2002 University of Colorado, Boulder
% This script calls the following internal functions/scripts:
%
% The function 'ecef2eci' is called to transform a vector in the ECEF frame to a
% vector in the ECI frame
%   input: vector in ECEF coordinates, Greenwich Sidereal Time in radians
%   output: vector in ECI coordinates
%
% The function 'alm_pos_vel' is called to solve Kepler's equation and
% compute the satellite's position and velocity
%   input: orbital elements from almanac file, current time, Greenwich Sideral time
%   output: 'sat_pos', the satellite's position vector in ECI (meters),
%           'sat_vel', the satellite's velocity vector in ECI (meters/sec)
%
% The function 'el_angle' is called to calculate the elevation angle of the
% satellite.
%   input: Local up vector at reference station, Line of Sight (LOS) unit vector
%   output: elevation angle of satellite in degrees
%
% The function 'eci2ecef' is called to transform a vector in the ECI frame to a
% vector in the ECEF frame
%   input: vector in ECI coordinates, Greenwich Sidereal Time in radians
%   output: vector in ECEF coordinates
%
% The function 'ECEF2NED' is calle dto transform the LOS unit vector from ECEF
% coordinates to NED (North, East, Down) coordinates.
%   input: Ref_LLA (lat., lon., alt. of reference station), LOS_unit (LOS unit vector)
%   output: LOS vector in NED coordinates
%
% The function 'az_angle' is called to calculate the azimuth angle of the satellite.
%   input: LOS_unit_NED (the line of sight unit vector from the reference station
%               to the satellite in the NED frame)
%	output: azimuth angle in degrees
%
% The function 'el_cut_off' is called to compare the satellite's elevation angle
% to a table of cut-off elevations to determine if the satellite is visible based on
% its azimuth.
%	input: az_angle (azimuth angle of the satellite, 0 to 359 degrees)
%	 		el_angle (elevation angle of the satellite, -90 to 90 degrees)
%           mask (structure containing cut-off elevations as a function of azimuth)
%	output: visible = 1 if satellite is visible (i.e., sat elev > cut-off elev.)
%			 visible = 0 if satellite is not visible
%
% The function 'get_beta_yaw_angles' calculates both the nominal and actual yaw angles
% and rates of the satellite based on the beta and alpha angles. Beta is the angle
% between the sun vector and the orbit plane; alpha is the angle between the sun
% vector and satellite position vector measured in the orbit plane.
%    input: sat_pos (ECI position vector of satellite in meters)
%           sat_vel (ECI velocity vector of satellite in meters/sec)
%           T (number of Julian centuries elapsed since J2000 epoch)
%           Theta_GST (Greenwich Sidereal Time, in radians)
%           j: index of current time in time vector
%           k: index of current PRN in PRN vector
%    output: beta angle (degrees), and global variables yaw, nominal yaw, yaw rate,
%               and nominal yaw rate (degrees and deg/sec)
%
% The function 'eci2ant' converts the ECI line of sight vector to local coordinates
% with respect to the satellite antenna
%    input: yaw - ideal yaw angle, degrees
%           sat_pos (ECI satellite position, in meters)
%           sat_vel (ECI satellite velocity, in meters/second)
%           LOS_unit (Line of sight unit vector, ECI frame, in meters)
%    output: LOS_ant (Line of sight unit vector, local antenna frame, in meters)
%
% The function 'get_axial_ratio' determines the satellite axial ratio as a function
% of the azimuth of the satellite antenna with respect to the line of sight vector
% and the off-boresite angle of the satellite.
%    input: sat_az (azimuth of the satellite, degrees)
%           OB_angle (off-boresite angle of the satellite, degrees)
%           sat_az_array (range of possible satellite azimuths, 0 to 359 degrees)
%           OB_angle_array (range of possible off-boresite angles, 0 to 15 degrees)
%           AR_table (table of axial ratios corresponding to each possible combination
%               of satellite antenna azimuths and off-boresite angles)
%    output: AR (axial ratio of the satellite)
%
% Finally, GDOP, HDOP, VDOP, PDOP, and TDOP are calculated for the visible satellites
% at each time, and a time history is plotted.
%************************************************************************************


% original code, now no longer used:
%global ISNOON ISSHADOW POSTSHADOW time yaw yaw_rate nom_yaw nom_yaw_rate

function [az,el,yaw] = alm_pos_az_el(epoch,time,rx_pos,tx_pos, tx_vel)



%% Get the receiver's position in the ECI reference frame

%Convert everything to meters
tx_pos = tx_pos*1000;
tx_vel = tx_vel*1000;
rx_pos = rx_pos*1000;


% Preallocate the required variables
rx_ECI_pos = zeros(size(rx_pos));
tx_ECI_pos = zeros(size(tx_pos));
tx_ECI_vel = zeros(size(tx_pos));
tx_ECI_pos_unit = zeros(size(tx_ECI_pos));
LOS_unit = zeros(size(tx_pos));

% Need a loop here, not sure how to do it otherwise

% Convert the receiver states from ECEF to ECI
w = [0; 0; JATConstant('wEarth')];
ecef2eci  = jatDCM('ecef2eci', (time/86400)+epoch);

for n=1:length(time)
  
   rx_ECI_pos(1:3,n) = ecef2eci(:,:,n)*rx_pos(1:3,n);
    
   % original code, now no longer used:
    %rx_ECI_pos = ecef2eci(gpsTime,rx_pos');
    %load tx_ECI_pos;
    %load tx_ECI_vel;
    
    for j = 1:1:size(tx_pos,3)
        
        % original code, now no longer used:
        %get the state into the correct form
        %tx_ecef_state = [tx_pos(:,:,j); tx_vel(:,:,j)];
        
        %for each transmitter, get the ECI position and velocity
        %[tmp_pos tmp_vel ] = ecef2eci([gpsTime(n,1) gpsTime(n,2)],tx_pos(1:3,n,j)',tx_vel(1:3,n,j)');
        
        %tx_ECI_pos(1:3,n,j) = tmp_pos';
        %tx_ECI_vel(1:3,n,j) = tmp_vel';
        
        tx_ECI_pos(1:3,n,j) = ecef2eci(:,:,n)*tx_pos(1:3,n,j);
        tx_ECI_vel(1:3,n,j) = ecef2eci(:,:,n)*(cross(w,tx_pos(1:3,n))+tx_vel(1:3,n));
        
        
        range = tx_ECI_pos(:,n,j)-rx_ECI_pos(1:3,n);
        LOS_unit(:,n,j) = range/norm(range);
        tx_ECI_pos_unit(:,n,j) = tx_ECI_pos(:,n,j)/norm(tx_ECI_pos(:,n,j));
        
        
    end
    
end

% original code, now no longer used:
%for n=1:length(time)
%    rx_ECI_pos(:,n) = ecef2eci(:,:,n)*rx_pos(:,n);
    
    % Now get the transmitters
%    for j = 1:1:32
%        tx_ECI_pos(:,n,j)  = ecef2eci(:,:,n)*tx_pos(:,n,j);
%        tx_ECI_vel(:,n,j)  = ecef2eci(:,:,n)*(cross(w,tx_pos(:,n,j))+tx_vel(:,n,j));
%        tx_ECI_pos_unit(:,n,j) = tx_ECI_pos(:,n,j)/norm(tx_ECI_pos(:,n,j));
%    end
%end



%% Need the elapsed julian centuries for the sun position computation

% Compute the initial Julian date
%dVec = datevec( convertTime('UTC','GPS',time/86400));
dVec = datevec(epoch+time/86400);
julianTime = 367*dVec(:,1) - floor(7*(dVec(:,1) + floor((dVec(:,2)+9)/12))/4) +...
    floor(275*dVec(:,2)/9) + dVec(:,3) + 1721013.5 + (((dVec(:,6)/60 + dVec(:,5))/60) + dVec(:,4))/24;

% Calculate the number of Julian centuries elapsed from the epoch J2000
T = (julianTime - 2451545)/36525;

% Calculate the beta angle and nominal yaw angle of the satellite in
% degrees.
[yaw] = get_beta_yaw_angles(tx_ECI_pos, tx_ECI_vel, T, time,tx_ECI_pos_unit);

% Transform the LOS vector to local coordinates with respect to the satellite antenna.
[LOS_ant] = eci2ant(yaw, tx_ECI_pos, tx_ECI_vel, LOS_unit);

% determine the off-boresite angle in degrees
el = (acos(dot(tx_ECI_pos_unit, LOS_unit))) * 180/pi;
if length(size(el)) > 2
    % only squeeze if we have 3 dimensions
    el = squeeze(el)';
    % correct the squeeze if needed:
    if size(el,1) ~= size(LOS_ant,3)
        el = el';
    end
end

% Compute the azimuth of the satellite antenna in degrees from 0 to 359
az = zeros(size(LOS_ant,3),size(LOS_ant,2));
for j = 1:1:size(LOS_ant,2) 
    for k = 1:1:size(LOS_ant,3)
        if (~isnan(LOS_ant(3,j,k)) && ~isnan(LOS_ant(2,j,k)))
            az(k,j) = (atan2(LOS_ant(3,j,k), LOS_ant(2,j,k))) * 180/pi;
        end
    end
end

az(az<0)=az(az<0)+360;

end


%% Function to compute the yaw angles
%*******************************************************************************
% Title: get_beta_yaw_angles
%
% Purpose: calculates the the vector from the center of the earth to the
% sun; the beta angle, which is the angle between the orbit plane of the
% satellite and the sun vector, and the yaw angle
%
% Input:
%       pos_vect: (3,n,m) position vector of the satellite (m, ECI)
%       vel_vect: (3,n,m) velocity vector of the satellite (m/s, ECI)
%       T: (1xn) number of elapsed Julian centuries (J2000 epoch)
%       time: (1xn) measurement times (secs from epoch)
%       tx_ECI_pos_unit: (3,n,m) unit vector of the satellite position
%       (ECI)
%           where n=number of time steps
%           and m=number of SVs, at each time step
%
% Output:
%       yaw: the yaw angle of the satellite in degrees
%
% Functions called:
%       sun_vector: computes the vector from the center of the Earth to the
%           center of the sun in meters in the ECI frame
%       get_nominal_yaw: computes the nominal yaw angle (degrees) for the
%           given PRN at the current time
%       compute_shadow: determines if the current PRN is in the earth's
%           penumbra at the current time (1 == yes, 0 == no)
%
% Reference: Bar-Sever, Yoaz E. "A new model for GPS yaw attitude," Journal
% of Geodesy (1996) 70:714-723
%
% Original function information:
% *************************************************************************
% Author: Lisa Reeh
% Date written: 8 Feb. 2002
% Date last modified: 27 May 2003
% Copyright 2002 University of Colorado, Boulder
%
% original parameter comments:
%       beta: angle in degrees
%       Theta: Greenwich Hour Angle in radians
%       time_tag: index of current time in time vector
%       PRN_tag: index of current PRN in PRN vector
%
% original Global Variables:
%       ISNOON: for each PRN, 0 (satellite is not in noon maneuver) or 1
%           (satellite is in noon maneuver)
%       ISSHADOW: for each PRN, 0 (satellite is not in shadow), 1 (satellite
%           is in shadow), or 2 (satellite is in post-shadow)
%       yaw, yaw_rate: actual yaw angle (degrees) and yaw rate (deg/sec)
%       nom_yaw, nom_yaw_rate: nominal yaw angle (degrees) and yaw rate
%           (deg/sec)
%       t_i, t_spin: shadow entry time and satellite spin up/down time
%       D: quantity used to determine the sign of the yaw rate during the
%           post-shadow maneuver
%*******************************************************************************

function yaw = get_beta_yaw_angles(pos_vect, vel_vect, T, time, tx_ECI_pos_unit)

% original code, now no longer used:
%global ISNOON ISSHADOW time dt yaw yaw_rate nom_yaw nom_yaw_rate
%global t_i t_spin shadow_entry_index shadow_exit_index D

gpssz = size(pos_vect,3);

yaw_rate = zeros(length(time),gpssz);
yaw = zeros(length(time),gpssz);


d2r = pi/180;       % convert from degrees to radians
r2d = 180/pi;       % convert from radians to degrees

% original code, now no longer used:
% compute satellite position unit vector and magnitude
%pos_vect_mag = norm(pos_vect);
%pos_vect_unit = pos_vect./pos_vect_mag;

% Compute the earth-sun unit vector (in ECI frame)
r_sun_vector = sun_vector(T);
mags = sqrt(dot(r_sun_vector,r_sun_vector,2));
mags = 1./mags;
magsMat = diag(mags);
r_sun_unit = r_sun_vector'*magsMat;

% Compute the orbit normal unit vector (in ECI frame)
orb_normal = cross(pos_vect, vel_vect);

% Not sure how to do this without a loop
tmp = zeros(length(T),gpssz);

orb_norm_unit = zeros(size(tx_ECI_pos_unit));
beta = zeros(length(time),gpssz);
h_cross_rsun  = zeros(size(tx_ECI_pos_unit));
h_cross_rsun_unit = zeros(size(tx_ECI_pos_unit));
for k = 1:1:gpssz 
    tmp(:,k) = sqrt(dot(orb_normal(:,:,k),orb_normal(:,:,k),1));
    
    %Build a 3xN matrix to divide against to avoid a second loop
    tmp2 = repmat(tmp(:,1),1,3);
    orb_norm_unit(:,:,k) = orb_normal(:,:,k)./tmp2';
    
    %Compute the beta angle in degrees
    % DEBUG!! the sin should be a cos, but not in the CO code, mabe a X-90
    % issue?
    beta(:,k) = asin(dot(orb_norm_unit(:,:,k),r_sun_unit,1))*180/pi;
    
    % determine the projection of the sun vector onto the orbit plane
    h_cross_rsun(:,:,k) = cross(orb_normal(:,:,k), r_sun_vector');
    h_cross_rsun_norm = sqrt(dot( h_cross_rsun(:,:,k), h_cross_rsun(:,:,k),1));
    h_cross_rsun_unit(:,:,k) = h_cross_rsun(:,:,k)./repmat(h_cross_rsun_norm,3,1);
    
    
end
clear tmp tmp2;

% This step can be completed out of the loop
rsun_in_plane = cross(h_cross_rsun, orb_normal);

rsun_in_plane_unit = zeros(size(tx_ECI_pos_unit));
alpha = zeros(length(time),gpssz);
for k = 1:1:gpssz 
    rsun_in_plane_mag = sqrt(dot(rsun_in_plane(:,:,k), rsun_in_plane(:,:,k),1));
    rsun_in_plane_unit(:,:,k) = rsun_in_plane(:,:,k)./repmat(rsun_in_plane_mag,3,1);

    % determine alpha (rad), the angle in the satellite's orbit plane between noon and
    % the satellite's position. "Noon" is defined as the direction of the sun vector
    alpha(:,k) = atan2(dot(tx_ECI_pos_unit(:,:,k), h_cross_rsun_unit(:,:,k)), dot(tx_ECI_pos_unit(:,:,k), rsun_in_plane_unit(:,:,k)));
end






% quadrant check: make sure 0 < alpha < 2*pi
alpha = ((alpha<0)*2*pi + (alpha<0).*alpha) + ((alpha >0).*alpha);

% compute the vector from the satellite to the sun
% Note: need to reverse the satellite position vector so it points from satellite to earth
r_sun_sat = rsun_in_plane + (-1)*pos_vect;

E = zeros(length(time),gpssz);
for k = 1:1:gpssz 
    r_sun_sat_mag = sqrt(dot(r_sun_sat(:,:,k), r_sun_sat(:,:,k),1));
    r_sun_sat_unit = r_sun_sat(:,:,k)./repmat(r_sun_sat_mag,3,1);
    
    % compute the earth-satellite-sun angle (in radians)
    E(:,k) = acos(dot((-1)*tx_ECI_pos_unit(:,:,k), r_sun_sat_unit));
end




% compute the actual yaw angle, B (in degrees), induced by the yaw bias, b, inserted in the
% satellite ACS
b = 0.5;      % degrees, set as of Nov. 1995
B = asin(0.0175*b*d2r./sin(E))*r2d;

% For E angles that are within +/- 5 degrees of noon or midnight, B becomes very large.
% An E angle of 5 degrees corresponds to a B of about 6 degrees. Therefore, since we
% don't know exactly what the ACS does in this range, we will set the maximum value of
% B at 6 degrees when E is within +/- 5 degrees of noon or midnight (see Bar-Sever).
B = min(abs(B),6).*sign(B);

% compute the nominal yaw angle
nom_yaw = get_nominal_yaw(beta*d2r, alpha, B);

% Compute nominal yaw rate
R = 0.13;     % maximal yaw rate, deg/sec,
% ftp://sideshow.jpl.nasa.gov/pub/GPS_yaw_attitude/nominal_yaw_rates
RR = 0.0017;  % maximal yaw rate rate, deg/sec^2 (Bar-Sever)

% "mu" is defined by Bar-Sever to be the orbit angle, measured from orbit
% midnight in the direction of motion. Thus, it is alpha + 180 deg. This
% angle is assumed to vary little over time, and its rate can be replaced
% by a constant.
mu_dot = 0.0083;     % deg/sec
% Compute B_dot (deg/sec), the rate of change of B
B_dot = -0.0175*b*cos(E).*cos(beta*d2r).*sin(alpha+pi)*mu_dot./(cos(B*d2r).*(sin(E)).^3);
% Compute the nominal yaw rate (deg/sec)
nom_yaw_rate = (mu_dot*tan(beta*d2r).*cos(alpha+pi)./(sin(alpha+pi).^2+tan(beta*d2r).^2)) + B_dot;


dt = 1;
ISNOON = zeros(gpssz,1);
ISSHADOW = zeros(gpssz,1);
shadow_entry_index = zeros(gpssz,1);
shadow_exit_index = zeros(gpssz,1);
t_i = zeros(gpssz,1);
t_spin = zeros(gpssz,1);
D = zeros(gpssz,1);

% Determine if satellite is in the "noon maneuver regime"; this occurs in
% the vicinity of orbit noon when the nominal yaw rate would be higher than
% the maximal yaw rate allowed by the satellite.
for time_tag = 1:1:length(T) 
    for PRN_tag = 1:1:gpssz 
        if time_tag == 1
            ISNOON(PRN_tag) = 0;
            ISSHADOW(PRN_tag) = 0;
            yaw(time_tag, PRN_tag) = nom_yaw(time_tag, PRN_tag);
            yaw_rate(time_tag, PRN_tag) = nom_yaw_rate(time_tag, PRN_tag);
        end
        
        
        
        if beta > -5 & beta < 5 & time_tag ~= 1
            if ISNOON(PRN_tag) == 0
                % satellite is NOT in the noon maneuver regime
                yaw(time_tag, PRN_tag) = nom_yaw(time_tag, PRN_tag);
                yaw_rate(time_tag, PRN_tag) = nom_yaw_rate(time_tag, PRN_tag);
                
                if abs(nom_yaw_rate(time_tag, PRN_tag)) > R
                    % satellite is beginning to enter the noon maneuver regime
                    yaw_rate(time_tag, PRN_tag) = R * sign(nom_yaw(time_tag, PRN_tag)-nom_yaw(time_tag-1, PRN_tag));
                    ISNOON(PRN_tag) = 1;
                end
            else
                % satellite is already in noon maneuver regime
                yaw(time_tag, PRN_tag) = yaw(time_tag-1, PRN_tag) + yaw_rate(time_tag-1, PRN_tag) * dt;
                yaw_rate(time_tag, PRN_tag) = R * sign(nom_yaw(time_tag, PRN_tag)-nom_yaw(time_tag-1, PRN_tag));
                
                if yaw(time_tag, PRN_tag) > nom_yaw(time_tag, PRN_tag) && (nom_yaw(time_tag-1, PRN_tag)-nom_yaw(time_tag, PRN_tag))<=0
                    % satellite is on the "upslope" and has overshot the nominal yaw angle;
                    % satellite returns to nominal regime
                    ISNOON(PRN_tag) = 0;
                elseif yaw(time_tag, PRN_tag) < nom_yaw(time_tag, PRN_tag) && (nom_yaw(time_tag-1, PRN_tag)-nom_yaw(time_tag, PRN_tag))>=0
                    % satellite is on the "downslope" and has undershot the nominal yaw angle;
                    % satellite returns to nominal regime
                    ISNOON(PRN_tag) = 0;
                end
            end   % end ISNOON loop
            
        else
            % beta not in "noon" range; satellite is in nominal regime
            yaw(time_tag, PRN_tag) = nom_yaw(time_tag, PRN_tag);
            yaw_rate(time_tag, PRN_tag) = nom_yaw_rate(time_tag, PRN_tag);
        end  % end beta loop
        
        switch ISSHADOW(PRN_tag)
            case 0
                if ISNOON(PRN_tag) == 0
                    % satellite is in nominal range
                    yaw(time_tag, PRN_tag) = nom_yaw(time_tag, PRN_tag);
                    yaw_rate(time_tag, PRN_tag) = nom_yaw_rate(time_tag, PRN_tag);
                    
                    % check to see if satellite is entering shadow
                    shadow_flag = compute_shadow(r_sun_vector(time_tag,:), pos_vect(:,time_tag,PRN_tag));
                    if shadow_flag == 1
                        ISSHADOW(PRN_tag) = 1;
                        shadow_entry_index(PRN_tag) = time_tag;
                        t_i(PRN_tag) = time(shadow_entry_index(PRN_tag));
                        
                        % determine yaw rate
                        if sign(beta) ~= sign(b)
                            yaw_rate(time_tag, PRN_tag) = R*-1;
                        else
                            yaw_rate(time_tag, PRN_tag) = R;
                        end
                        % determine spin-up/down time
                        t_spin(PRN_tag) = (R - yaw_rate(shadow_entry_index(PRN_tag),PRN_tag))/RR;
                    end
                end   % end ISNOON check
                
            case 1
                % satellite is in shadow-crossing regime
                if time(time_tag) <= (t_i(PRN_tag) + t_spin(PRN_tag))
                    yaw(time_tag, PRN_tag) = yaw(PRN_tag, shadow_entry_index(PRN_tag)) + yaw_rate(PRN_tag, shadow_entry_index(PRN_tag))*...
                        (time(time_tag) - t_i(PRN_tag)) + 0.5*RR*(time(time_tag) - t_i(PRN_tag))^2;
                    yaw(time_tag, PRN_tag) = mod(yaw(time_tag, PRN_tag),360);
                else
                    yaw(time_tag, PRN_tag) = yaw(shadow_entry_index(PRN_tag),PRN_tag) + yaw_rate(shadow_entry_index(PRN_tag),PRN_tag )*...
                        t_spin(PRN_tag) + 0.5*RR*t_spin(PRN_tag)^2 + R*(time(time_tag) - t_i(PRN_tag) - t_spin(PRN_tag));
                    yaw(time_tag, PRN_tag) = mod(yaw(time_tag, PRN_tag),360);
                end
                
                % check to see if satellite is leaving shadow regime
                shadow_flag = compute_shadow(r_sun_vector(time_tag,:), pos_vect(:,time_tag,PRN_tag));
                if shadow_flag == 0
                    % satellite is entering post-shadow maneuver
                    ISSHADOW(PRN_tag) = 2;
                    shadow_exit_index(PRN_tag) = time_tag;
                    D(PRN_tag) = nom_yaw(shadow_exit_index(PRN_tag),PRN_tag) - yaw(shadow_exit_index(PRN_tag),PRN_tag) -...
                        round((nom_yaw(shadow_exit_index(PRN_tag),PRN_tag) - yaw(shadow_exit_index(PRN_tag),PRN_tag))/360)*360;
                end
            case 2
                % satellite is in post-shadow maneuver
                yaw_rate(PRN_tag,time_tag) = R * sign(D(PRN_tag));
                t_1 = (yaw_rate(PRN_tag,time_tag) - R)/(RR*sign(D(PRN_tag)));
                % compute the yaw angle
                if time(time_tag) < (time(shadow_exit_index(PRN_tag)) + t_1)
                    yaw(time_tag, PRN_tag) = yaw(shadow_exit_index(PRN_tag),PRN_tag) + R*(time(time_tag)-time(shadow_exit_index(PRN_tag)))+...
                        0.5*RR*sign(D(PRN_tag))*(time(time_tag)-time(shadow_exit_index(PRN_tag)))^2;
                    yaw(time_tag, PRN_tag) = mod(yaw(time_tag, PRN_tag),360);
                else
                    yaw(time_tag, PRN_tag) = yaw(shadow_exit_index(PRN_tag),PRN_tag) + R*t_1+...
                        0.5*RR*sign(D(PRN_tag))*t_1^2 + R*sign(D(PRN_tag))*(time(time_tag)-time(shadow_exit_index(PRN_tag))-t_1);
                    yaw(time_tag, PRN_tag) = mod(yaw(time_tag, PRN_tag),360);
                end
                
                % check to see if satellite is leaving post-shadow maneuver
                d = yaw(time_tag, PRN_tag) - nom_yaw(time_tag, PRN_tag);
                % make sure -180 < d < +180 degrees
                if d > 180
                    d = d - 360;
                elseif d < -180
                    d = d + 360;
                end
                if sign(d) == sign(yaw_rate(time_tag, PRN_tag))
                    % satellite has completed post-shadow maneuver; return to nominal
                    ISSHADOW(PRN_tag) = 0;
                end
        end
    end   % end switch/case
end
end


%*****************************************************************
% Title: sun_vector
% Author: Lisa Reeh
% Date written: 14 Feb. 2002
% Date last modified: 27 May 2003
% Copyright 2002 Univesity of Colorado, Boulder
%
% Purpose: computes the vector from the center of the earth
% to the center of the sun
%
% Input:
%       T: number of elapsed Julian centuries (J2000 epoch)
%
% Output:
%       r_sun_vector: ECI sun vector (meters)
%
% Reference: Vallado, p 183 (1st ed.), OR p.263 (2nd ed.)
%*****************************************************************

function r_sun_vector = sun_vector(T)

d2r = pi/180;
AU = JATConstant('AU');  % in meters

% determine the mean longitude of the sun in degrees
lambda = mod( (280.4606184 + 36000.77005361 * T), 360 );

% determine the mean anomaly of the sun in degrees
M = mod( (357.5277233 + 35999.05034 * T), 360 );

% Calculate the ecliptic longitude in degrees
lambda_ecl = mod(lambda + 1.914666471 * sin(M*d2r) + 0.019994643 * sin(2*M*d2r), 360);

% compute the magnitude of the sun vector
r_sun_mag = 1.000140612 - 0.016708617*cos(M*d2r) - 0.000139589*cos(2*M*d2r);

% approximate the obliquity of the ecliptic
ob_ecl = mod(23.439291 - 0.0130042 * T, 360);

% Compute sun vector in ECI frame
rm = repmat(r_sun_mag,1,3);  %use this to remove the requirement of a loop

% Lots of transposes to keep matrix dimensions consistant
r_sun_vector = AU*rm.*([cos(lambda_ecl.*d2r)'; cos(ob_ecl.*d2r)'.*sin(lambda_ecl.*d2r)';...
    sin(ob_ecl.*d2r)'.*sin(lambda_ecl.*d2r)'])';


end


%***********************************************************************************
% Title: get_nominal_yaw
% Author: Lisa Reeh
% Date written: 8 March 2002
% Date last modified: 8 March 2002
% Copyright 2002 University of Colorado, Boulder
%
% Purpose: computes nominal yaw angle of satellite in degrees
% Inputs: beta angle - angle between sun vector and orbit plane, measured in radians
%         alpha angle - angle in the satellite's orbit plane between noon and
%             the satellite's position, measured in radians
%         B - actual yaw angle induced by the yaw bias inserted in the satellite
%             ACS, measured in degrees
%
% Output: yaw angle, in degrees
%
% Reference: Bar-Sever, Yoaz E. "A new model for GPS yaw attitude," Journal
% of Geodesy (1996) 70:714-723
%***********************************************************************************

function nom_yaw = get_nominal_yaw(beta, alpha, B)

r2d = 180/pi;    % convert radians to degrees

% compute the nominal yaw angle

%if beta > 0
tmp = atan((tan(beta)./sin(alpha)))*r2d;
nom_yaw = tmp.*(beta>0);

%if alpha is between pi and 2*pi, add 180 degrees + B
V = (alpha > pi & alpha <= 2*pi & beta > 0);
nom_yaw = nom_yaw + V*180 + B.*V;

    
%if beta < 0
clear tmp;
tmp = mod( (atan( tan(beta)./sin(alpha) ))*r2d, 360);
nom_yaw2 = tmp.*(beta<0);    

%if alpha > than pi, add 180 + B;
V = (alpha > pi &  beta < 0);
nom_yaw2 = nom_yaw2 + V*180 + B.*V;

%Now sum the two cases
nom_yaw = nom_yaw + nom_yaw2;
end


%****************************************************************************
% Title: compute_shadow
% Author: Lisa Reeh
% Date written: Jan. 2003
% Date last modified: 27 May 2003
% Copyright 2002 University of Colorado, Boulder
%
% Purpose: Determines if the current GPS satellite is in the earth's
%       penumbra region.
%
% Input: 
%       r_sun_vector: vector from earth to sun (m, ECI)
%       sat_pos_vect: position vector of the satellite (m, ECI)
%
% Output:
%       shadow: (1) if satellite is in earth penumbra
%               (0) if satellite is not in earth penumbra
%
% Reference: Vallado, 2nd ed., Algorithm 34, with slight modifications
%****************************************************************************

function shadow = compute_shadow(r_sun_vector, sat_pos_vect)

R_sun = 695990000;      % sun's radius, m (Conway & Prussing)    
R_earth = 6378136.3;    % earth's radius, m (Vallado, 2nd ed., JGM-2)
r_sun_mag = norm(r_sun_vector);
sat_pos_mag = norm(sat_pos_vect);
angle_penumbra = atan2((R_sun + R_earth),r_sun_mag);

if dot(r_sun_vector, sat_pos_vect) < 0
    % determine the angle between r_sun_vector and sat_pos_vect
    angle = acos(dot(r_sun_vector, sat_pos_vect)/(r_sun_mag * sat_pos_mag));
    
    sat_horiz = sat_pos_mag * cos(angle);
    sat_vert = sat_pos_mag * sin(angle);
    pen_vert = R_earth + tan(angle_penumbra)*sat_horiz;
    
    if sat_vert <= pen_vert
        shadow = 1;
    else
        shadow = 0;
    end
else
    shadow = 0;
end
end


%********************************************************************************
% Title: eci2ant
% 
% Purpose: convert the LOS ECI vector to local coordinates with
% respect to the satellite antenna.
% 
% Input:    yaw - ideal yaw angle, degrees
%           tx_sat_pos - ECI satellite position, in meters
%           tx_sat_vel - ECI satellite velocity, in meters/second
%           LOS_unit - Line of sight unit vector, ECI frame, in meters
% Output:   LOS_ant - Line of sight unit vector, local antenna frame, in meters
%
% Original function information:
% Author: Lisa Reeh
% Date written: 16 May 2002
% Date last modified: 16 May 2002
% Copyright 2002 University of Colorado, Boulder
%********************************************************************************

function [LOS_ant] = eci2ant(yaw, tx_sat_pos, tx_sat_vel, LOS_unit)

LOS_ant = zeros(size(tx_sat_pos));

%Need to unwrap this loop, but not sure how

for j = 1:1:size(tx_sat_pos,2)
    
    for k = 1:1:size(yaw,2)
        yaw_rad = yaw(j,k)*pi/180;
        
        ROT1 = [1, 0, 0; 0, cos(yaw_rad), sin(yaw_rad); 0, -sin(yaw_rad), cos(yaw_rad)];
        
        R_hat = tx_sat_pos(1:3,j,k)./norm(tx_sat_pos(1:3,j,k));
        C_hat = cross(tx_sat_pos(1:3,j,k), tx_sat_vel(1:3,j,k))./norm(cross(tx_sat_pos(1:3,j,k), tx_sat_vel(1:3,j,k)));
        I_hat = cross(C_hat, R_hat);
        
        RIC_matrix = [R_hat, I_hat, C_hat];
        
        LOS_ant(1:3,j,k) = ROT1 * RIC_matrix' * LOS_unit(1:3,j,k);
    end
end
end


