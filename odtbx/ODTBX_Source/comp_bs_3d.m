% COMP_BS_3D    Computes unit vectors describing antenna boresite 
%
%  Computes unit vectors corresponding to antenna boresites corresponding
%  to the attitude profile specified in the _gpsdef.m file
%  Called by GPSSIM1, STATICVIS1, GPSSIM_AO40
%
%  [boresite] = comp_bs(num_ant,time,r_sat,v_sat,pointing_ref,frame)
%
%  Inputs:    num_ant       Number of antennas, 0-4
%             time          Time vector in total seconds since start GPS time
%             r_sat         [3,nn] vector containing host vehicle position vectors 
%                           (from the variable sat_pos).  Position vectors given in 
%                           the refrence frame specified by 'frame'
%             v_sat         [3,nn] vector containing host vehicle velocity vectors 
%                           (from the variable sat_vel in "frame").  
%             pointing_ref  (1) zenith pointing or (-1) nadir pointing 
%                           (2) parallel or (-2) antiparallel to the Earth-Sun vector
%                           (3) perpendicular to the ecliptic plane, +Z or -Z (-3)
%                           (4) fixed with respect to apogee zenith, or (-4) nadir vector
%                               apogee selected as the point of highest altitude (ephemeris 
%                               must include apogee)
%                           (5) body fore and (-5) aft directions relative to geocentric
%                               LVLH 
%                           (6) body port and (-6) starboard directions relative to 
%                               geocentric LVLH 
%             frame         Ref. Frame: 0-ECI, 1-ECEF
%
%  Outputs:   boresite      [3,nn,num_ant] matrix containing bs att.
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
%   Mike Moreau         12/26/2001         Original
%   Mike Moreau         11/19/2002         Updated
%   Mike Moreau         12/09/2002         made boresite a 3D array
%   Kevin Berry         09/21/2007         Modified for ODTBX
%   Keith Speckman      03/13/2008         Minor comment change
%   Sun Hur-Diaz        02/07/2011         Fixed bugs in pointing opts 5 and 6

function [boresite] = comp_bs_3d(num_ant,time,r_sat,v_sat,pointing_ref,frame)

%  Check for valid number of antennas
if ~(num_ant==0 || num_ant==1 || num_ant==2 || num_ant==3 || num_ant==4)
   error('Number of antennas must be an integer between 0-4')
end

% Convert GPS Seconds to other formats
gpsDays = time'/86400;
% gpsWeeks(:,1) = floor(time'./(86400*7));
% gpsWeeks(:,2) = time' - gpsWeeks(:,1).*86400*7;

%  Initialize boresite matrix
boresite = zeros(3,length(time),num_ant);

if num_ant == 0
   %  Host vehicle position vectors (r_sat) passed as function input (sat_pos)
   %  Compute magnitude of position vector, the radial distance to host vehicle
   mag_r_sat = sqrt(sum(r_sat.^2));                         % [1,nn]
   %  Compute unit vectors corresponding to r_sat
   unit_r_sat = r_sat./[mag_r_sat;mag_r_sat;mag_r_sat;];    % [3,nn]
   
   %  Set boresites parallel to r_sat direction (zenith direction)
   ant_bs = unit_r_sat;                                   % [3,nn]
   
   %  Boresite is based on r_sat vector, so boresite is automatically in the
   %  correct reference frame
   boresite = ant_bs;   
end

for ii = 1:num_ant
   if abs(pointing_ref(ii)) == 1   % Earth pointing s/c
      %  Host vehicle position vectors (r_sat) passed as function input (sat_pos)
      %  Compute magnitude of position vector, the radial distance to host vehicle
      mag_r_sat = sqrt(sum(r_sat.^2));                      % [1,nn]
      %  Compute unit vectors corresponding to r_sat
      unit_r_sat = r_sat./[mag_r_sat;mag_r_sat;mag_r_sat;]; % [3,nn]
      
      %  Set boresites parallel to r_sat direction if pointing_ref = 1
      %  Set boresites anti-parallel to r_sat direction if pointing_ref = -1
      ant_bs = unit_r_sat.*(pointing_ref(ii)/abs(pointing_ref(ii)));
      
      %  Boresite is based on r_sat vector, so boresite is automatically in the
      %  correct reference frame
   end
   
   if abs(pointing_ref(ii)) == 2   % Sun pointing s/c
      % Compute sun vector
      [sunv] = sun_vec(time./(24*3600));   % [3,nn]
      
      %  Set boresites parallel to Earth-Sun vector if pointing_ref = 2
      %  Set boresites anti-parallel to Earth-Sun vector if pointing_ref = -2
      bs = sunv.*(pointing_ref(ii)/abs(pointing_ref(ii)));
      
      %  Boresite is referenced to ECI, if frame = 1, must rotate to ECEF
      if frame == 1
         ant_bs = zeros(size(bs));
         utc = convertTime('UTC','GPS', gpsDays);
         xform = jatDCM('eci2ecef', utc, 14);
         for jj = 1:size(time,2)
             ant_bs(:,jj) = xform(:,:,jj)*bs(:,jj);
         end
      else
         ant_bs = bs;
      end
   end
   
   if abs(pointing_ref(ii)) == 3   % s/c alligned perpendicular to ecliptic

      t_mat = convertTime('UTC','GPS', gpsDays);
      if(frame == 1)
         xform = jatDCM('eci2ecef', t_mat, 14);
      end

      % Compute the list of times as Julian Dates.
      t_JD = matlabTime2JD(t_mat);

      bs = zeros(3, size(t_JD, 2));

      ant_bs = zeros(3, size(t_JD, 2));

      for jj=1:size(t_JD, 2)

          % Calculate the time in Julian Centuries since the ephemeris epoch of 2000.
          [T_JC] = calc_JulianCenturies2000(t_JD(jj));

          % Calculate the obliquity of the ecliptic.
          [ob] = calc_obliquity(T_JC);

          %  Create the unit vector normal to the ecliptic plane, written in ECI
          eclipt_norm = [1 0 0; 0 cos(ob) sin(ob); 0 -sin(ob) cos(ob)]'*[0;0;1];

          %  Set boresites parallel to eclipt_norm if pointing_ref = 3
          %  Set boresites anti-parallel to eclipt_norm if pointing_ref = -3
          bs(:,jj) = sign(pointing_ref(ii))*eclipt_norm;

          %  Boresite is referenced to ECI, if frame = 1, must rotate to ECEF
          if(frame == 1)
              ant_bs(:,jj) = xform(:,:,jj)*bs(:,jj);
          else
              ant_bs(:,jj) = bs(:,jj);
          end

      end

   end
   
   if abs(pointing_ref(ii)) == 4   % inertially fixed wrt apogee nadir vector
     
       if frame  % working in ecef frame
           %  r_sat is in ECEF, so compute bs in ECI then rotate to ECEF
           %  Compute r_sat in ECI frame
           utc = convertTime('UTC', 'GPS', gpsDays);
           r_sat_eci = zeros(size(r_sat));
           xform = jatDCM('ecef2eci', utc, 14);
           for jj=1:size(utc, 2)
               r_sat_eci(:,jj) = xform(:,:,jj)*r_sat(:,jj);
           end
           
           %  Compute magnitude of position vector, the radial distance to host vehicle
           mag_r_sat = sqrt(sum(r_sat_eci.^2));               % [1,nn]
           
           % find the apogee point (highest altitude) of the current data set
           ind = find(mag_r_sat == max(mag_r_sat));
           % If the maximum is either the first or last point, data may not contain apogee point
           if (ind == 1) || (ind == length(time))
               warning('Warning:OrbitLength','Apogee at end point, make sure ephemeris covers >= one full orbit!')
           end
           
           %  Compute unit vector corresponding to r_sat at apogee point
           unit_r_sat = r_sat_eci(:,ind)./mag_r_sat(ind);  % [3,1]
           
           %  Set boresites parallel to r_sat direction if pointing_ref = 1
           %  Set boresites anti-parallel to r_sat direction if pointing_ref = -1
           bs_eci = ( unit_r_sat*ones(1,length(time)) ) .* (pointing_ref(ii)/abs(pointing_ref(ii)));
           
           % convert boresite unit vectors back to ECEF
           ant_bs = zeros(size(bs_eci));
           xform = jatDCM('eci2ecef', utc, 14);
           for jj=1:size(utc, 2)
               ant_bs(:,jj) = xform(:,:,jj)*bs_eci(:,jj);
           end
           
       else % working in eci frame
           %  r_sat is in ECI, so boresite will be returned as ECI
           %  Compute magnitude of position vector, the radial distance to host vehicle
           mag_r_sat = sqrt(sum(r_sat.^2));                      % [1,nn]
           
           % find the apogee point (highest altitude) of the current data set
           ind = find(mag_r_sat == max(mag_r_sat));
           % If the maximum is either the first or last point, data may not contain apogee point
           if (ind == 1) || (ind == length(time))
               warning('Warning:OrbitLength','Apogee at end point, make sure ephemeris covers >= one full orbit!')
           end

           %  Compute unit vector corresponding to r_sat at apogee point
           unit_r_sat = r_sat(:,ind)./mag_r_sat(ind);  % [3,1]
           
           %  Set boresites parallel to r_sat direction if pointing_ref = 1
           %  Set boresites anti-parallel to r_sat direction if pointing_ref = -1
           ant_bs = ( unit_r_sat*ones(1,length(time)) ) .* (pointing_ref(ii)/abs(pointing_ref(ii)));
       end
   end
   
   if abs(pointing_ref(ii)) == 5 || abs(pointing_ref(ii)) == 6
       % 5 = antennas in fore and aft directions wrt nadir pointing body frame in LVLH
       % 6 = antennas in port and starboard directions wrt nadir pointing body frame in LVLH
       
       if frame  % Vectors given in ECEF; 
           % Need to convert the relative velocity to total velocity
           
           w = repmat([0; 0; JATConstant('wEarth')],1,size(r_sat,2));
           v = v_sat + cross(w,r_sat,1); % Add the omega cross term
       
       else  % vectors given in eci, velocity already total
           v = v_sat;
       end
       
       r = r_sat;
       
       R = dcm('ric', r,v);
       
       % Select appropriate row from the transformation matrix
       % Row 2 corresponds to in-track (fore) and row 3 corresponds to 
       % cross-track (port); therefore -3+abs(pointing_ref(ii)) selects the
       % right row.
       
       ant_bs = sign(pointing_ref(ii)) * squeeze(R(-3+abs(pointing_ref(ii)),:,:));
       
       clear r v
   end  % if '5'

   %  Save boresites in one composite boresite matrix [3*num_ant,nn]
   boresite(:,:,ii)=ant_bs;
end
