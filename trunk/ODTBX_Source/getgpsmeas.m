function out = getgpsmeas(t,x,options,qatt,params)
% GETGPSMEAS  Computes physical parameters required for GPS based measurements 
%
%   out = GETGPSMEAS(t,x,options) computes physical parameters such as LOS 
% angles, vectors, range, and range rate required for GPS measurements.  
% These are based on the information in OPTIONS. See the OD Toolbox 
% function of the same name for details of each measurement type.  This
% function is called by the functions GPSMEAS and gps_phys_params.  
%
% Note that this function can be called in two different ways based on the
% params.GPS_SIZE and params.PRN.  If called to compute all satellites in
% the almanac set parmas.GPS_SIZE to 32 and leave params.PRN empty.  The
% outputs are sized to 32 and indexed via the GPS PRN value.  However, if
% only one satellite computation is desired, set params.PRN to the desired
% GPS PRN and set parms.GPS_SIZE to 1.
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are used in this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   epoch               datenum                Time associated with start
%                                              of simulation, UTC
%   YumaFile            filename               GPS almanac file
%   Rotation2ECI        @functionname          A function that allows for transformation of the state
%                                              vector, x, from an arbitrary coodinate system to ECI at a
%                                              given time. Returns a 6x6 matrix.  Input must be a pointer to a
%                                              rotation function that uses time as an input. See the
%                                              embedded IdentRot() function for details.
%   AntennaPointing     1 x num_ant            Specify attitude profile for each antenna
%                                              (1) zenith pointing or (-1) nadir pointing wrt geocentric LVLH
%                                              WARNING: "1" and "-1"  values work, these other
%                                              types are not yet implemented:
%                                              (2) parallel or (-2) antiparallel to the Earth-Sun vector
%                                              (3) ecliptic north or south (-3)
%                                              (4) fixed with respect to apogee zenith, or (-4) nadir vector apogee 
%                                                selected as the point of highest altitude (ephemeris must include apogee)
%                                              (5) body fore and (-5) aft directions relative to geocentric LVLH
%                                              (6) body port and (-6) starboard directions relative to geocentric LVLH
%                                              NOTE: The AntennaPointing
%                                              parameter is used only if
%                                              the receiver antenna pattern
%                                              specified by AntennaPattern
%                                              is 1-D.  If the pattern is
%                                              2-D, then the attitude of the
%                                              antennas is determined from
%                                              the spacecraft body attitude
%                                              quaternion specified in the
%                                              input variable x as well as
%                                              the AntennaOrientation
%                                              parameter in this options
%                                              structure.
%   PrecnNutnExpire     days                   For how long, in days, a
%                                              computed value for
%                                              precession and nutation of
%                                              the Earth can be used before
%                                              new values need to be
%                                              recomputed.  Reusing values
%                                              improves performance of
%                                              ECI-to-ECEF conversions
%                                              which are done inside the
%                                              GPS computations.  The
%                                              default value is 1 day.
%   AntennaOrientation  3 x 3 x num_ant        unit orthogonal matrices
%                                              representing the rotation of
%                                              each antenna with respect to
%                                              the body.  The Z-axis of the 
%                                              antenna frame is the 
%                                              boresite.  This option
%                                              parameter is used only if
%                                              the corresponding antenna
%                                              pattern is 2-D.  If it's 1-D
%                                              then the antenna boresite is
%                                              computed based on the
%                                              antenna pointing specified
%                                              in AntennaPointing
%                                              parameter. Default is
%                                              identity matrix.
%
%
%   INPUTS:
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch in options)
%      x            (6xN)   spacecraft state [position;velocity] (km)
%      options      (1x1)   data structure (see above description)
%      qatt         (4xN)   spacecraft body quaternion history where the
%                           quaternion specifies the user spacecraft
%                           attitude with respect to the same coordinate
%                           frame that the position and velocity are given
%                           in.  The quaternion is specified as
%                           [sin(theta/2)*e_vec;cos(theta/2)] where e_vec is the
%                           unit vector representing the axis of rotation
%                           and theta is the angle of rotation about that
%                           axis.  The quaternion information is used IF a 
%                           2-D antenna is specified for the user satellite.  
%                           If the quaternion is not supplied and the 2-D 
%                           antenna is specified, then it will assume that 
%                           spacecraft attitude is aligned with the 
%                           coordinate frame that the position and velocity 
%                           are given in, i.e., the default is [0 0 0 1]'.
%                           Additionally, the orientation matrix
%                           of each antenna wrt to the spacecraft body can 
%                           be specified through the options strucutre.  If
%                           the antennae orientation is not specified, the
%                           default is identity, i.e., aligned with the body
%                           axes.
%     params        (1x1)   Input structure containing the following
%                           fields:
%                              params.num_ant           (1 x 1)       
%                                  Number of receiver antennas
%                              params.xmit_pattern_dim  (1 x 1)       
%                                  1 or 2 for dimension of the transmit antenna pattern
%                              params.rec_pattern_dim   (num_ant x 1) 
%                                  1 or 2 for dimension of the receiver antenna patterns
%                              params.GPS_SIZE          (1 x 1)       
%                                  Number of GPS satellites to calculate,
%                                  allowable values: 1 or 32
%                              params.doH               (1 x 1)       
%                                  Flag indicating whether the parameters
%                                  required for the computation of the 
%                                  measurement partials are required
%                              params.PRN               (1 x 1)       
%                                  OPTIONAL, forces calculation for only one
%                                  GPS PRN, and designates the index of the
%                                  desired PRN to be calculated, note
%                                  GPS_SIZE must be one.
% 
%   OUTPUTS
%      out           (1x1)  Output structure containing the following
%                           fields:
%                              out.epoch      (1 x 1) 
%                                  epoch in Matlab datenum
%                              out.TX_az      (N x GPS_SIZE) 
%                                  transmitter azimuth [deg]
%                              out.TX_el      (N x GPS_SIZE) 
%                                  transmitter elevation [deg]
%                              out.RX_az      (N x GPS_SIZExnum_ant) 
%                                  receiver azimuth [deg]
%                              out.RX_el      (N x GPS_SIZExnum_ant) 
%                                  receiver elevation [deg]
%                              out.range      (GPS_SIZE x N) 
%                                  range [km]
%                              out.rrate      (GPS_SIZE x N) 
%                                  range rate [km/sec]
%                              out.GPS_yaw    (N x GPS_SIZE) 
%                                  GPS SV yaw angle
%                              out.rgps_mag   (N x GPS_SIZE) 
%                                  magnitude of the GPS SV position [km]
%                              out.health     (N x GPS_SIZE) 
%                                  health flag of each GPS SV
%                              out.prn        (1 x H) 
%                                  PRN ID of each healthy GPS SV in this out struct
%
%                              If the params.doH == 1, then the following 
%                              fields are also set:
%                                  out.eciRotation   (9 x 9 x N) 
%                                     rotation matrix from ECI to ECEF
%                                  out.los_3d        (3 x N x GPS_SIZE) 
%                                     LOS 3D vector
%                                  out.gps_vel_tot   (3 x N x GPS_SIZE) 
%                                     total velocity of GPS SVs in ECEF
%                                     frame [km/s]
%                                  out.sat_vel_tot   (3 x N) 
%                                     total velocity of the satellite in
%                                     ECEF frame [km/s]
%                                  out.Rotation2ECI  (1 x 1) 
%                                     name of function that computes 
%                                     rotation from input frame to ECI
%
%
%   keyword: measurement, physical parameters, GPS
%   See also odtbxOptions, gpsmeas, gps_phys_params
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
%   Sun Hur-Diaz        03/02/2011      Original refactoring from the
%                                       original gpsmeas.m

%% Persistent variables for file data caching
persistent YumaCache;

%% If required, initialize any persistent variables
if isempty(YumaCache)
    YumaCache = dataCache('create');
end

%% Set some constants for this call
r2d = 180/pi;
nn = length(t);
loop = max([1,params.num_ant]);
GPS_SIZE = params.GPS_SIZE;

%% check arguments
if ~isfield(params,'GPS_SIZE') || isempty(params.GPS_SIZE)
    error('getgpsmeas was passed missing input: params.GPS_SIZE');
end
if params.GPS_SIZE ~= 1 && params.GPS_SIZE ~= 32
    error('getgpsmeas was passed invalid input: params.GPS_SIZE = %d (must be 1 or 32)',params.GPS_SIZE);
end

if isfield(params,'PRN') && ~isempty(params.PRN)
    if (GPS_SIZE ~= 1)
        error('getgpsmeas was passed conflicting inputs: params.GPS_SIZE=%d vs. params.PRN=%d.',GPS_SIZE,params.PRN);
    end
else
    % no usable params.PRN, always set for checks below
    params.PRN = [];
end

%% Get values from options
epoch           = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') );
precnNutnExp    = getOdtbxOptions(options, 'PrecnNutnExpire', 1);
Rotation2ECI    = getOdtbxOptions(options, 'Rotation2ECI', @IdentRot );
                %For coordinate frames other than ECI. This value allows
                %the user to have x in any coordinate frame as long as the
                %rotation to ECI goes with it. Input must be a pointer to a
                %rotation function with time as the input.
YumaFile        = getOdtbxOptions(options, 'YumaFile', 'Yuma1356.txt' );
pointing_ref    = getOdtbxOptions(options, 'AntennaPointing', [-1,1] );
                % Specify attitude profile for each antenna
                %   (1) zenith pointing or (-1) nadir pointing wrt geocentric LVLH
                %   (2) parallel or (-2) antiparallel to the Earth-Sun vector
                %   (3) ecliptic north or south (-3)
                %   (4) fixed with respect to apogee zenith, or (-4) nadir
                %   vector apogee 
                %       selected as the point of highest altitude (ephemeris must include apogee)
                %   (5) body fore and (-5) aft directions relative to geocentric LVLH
                %   (6) body port and (-6) starboard directions relative to geocentric LVLH
ant_body        = getOdtbxOptions(options, 'AntennaOrientation',reshape(repmat(eye(3),1,params.num_ant),3,3,params.num_ant));

%% Create the spacecraft position and velocity
% For GPS calculations, it is best to Rotate to ECEF
init2fixed  = jatDCM('eci2ecef', (t/86400)+epoch, precnNutnExp);
w           = [0;0;JATConstant('wEarth')];
eciRotation = rotransf(repmat(w,[1,nn]),init2fixed); % Removes omega cross term
sat_pos = zeros(3,nn);
sat_vel = zeros(3,nn);
sat_vel_tot = sat_vel;          % Total velocity in ECEF frame
for n=1:nn
    R2ECI = Rotation2ECI((t(n)/86400)+epoch);
    xf             = eciRotation(1:6,1:6,n) * R2ECI * x(1:6,n); %km
    sat_pos(1:3,n) = xf(1:3);
    sat_vel(1:3,n) = xf(4:6);
    sat_vel_tot(:,n) = init2fixed(:,:,n) * R2ECI(4:6,:) * x(1:6,n);
end


%% Get the positions of the GPS satellites
%  Convert simulation start time to GPS time in seconds
time = 86400 * convertTime('GPS','UTC', epoch+t/86400);

% Read the input file
gps_alm = dataCache('get',YumaCache,YumaFile);
if isempty(gps_alm)
    % Not already in the cache.  Read it now and save it in the cache.
    gps_alm = read_yuma(YumaFile,epoch,0);
    YumaCache = dataCache('add', YumaCache, YumaFile, gps_alm);
end

% Find the list of PRNs with valid almanac health flag
healthy = gps_alm(:,2)==0;
healthy_ind = find(healthy);

% check against params.PRN request
if ~isempty(params.PRN) 
    if any(healthy_ind == params.PRN)
        healthy_ind = params.PRN; %override and select only this PRN
    else
        healthy_ind = []; % bad selection vs almanac health
    end
end

prns = gps_alm(healthy_ind,1);

% No healthy SVs in the healthy_ind, early exit w/ no computation.
% This is due to the almanac file, the GPS_SIZE, or by selecting a PRN
% that isn't healthy.
if isempty(healthy_ind)
    out.epoch = 0;
    out.TX_az = zeros(nn,GPS_SIZE);
    out.TX_el = zeros(nn,GPS_SIZE);
    out.RX_az = zeros(nn,GPS_SIZE,loop);
    out.RX_el = zeros(nn,GPS_SIZE,loop);
    out.range = zeros(GPS_SIZE,nn);
    out.rrate = zeros(GPS_SIZE,nn);
    out.GPS_yaw = zeros(nn,GPS_SIZE);
    out.rgps_mag = zeros(nn,GPS_SIZE);
    out.health = zeros(nn,GPS_SIZE);
    out.prn = [];
    if params.doH == 1 % Following only needed for Jacobian computation
        out.eciRotation = zeros(9,9,nn);
        out.los_3d = zeros(3,nn,GPS_SIZE);
        out.gps_vel_tot = zeros(3,nn,GPS_SIZE);
        out.sat_vel_tot = zeros(3,nn);
        out.Rotation2ECI = Rotation2ECI;
    end
    return;
end

% Compute satellite orbits from almanac parameters
% Only hand in almanac parameters for healthy SVs
% Initialize variables
gps_pos = zeros(3,nn,GPS_SIZE);
gps_vel = zeros(3,nn,GPS_SIZE);
gps_vel_tot = gps_vel;

% GPS position and velocity in and relative to ECEF frame
[pos,vel,vel_tot] = alm2xyz(time(:)',gps_alm(healthy_ind,:),1);	% [3,nn,mm] arrays
if GPS_SIZE == 1
    % one SV calc only, PRN is in out.prn
    gps_pos(:,:,1) = pos; %km
    gps_vel(:,:,1) = vel; %km, omega-cross term removed
    gps_vel_tot(:,:,1) = vel_tot;
else
    % all 32 SVs, index PRN via location in GPS_SIZE
    gps_pos(:,:,prns) = pos; %km
    gps_vel(:,:,prns) = vel; %km, omega-cross term removed
    gps_vel_tot(:,:,prns) = vel_tot;
end

%% Compute Visibility

health = reshape(max(gps_pos),nn,GPS_SIZE) ~= 0;        % [nn,GPS_SIZE]
rgps_mag = reshape(sqrt(sum(gps_pos.^2)),nn,GPS_SIZE);    % [nn,GPS_SIZE]

ruser_3d = zeros(3,nn,GPS_SIZE);
vuser_3d = zeros(3,nn,GPS_SIZE);
vuser_3d_tot = zeros(3,nn,GPS_SIZE);
for i=1:GPS_SIZE
    ruser_3d(:,:,i) = sat_pos;                     % [3,nn,GPS_SIZE]
    vuser_3d(:,:,i) = sat_vel;                     % [3,nn,GPS_SIZE]
    vuser_3d_tot(:,:,i) = sat_vel_tot;                     % [3,nn,GPS_SIZE]
end

los_3d = gps_pos - ruser_3d;                          % [3,nn,GPS_SIZE]
los_mag = reshape(sqrt(sum(los_3d.^2)),nn,GPS_SIZE);     % [nn,GPS_SIZE]
los_mag_3d = zeros(3,nn,GPS_SIZE);
rgps_mag_3d = zeros(3,nn,GPS_SIZE);
for i=1:3
    los_mag_3d(i,:,:) = los_mag;                   % [3,nn,GPS_SIZE]
    rgps_mag_3d(i,:,:) = rgps_mag;                 % [3,nn,GPS_SIZE]
end
los_unit_3d = los_3d./los_mag_3d;                  % [3,nn,GPS_SIZE]

%  Range Rate in km/s (relative velocity projected along the LOS)
%  DOES NOT ACCOUNT FOR SECOND ORDER DOPPLER EFFECT DUE TO RELATIVITY
rrate = reshape(dot(los_unit_3d,(gps_vel_tot - vuser_3d_tot)),nn,GPS_SIZE); % (nn,GPS_SIZE)

%% GPS satellites
if(params.xmit_pattern_dim == 2)
    
    % transmitter azimuth and elevation based on a GPS SV yaw model
    [TX_az,TX_el,GPS_yaw] = alm_pos_az_el(epoch,t,sat_pos,gps_pos, gps_vel);
    TX_az = TX_az'; % (nn,GPS_SIZE)
    TX_el = TX_el'; % (nn,GPS_SIZE)
    
else
    
    TX_az = [];
    GPS_yaw = [];
    %--- Compute elevation angles wrt GPS SV antenna
    % SV off-boresite angle - this is the transmitter 'elevation' angle
    % measured between los and the SV antenna boresite (nadir)
    denom=rgps_mag.*los_mag;
    denom(denom==0)=NaN;
    TX_el = acos(reshape(dot(gps_pos,los_3d),nn,GPS_SIZE)./denom)*r2d;    % (nn,GPS_SIZE)

end

if ~exist('qatt','var') || isempty(qatt)
    wstr = ['S/C body quaternion is not specified in the input.  ',...
        'Body will be assumed to be aligned with the same coordinate ',...
        'frame as the position and velocity states.'];
    warning(wstr); %#ok<WNTAG>
    qatt = repmat([0;0;0;1],1,nn);
end
if ~exist('ref2body','var')
    [mq,nq] = size(qatt);
    if mq ~= 4 || nq ~= nn
        error('In GPSMEAS, the specified spacecraft attitude quaternion does not have the right dimension(s).');
    end
    ref2body = q2dcm(qatt);
end

RX_az   = NaN(nn,GPS_SIZE,loop);
RX_el   = NaN(nn,GPS_SIZE,loop);

% Execute loop at least once
for ANT=1:loop
    
    if(params.rec_pattern_dim(ANT) == 2) % 2D antenna
        
        % Transform los from ecef frame to receiver antenna frame
        for j = 1:GPS_SIZE
            los_ant = zeros(3,nn);
            for i = 1:nn
                R2ECI = Rotation2ECI(t(i)+epoch);
                % los_ant = (antenna <- body <- state frame <- ECI <- ECEF) * los_ECEF
                los_ant(:,i) = ant_body(:,:,ANT) * ref2body(:,:,i) * R2ECI(1:3,1:3)' * ...
                    init2fixed(:,:,i)' * los_unit_3d(:,i,j);
            end
            RX_az(:,j,ANT) = atan2(los_ant(2,:),los_ant(1,:))*r2d;
            RX_el(:,j,ANT) = 90 - asin(los_ant(3,:));
            
        end
        
    else % 1D antenna
        
        % Compute matrix describing antenna boresite(s) [3,nn]
        boresite = comp_bs_3d(1, time, sat_pos, sat_vel, pointing_ref(ANT), 1);  % [3,nn]
        boresite_3d = reshape(repmat(boresite,1,GPS_SIZE),3,nn,GPS_SIZE);
        RX_el(:,:,ANT) = (abs(acos(reshape(dot(boresite_3d,los_unit_3d),nn,GPS_SIZE))))*r2d;  % (nn,GPS_SIZE)

    end
    
end

%% Set output structure
out.epoch = epoch;
out.TX_az = TX_az;
out.TX_el = TX_el;
out.RX_az = RX_az;
out.RX_el = RX_el;
out.range = los_mag';
out.rrate = rrate';
out.GPS_yaw = GPS_yaw;
out.rgps_mag = rgps_mag;
out.health = health;
if GPS_SIZE == 1
    out.prn = prns;
else
    out.prn = zeros(1,32);
    out.prn(prns) = 1;
end
if params.doH == 1 % Following only needed for Jacobian computation
    out.eciRotation = eciRotation;
    out.los_3d = los_3d;
    out.gps_vel_tot = gps_vel_tot;
    out.sat_vel_tot = sat_vel_tot;
    out.Rotation2ECI = Rotation2ECI;
end

% Generic template for Rotation2ECI(t).
% Returns a 6x6 transformation matrix that transforms a position and
% velocity in a given frame to ECI at a specified time, t.  Time is
% specified in MATLAB datenum format.
function I=IdentRot(t) %#ok<INUSD>
    % This template function always returns an identity transformation.
    I = eye(6);



