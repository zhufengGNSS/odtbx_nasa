function [ AntLB, HVIS ] = get_linkbudget(t, x1, x2, options, qatt )
% Calculate a link budget between the receiver x1 with num_ant antennas and
% the transmitter x2.  
%
%[ AntLB, HVIS ] = get_linkbudget(t, x1, x2, options, qatt )
%
% INPUTS
%   t           = (vector of length N) times of the measurements
%   x1          = (6xN) receiver position and velocity in ECI (km or km/s)
%   x2          = (6xN) transmitter position and velocity in ECI (km or km/s)
%   options     = (1x1) OD Toolbox Measurement Options data structure. See
%                   ODTBXOPTIONS for all available options settings.
%   qatt        = (4xN)   spacecraft (x1) body quaternion history where the
%                   quaternion specifies the user spacecraft
%                   attitude with respect to the same coordinate
%                   frame that the position and velocity are given
%                   in.  The quaternion is specified as
%                   [sin(theta/2)*e_vec;cos(theta/2)] where e_vec is the
%                   unit vector representing the axis of rotation
%                   and theta is the angle of rotation about that
%                   axis.  The quaternion information is used IF a 
%                   2-D antenna is specified for the user satellite.  
%                   Additionally, the orientation matrix
%                   of each antenna wrt to the spacecraft body can 
%                   be specified through the options strucutre.  If
%                   the antennae orientation is not specified, the
%                   default is identity, i.e., aligned with the body
%                   axes.
%   OUTPUTS
%      AntLB        {num_antx1}  Cell array containing link budget for each
%                           antenna.  Note that "masked" refers to an applied
%                           bias due to visibility or blockage.  The fields
%                           are:
%           Halpha_r    (MxN)   The receiver elevation angle (rad)
%           Hvis_beta   (MxN)   Logical array where both transmit and
%                               receive antenna elevation angles are
%                               within antenna mask angle limits
%           Hvis_CN0    (MxN)   Logical array where CN0 is above the
%                               acquisition/tracking threshold
%           HCN0        (MxN)   Masked Signal carrier to noise ratio
%           HAd         (MxN)   Masked attenuation from R^2 losses (dB)
%           HAr         (MxN)   Receive antenna gain (dB)
%           HAP         (MxN)   Masked budget gain before receiver antenna (dB)
%           HRP         (MxN)   Masked budget gain before receiver amplifiers
%                               and conversion (dB)
%           HAt         (MxN)   Masked transmit antenna gain (dB)
%       HVIS
%
%  keywords
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
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
%   Jason Schmidt      11/21/2013       Initial Version

%% Extract values from inputs
N=length(t);
epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
link_budget = getOdtbxOptions(options, 'linkbudget',[]);
if isempty(link_budget)
    error('link_budget must be defined')
end
num_ant = length(link_budget.AntennaPattern);
%define DCM for x1
[mq,nq] = size(qatt);
if mq ~= 4 || nq ~= N
    error('ODTBX:GSMEAS:badQuat','In GPSMEAS, the specified spacecraft attitude quaternion does not have the right dimension(s).');
end
ref2body = q2dcm(qatt); % Used by 2D antenna
ant_body = getOdtbxOptions(options, 'AntennaOrientation',reshape(repmat(eye(3),1,num_ant),3,3,num_ant)); % Used by 2D antennas
pointing_ref    = getOdtbxOptions(options, 'AntennaPointing', [-1 1] ); % Used by 1D antennas
% Specify attitude profile for each antenna
%   (1) zenith pointing or (-1) nadir pointing wrt geocentric LVLH
%   (2) parallel or (-2) antiparallel to the Earth-Sun vector
%   (3) ecliptic north or south (-3)
%   (4) fixed with respect to apogee zenith, or (-4) nadir
%   vector apogee
%       selected as the point of highest altitude (ephemeris must include apogee)
%   (5) body fore and (-5) aft directions relative to geocentric LVLH
%   (6) body port and (-6) starboard directions relative to geocentric LVLH

%% Calculate x1 position in ECEF
[x1ECEF,init2fixed] = getECEFstates(x1,epoch,t);

%% Calculate x2 position in ECEF
[x2ECEF,~] = getECEFstates(x2,epoch,t);

%% Define Transmitter antenna Pattern and other attributes
TX_link.pattern = load(link_budget.TXpattern);
TX_link.P_sv = link_budget.TransPowerOffset;

%% Define the Receiver antenna Pattern and other attributes
% Set loop to num_ant or 1
loop = max([1,num_ant]);
rec_pattern_dim = ones(loop,1);
for ANT = 1:loop
    RX_link.pattern{ANT} = load(link_budget.AntennaPattern{ANT});
    if size(RX_link.pattern{ANT},2) > 2
        rec_pattern_dim(ANT) = 2;
    end
end 

%% Calculate the az and el of the signal wrt the transmitter antenna
LOS_TX2RX_3d   = x1ECEF(1:3,:)-x2ECEF(1:3,:);
[LOS_RX2TX_unit,out.range] = unit(-LOS_TX2RX_3d);
if link_budget.TX_AntennaPointing == 1 % zenith pointing TX
    [out.TX_el] = angleBetweenVectors(x2ECEF(1:3,:),LOS_TX2RX_3d)';
    el_local_horizon = 90-out.TX_el;
elseif link_budget.TX_AntennaPointing == -1 % nadir pointing TX
    [out.TX_el] = angleBetweenVectors(-x2ECEF(1:3,:),LOS_TX2RX_3d)';
    el_local_horizon = out.TX_el-90;
else
    error('link_budget.TX_AntennaPointing must be either 1 for zenith pointing or -1 for nadir pointing');
end
% calculate azimuth
% get Topocentric Horizon Coordinate Rotations for x2ECEF
[~,T_ECEF2TopoHoriz] = getTopoHorizStates(x2ECEF(1:3,:));
% get Topocentric Horizon version of LOS_RX2TX_unit
LOS_RX2TX_unit_TopoHoriz = zeros(size(LOS_RX2TX_unit));
for ii = 1:N
    LOS_RX2TX_unit_TopoHoriz(:,ii) = T_ECEF2TopoHoriz(:,:,ii)*LOS_RX2TX_unit(:,ii);
end

% Calculate azimuth assuming az if from north going clockwise (N,E,S,W)
sinAz = -LOS_RX2TX_unit_TopoHoriz(1,:)./cosd(el_local_horizon)';
out.TX_az = acosd(-LOS_RX2TX_unit_TopoHoriz(2,:)./cosd(el_local_horizon)');
out.TX_az(sinAz<0) = -out.TX_az(sinAz<0);


%% Calculate the az and el of the signal wrt the receiver antenna
out.RX_az    = NaN(N,1,num_ant);
out.RX_el    = NaN(N,1,num_ant);

for ANT=1:num_ant

    if(size(RX_link.pattern{ANT},2) > 2) % 2D antenna

        % Transform los from ecef frame to receiver antenna frame
        los_ant = zeros(3,N);
        for i = 1:N
            % los_ant = (antenna <- body <- state frame <- ECI <- ECEF) * los_ECEF
            los_ant(:,i) = ant_body(:,:,ANT) * ref2body(:,:,i) * ...
                init2fixed(:,:,i)' * LOS_RX2TX_unit(:,i);
        end
        out.RX_az(:,1,ANT) = atan2(los_ant(2,:),los_ant(1,:))*(180/pi);
        out.RX_el(:,1,ANT) = 90 - asin(los_ant(3,:))*(180/pi);

    else % 1D antenna
        %  Convert simulation start time from UTC to GPS time in seconds
        time = 86400 * convertTime('GPS','UTC', epoch+t/86400);
        % Compute matrix describing antenna boresite(s) [3,N]
        boresite = comp_bs_3d(1, time, x1ECEF(1:3,:), x1ECEF(4:6,:), pointing_ref(ANT), 1);  % [3,N]
        boresite_3d = reshape(repmat(boresite,1,1),3,N,1);
        out.RX_el(:,:,ANT) = (abs(acos(reshape(dot(boresite_3d,LOS_RX2TX_unit),N,1))))*(180/pi);  % (N,numGS)

    end

end

%% Calculate other required values (range, health, etc)
[~,out.rgps_mag] = unit(x2ECEF(1:3,:));
out.rgps_mag = out.rgps_mag';
out.health = reshape(max(x2ECEF),N,1) ~= 0;
%% Calculate Link Budget
[AntLB, HVIS] = calc_linkbudgets(out, options, RX_link, TX_link); 
end


%% Subfunction for getting ECEF states of sat at times t
function [xECEF,D] = getECEFstates(xECI,epoch,t)

    D = jatDCM('eci2ecef', epoch+t/86400);
    w = [0;0;JATConstant('wEarth')];

    xECEF = zeros(6,length(t));
    for nt = 1:length(t);
        M          = rotransf(D(:,:,nt)*w,D(:,:,nt));
        xECEF(:,nt) = M(1:6,1:6)*xECI(:,nt);
    end
end

%% Subfunction for getting angle between vectors
function [angle] = angleBetweenVectors(A,B)
% output in degrees
    [m,n]=size(A);
    if m~=size(B,1) || n~=size(B,2)
        error('The two vectors/matrices must be the same size')
    end

    if xor(m==1,n==1) && xor(m==3, n==3) % Treat like two unit vectors 
        angle = atan2d(norm(cross(A,B)),dot(A,B));
    elseif xor(m==3,n==3) % Treat like a collection of vectors
        crossAB = cross(A,B);
        normcrossAB = sqrt(sum(crossAB.^2));
        angle = atan2d(normcrossAB,dot(A,B));
    elseif m==3 && n==3
        warning('Both dimensions of A and B are 3. Assuming matrix of column vectors')
        crossAB = cross(A,B);
        normcrossAB = sqrt(sum(crossAB.^2));
        angle = atan2d(normcrossAB,dot(A,B));
    else
        error('The two vectors/matrices must be unit vectors or a collection of unit vectors');
    end

end

%% Subfunction for getting Topocentric horizon states of sat at times t
function [xTopoHoriz,T_ECEF2TopoHoriz] = getTopoHorizStates(xECEF)

% convert x2ECEF into Geocentric Lat/Long 
[xECEF_unit, ~] = unit(xECEF(1:3,:));
longitude = atan2(xECEF_unit(2,:),xECEF_unit(1,:));
latitude = asin(xECEF_unit(3,:));

% Pre calculate sin and cosines
sinLong = sin(longitude);
cosLong = cos(longitude);
sinLat = sin(latitude);
cosLat = cos(latitude);

N = size(xECEF,2);
xTopoHoriz = zeros(3,N);
T_ECEF2TopoHoriz = zeros(3,3,N);
for ii=1:N
    
    T_ECEF2TopoHoriz(:,:,ii) = [-sinLong(ii)                cosLong(ii)                 0
                        -sinLat(ii)*cosLong(ii)     -sinLat(ii)*sinLong(ii)     cosLat(ii)
                        cosLat(ii)*cosLong(ii)      cosLat(ii)*sinLong(ii)      sinLat(ii)];
    xTopoHoriz(:,ii) = T_ECEF2TopoHoriz(:,:,ii)*xECEF(:,ii);
end
          
end