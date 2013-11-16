function [y,H,R,AntLB] = gsmeas(t,x,options,qatt)
% GSMEAS  Makes ground station based measurements.
%
% [y,H,R] = GSMEAS(tspan,x,options) creates ground station measurements
% based on the information in OPTIONS. The measurement types that can be 
% returned are LOSRANGE, LOSRANGERATE, LOSDOPPLER. See the OD Toolbox 
% function of the same name for details of each measurement type.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch)
%      x            (6xN)   ECI J2000 spacecraft state [pos;vel] (km)
%      options      (1x1)   data structure, see below
%
%   OUTPUTS
%      y            (MxN)   measurements
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each groundstation are grouped together in
% rows. Thus, using the default options settings and passing in 3 ground
% stations, the output y will look like:
%   y = [   range_gs1(t1)       range_gs1(t2)...;
%           rangeRate_gs1(t1)   rangeRate_gs1(t2)...;
%           range_gs2(t1)       range_gs2(t2)...;
%           rangeRate_gs2(t1)   rangeRate_gs2(t2)...;
%           range_gs3(t1)       range_gs3(t2)...;
%           rangeRate_gs3(t1)   rangeRate_gs3(t2)  ]
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
% For station locations based on information in the NASA Directory Of
% Station Locations (NDOSL), set these parameters:
%   PARAMETER           VALID VALUES         NOTES
%   gsList                JAT object         See CREATEGROUNDSTATIONLIST
%   gsID                  See NDOSL file     Cell array of ground station
%                                             IDs in the  NDOSL
%
% For faster perfomance (or for unlisted ground stations), predefine the 
% ECEF station locations as columns of a matrix using this parameter:
%   PARAMETER           VALID VALUES         NOTES
%   gsECEF                [(3xNs) matrix]    ECEF locations of ground
%                                             stations (km). Ns is the 
%                                             number of stations.
%
% The following parameters apply to both location input types:
%   PARAMETER           VALID VALUES         NOTES
%   epoch                 datenum            UTC time associated with 
%                                             start of simulation
%   gsElevationConstraint degs               Spacecraft elevation must be
%                                             above this value to return
%                                             measurement
%   useRange          {true(default), false} Return range measurement
%   rangeType         {'1way','1wayFWD','1wayRTN','2way'(default)} 
%                                             Note: 1way=1wayFWD
%   useRangeRate      {true(default), false} Return range rate measurement
%   useDoppler        {true, false(default)} Return doppler measurement
%   useUnit           {true, false(default)} Return unit vector measuerment
%   frequencyTransmit {scalar>0, 1.57542e9}  Hz, Only used for Tropo and 
%                                             for Doppler
%   rSigma            {(1xM),ones(1,M)*1e-3) Measurement covariance
%   useLightTime      {true, false(default)} Include light time delay
%   useGPSIonosphere  {true, false(default)} Includes GPS Ionospheric delay
%   useIonosphere     {true, false(default)} Includes Ionospheric delay
%   useTroposphere    {true, false(default)} Includes Tropospheric delay
%   Schedule          [(ncx3) matrix]        Ground Tracking Schedule. 
%      Schedule restricts the gsmeas measurement model to only provide  
%      measurements for specific ground stations during specific time  
%      intervals. Theformat for each row is:
%      [ gs_index start_time stop_time ].
%      The gs_index corresponds to index of the ground stations in gsID or 
%      gsECEF. The start_time and stop_time must be in seconds from epoch.
%      The matrix can be any length (nc = number of contacts)
%      Example: 
%       TrackSched = {1 '28-Jan-2010 06:56:35' '28-Jan-2010 12:56:35'
%                     2 '28-Jan-2010 12:56:35' '28-Jan-2010 19:56:35'
%                     1 '28-Jan-2010 19:56:35' '29-Jan-2010 02:56:35'};
%       Sched = cell2mat(TrackSched(:,1)); %ground station numbers
%       Sched(:,2) = (datenum(TrackSched(:,2))-epoch)*86400; %start(epsecs)
%       Sched(:,3) = (datenum(TrackSched(:,3))-epoch)*86400; %end(epsecs)
%
% The options parameters associated with each of the desired measurement
% types are passed to the appropriate function.
%
% The ground stations to use are input as part of OPTIONS. The
% groundstation list must also be input as part of OPTIONS. This list can
% be created by CREATEGROUNDSTATIONLIST.
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, ODTBXOPTIONS,
% CREATEGROUNDSTATIONLIST
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
%   Kevin Berry         09/07/2007      Original
%   Derek Surka         09/17/2007      Modified to use options input and 
%                                        elevation constraint
%   Kevin Berry         04/10/2008      Added Doppler functionality
%   Kevin Berry         05/21/2008      Simplified with my new rrdot and
%                                        rrdotlt functions
%   Allen Brown         02/25/2009      Updated documentation
%   Kevin Berry         06/25/2009      Added time scale comments
%   Kevin Berry         07/01/2009      Added angle based measurements
%   Kevin Berry         09/03/2009      Changed angle measurement to Unit
%   Kevin Berry         02/02/2010      Added the option to input station
%                                         locations in ECEF instead of 
%                                         using the NDOSL list at every 
%                                         time step
%   Russell Carpenter   02/11/2011      Added useAngles option

d2r          = pi/180;
r2d          = 180/pi;
%% Get values from options
gsID         = getOdtbxOptions(options, 'gsID', [] );
gsList       = getOdtbxOptions(options, 'gsList', []);
gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
elMin        = getOdtbxOptions(options, 'gsElevationConstraint', 10)*pi/180; % convert from degs to rads
uselt        = getOdtbxOptions(options, 'useLightTime', false);
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
useUnit      = getOdtbxOptions(options, 'useUnit', false );
useAngles    = getOdtbxOptions(options, 'useAngles', false );
Sched        = getOdtbxOptions(options, 'Schedule',[]); %Tracking Schedule
numtypes     = useRange + useRangeRate + useDoppler+3*useUnit+2*useAngles;

numGS = length(gsID);
if isempty(gsECEF)
    if( isempty(gsList) && ~isempty(gsID) ); gsList = createGroundStationList(); end
    gsECEF = zeros(3,length(gsID));
    for n=1:numGS
            gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
    end
end
M            = numGS * numtypes;
N            = length(t);
if size(t,1)==N, t=t'; end

if isnan(epoch); error('An epoch must be set in the options structure.'); end

if isfield(options,'linkbudget') && ~isempty(options.linkbudget)
    dolinkbudget = true;
    link_budget = getOdtbxOptions(options, 'linkbudget',[]);

    % Set some default values
    link_budget = linkbudget_default(link_budget, 'AntennaPattern', {'omni.txt'});
        %  Specify antenna pattern for each antenna, existing antennas are:
        %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
        %     omni.txt                  - zero dB gain,  180 degree half beamwidth
        %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
        %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
        %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
    num_ant = length(link_budget.AntennaPattern); %hasn't been tested for >4 antennas
    link_budget = linkbudget_default(link_budget, 'RXAntennaMask', 180*d2r);
    link_budget = linkbudget_default(link_budget, 'AtmosphereMask', 0); % km
        %  Troposphere mask radius ~50 km
        %  Ionosphere mask radius ~(500-1000 km)
    link_budget = linkbudget_default(link_budget, 'NoiseTemp', 300); % K
        % System noise temp [K], space pointing antenna = 290
        % System noise temp [K], earth pointing antenna = 300
    link_budget = linkbudget_default(link_budget, 'AtmAttenuation', 0.0); % dB
        % attenuation due to atmosphere (should be negative) [dB]
    link_budget = linkbudget_default(link_budget, 'TransPowerLevel', 2); % 1-minimum, 2-typical, 3-max
    link_budget = linkbudget_default(link_budget, 'TransPowerOffset', 0.0); % dB, global offset
    link_budget = linkbudget_default(link_budget, 'TXAntennaMask', 70*d2r );  % in rad
        %  The actual mask used is the lesser of this mask and the limit of the defined pattern
        %  Note:  mask = 70 deg includes entire defined pattern
        %         mask = 42 deg includes only main and first side lobes
        %         mask = 26 deg includes only main lobe
    link_budget = linkbudget_default(link_budget, 'ReceiverNoise', -3 );  % dB, Noise figure of receiver/LNA
    link_budget = linkbudget_default(link_budget, 'RecConversionLoss', -1.5 );  % dB
        % Receiver implementation, A/D conversion losses [dB]
        %   Novatel: L = -4.0 	
        %   Plessey: L = -1.5		
    link_budget = linkbudget_default(link_budget, 'SystemLoss', 0 ); % dB, System losses, in front of LNA
    link_budget = linkbudget_default(link_budget, 'LNAGain', 40 ); % dB, LNA gain (trimble pre-amp spec = 42-48)
    link_budget = linkbudget_default(link_budget, 'CableLoss', -2 ); % dB, Cable losses (after LNA)
    link_budget = linkbudget_default(link_budget, 'RecAcqThresh', 32 ); % dB-Hz, Receiver acquisition threshold
    link_budget = linkbudget_default(link_budget, 'RecTrackThresh', link_budget.RecAcqThresh ); % dB-Hz, Receiver tracking threshold
    link_budget = linkbudget_default(link_budget, 'DynamicTrackRange', 15 ); % dB
        % Dynamic tracking range of receiver, or maximum difference  
        % in power levels tracked simultaneously. If the difference
        % in snrs between two satellites is more that link_budget.DynamicTrackRange,
        % the weaker of the two will not be considered visible. 
    % Reassign the options structure with any changed/default link budget values
    options = setOdtbxOptions(options, 'linkbudget', link_budget);
else
    dolinkbudget = false;
end


%% call rrdot with all the options passing straight through
y    = nan(M,N);
H    = zeros(M,6,N);  
if uselt
    %x2 needs states before and after each state in x1 at time steps
    %comparable to the light time delay in order for the interpolation
    %within the lightTimeCorrection function to have a high level of
    %accuracy.
    c    = JATConstant('c')/1000;
    ltDT = sqrt(sum(x(1:3,:).^2))/c; %
    t2   =  unique([t-ltDT, t, t+ltDT]);

    % Get the ground station positions at times t2 (see subfunction below)
    x2 = getECIstates([gsECEF;zeros(size(gsECEF))],epoch,t2);

    % Run rrdotlt for each ground station
    for n=1:size(gsECEF,2)
        % Check schedule for times when tracking is done
        if ~isempty(Sched)
            gSch = Sched(n==Sched(:,1),2:3);
            tind = [];
            for m=1:size(gSch,1)
                tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
            end
        else
            tind = 1:length(t);
        end
        t1  = t(tind);
        x1  = x(:,tind);
        
        if ~isempty(t1)

            [y1,H1,R,t2_lt,x2_lt] = rrdotlt(t1,x1,t2,x2(:,:,n),options);
            
            % apply the elevation constraint
            Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
            Ephem.SatCoords   = 'ECI';
            Ephem.Epoch       = epoch+t2_lt{1}/86400; %UTC
            Ephem.StationInfo = 'ECEF';
            Ephem.staPos      = gsECEF(:,n)*1000;
            [az,el]            = jatStaAzEl(Ephem);
            index0            = find( el < elMin );
            if dolinkbudget
                %add TX_az, TX_el to out
                if n==1
                    out.TX_az = zeros(N,numGS);
                    out.TX_el = zeros(N,numGS);
                end
                out.TX_az(:,n) = az*180/pi;
                out.TX_el(:,n) = -(el*180/pi-90);
            end
            if length(x2_lt)==2 %then it was a 2way measurement
                Ephem.Epoch = epoch+t2_lt{2}/86400; %UTC
                [~,el]      = jatStaAzEl(Ephem);
                index1      = find( el < elMin );
                index0      = union(index0,index1);
            end
            y1(:,index0) = NaN;

            % combine with results from previous stations
            indstart                     = 1 + numtypes*(n-1);
            indstop                      = numtypes*n;
            y(indstart:indstop, tind)    = y1;
            H(indstart:indstop, :, tind) = H1;
        end
    end
    clear R;
else
    % Get the ground station positions at times t (see subfunction below)
    gx = getECIstates([gsECEF;zeros(size(gsECEF))],epoch,t);

    % Run rrdot for each ground station
    for n=1:size(gsECEF,2)
        % Check schedule for times when tracking is done
        if ~isempty(Sched)
            gSch = Sched(n==Sched(:,1),2:3);
            tind = [];
            for m=1:size(gSch,1)
                tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
            end
        else
            tind = 1:length(t);
        end
        if isempty(tind),continue,end
        t1 = t(tind);
        x1 = x(:,tind);
        x2 = gx(:,tind,n);

        [y1,H1] = rrdotang(t1,x1,x2,options);

        % apply the elevation constraint
        Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
        Ephem.SatCoords   = 'ECI';
        Ephem.Epoch       = epoch+t1/86400;%UTC
        Ephem.StationInfo = 'ECEF';
        Ephem.staPos      = gsECEF(:,n)*1000;
        [az,el]            = jatStaAzEl(Ephem);        
        y1(:,el<elMin)    = NaN;
        if dolinkbudget
            %add TX_az, TX_el to out
            if n==1
                out.TX_az = zeros(N,numGS);
                out.TX_el = zeros(N,numGS);
            end
            out.TX_az(:,n) = az*180/pi;
            out.TX_el(:,n) = -(el*180/pi-90);
        end
        % combine with results from previous stations
        indstart                     = 1 + numtypes*(n-1);
        indstop                      = numtypes*n;
        y(indstart:indstop, tind)    = y1;
        H(indstart:indstop, :, tind) = H1;
    end
end

%% Set the measurement covariance output
if nargout > 2,
    %get sigma out of the options
    sigmaDefault = ones(1,size(y,1))*1e-3;
    sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
    if( length(sigma)~=length(sigmaDefault) )
        disp('WARNING: In gsmeas, length(rSigma) does not match total number of measurements');
        disp('   Setting rSigma = default (1e-3) for all measurements');
        sigma = sigmaDefault;
    end
    sigma           = diag(sigma);
    R = repmat(sigma.^2,[1,1,N]);
end

%% Perform link budget analysis
if dolinkbudget
    % Generate link budget structures
    pointing_ref    = getOdtbxOptions(options, 'AntennaPointing', [-1 1] );
    % Specify attitude profile for each antenna
    %   (1) zenith pointing or (-1) nadir pointing wrt geocentric LVLH
    %   (2) parallel or (-2) antiparallel to the Earth-Sun vector
    %   (3) ecliptic north or south (-3)
    %   (4) fixed with respect to apogee zenith, or (-4) nadir
    %   vector apogee
    %       selected as the point of highest altitude (ephemeris must include apogee)
    %   (5) body fore and (-5) aft directions relative to geocentric LVLH
    %   (6) body port and (-6) starboard directions relative to geocentric LVLH
    
    % set link budget params in structs
    TX_link.P_sv = link_budget.TransPowerOffset;
    % Set loop to num_ant or 1
    loop = max([1,num_ant]);
    % Transmitter and Receiver antenna patterns
    TX_link.pattern = load(link_budget.TXpattern);
    RX_link.pattern = cell(loop,1);
    rec_pattern_dim = ones(loop,1);
    for ANT = 1:loop
        RX_link.pattern{ANT} = load(link_budget.AntennaPattern{ANT});
        if size(RX_link.pattern{ANT},2) > 2
            rec_pattern_dim(ANT) = 2;
        end
    end 
    
    ant_body        = getOdtbxOptions(options, 'AntennaOrientation',reshape(repmat(eye(3),1,num_ant),3,3,num_ant));
    


    %add rgps_mag, health to out
    gsECEF_overTime=repmat(reshape(gsECEF,3,1,numGS),[1 N 1]);
    out.rgps_mag = reshape(sqrt(sum(gsECEF_overTime.^2)),N,numGS);
    out.health = reshape(max(gsECEF_overTime),N,numGS) ~= 0;
    
    %add RX_az, RX_el to out
    
        % Calculate x position in ECEF
        [xECEF,init2fixed] = getECEFstates(x,epoch,t); 
        sat_vel_tot = zeros(3,N);
        for ii=1:N
            R2ECI = eye(6);  % Unused, set at Identity
            sat_vel_tot(:,ii) = init2fixed(:,:,ii) * R2ECI(4:6,:) *  x(1:6,ii);
        end
        %define attitude quaternion and DCM
        if nargin < 4
            wstr = ['S/C body quaternion is not specified in the input.  ',...
                'Body will be assumed to be aligned with the same coordinate ',...
                'frame as the position and velocity states.'];
            warning('ODTBX:GSMEAS:noBodyQuat',wstr);
            qatt = repmat([0;0;0;1],1,N);
        end
        [mq,nq] = size(qatt);
        if mq ~= 4 || nq ~= N
            error('In GPSMEAS, the specified spacecraft attitude quaternion does not have the right dimension(s).');
        end
        ref2body = q2dcm(qatt);
    
        % Define LOS unit vector
        ruser_3d = zeros(3,N,numGS);
        vuser_3d = zeros(3,N,numGS);
        vuser_3d_tot = zeros(3,N,numGS);

        for i=1:numGS
            ruser_3d(:,:,i) = xECEF(1:3,:);                     % [3,N,numGS]
            vuser_3d(:,:,i) = xECEF(4:6,:);                     % [3,N,numGS]
            vuser_3d_tot(:,:,i) = sat_vel_tot;                     % [3,N,numGS]
        end

        los_3d = gsECEF_overTime - ruser_3d;                 % [3,N,numGS]
        los_mag = reshape(sqrt(sum(los_3d.^2)),N,numGS);     % [N,numGS]
        los_mag_3d = zeros(3,N,numGS);
        rgps_mag_3d = zeros(3,N,numGS);
        for i=1:3
            los_mag_3d(i,:,:) = los_mag;                   % [3,N,numGS]
            rgps_mag_3d(i,:,:) = out.rgps_mag;                 % [3,N,numGS]
        end
        los_unit_3d = los_3d./los_mag_3d;                  % [3,N,numGS]

    % Define range, az and el
    out.range = los_mag';
    out.RX_az    = NaN(N,numGS,num_ant);
    out.RX_el    = NaN(N,numGS,num_ant);

    for ANT=1:num_ant

        if(size(RX_link.pattern,2) >= 2) % 2D antenna

            % Transform los from ecef frame to receiver antenna frame
            for j = 1:numGS
                los_ant = zeros(3,N);
                for i = 1:N
                    % los_ant = (antenna <- body <- state frame <- ECI <- ECEF) * los_ECEF
                    los_ant(:,i) = ant_body(:,:,ANT) * ref2body(:,:,i) * R2ECI(1:3,1:3)' * ...
                        init2fixed(:,:,i)' * los_unit_3d(:,i,j);
                end
                out.RX_az(:,j,ANT) = atan2(los_ant(2,:),los_ant(1,:))*r2d;
                out.RX_el(:,j,ANT) = 90 - asin(los_ant(3,:))*r2d;

            end

        else % 1D antenna
            %  Convert simulation start time from UTC to GPS time in seconds
            time = 86400 * convertTime('GPS','UTC', epoch+t/86400);
            % Compute matrix describing antenna boresite(s) [3,N]
            boresite = comp_bs_3d(1, time, xECEF(1:3,:), xECEF(4:6,:), pointing_ref(ANT), 1);  % [3,N]
            boresite_3d = reshape(repmat(boresite,1,numGS),3,N,numGS);
            out.RX_el(:,:,ANT) = (abs(acos(reshape(dot(boresite_3d,los_unit_3d),N,numGS))))*r2d;  % (N,numGS)

        end

    end
    % Calculate Link Budget
    [AntLB, HVIS] = calc_linkbudgets(out, options, RX_link, TX_link);

    % Make measurements that are not visible NaNs
    kk=1;
    for ii=1:size(HVIS,1)
        for jj=1:size(HVIS,2)
            if HVIS(ii,jj)==0
                y(kk:kk+numtypes-1,jj) = NaN;
            end
        end
        kk=kk+numtypes;
    end
end

end


%% Subfunction for getting ECI states of ground stations at times t
function [xECI,D] = getECIstates(xECEF,epoch,t)

D = jatDCM('ecef2eci', epoch+t/86400);
w = [0;0;JATConstant('wEarth')];

xECI = zeros(6,length(t),size(xECEF,2));
for nt = 1:length(t);
    M          = rotransf(-D(:,:,nt)*w,D(:,:,nt));
    xECI(:,nt,:) = M(1:6,1:6)*xECEF;
end
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