function [AntLB, HVIS] = getlinkbudget(out, options, RX_link, TX_link)
%% Process incoming data
d2r             = pi/180;

% Get Some Constants from JAT
EARTH_RADIUS = JATConstant('rEarth','WGS84') / 1000;  % km Equatorial radius of Earth
% c            = JATConstant('c') / 1000;               % km/sec speed of light

% Data from odtbxOptions structure
rec_pattern     = getOdtbxOptions(options, 'AntennaPattern', {'sensysmeas_ant.txt','sensysmeas_ant.txt'} );
                %  Specify antenna pattern for each antenna, existing antennas are:
                %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
                %     omni.txt                  - zero dB gain,  180 degree half beamwidth
                %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
                %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
                %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
num_ant         = length(rec_pattern); %hasn't been tested for >4 antennas
% Ts              = getOdtbxOptions(options, 'NoiseTemp', 300 ); % K
AtmMask         = getOdtbxOptions(options, 'AtmosphereMask', 50 ); % km 
                %  Troposphere mask radius ~50 km
                %  Ionosphere mask radius ~(500-1000 km)

                % System noise temp [K], space pointing antenna = 290
                % System noise temp [K], earth pointing antenna = 300
% Ae              = getOdtbxOptions(options, 'AtmAttenuation', 0.0 ); % dB
%                 % attenuation due to atmosphere (should be negative) [dB]      
% xmit_ant_mask   = getOdtbxOptions(options, 'TransAntMask', pi );  % in rad
%                 %  The actual mask used is the lesser of this mask and the limit of the defined pattern
%                 %  Note:  mask = 70 deg includes entire defined pattern
%                 %         mask = 42 deg includes only main and first side lobes
%                 %         mask = 26 deg includes only main lobe
% Nf              = getOdtbxOptions(options, 'ReceiverNoise', -3 );  % dB, Noise figure of receiver/LNA
% L               = getOdtbxOptions(options, 'RecConversionLoss', -1.5 );  % dB
% 				% Receiver implementation, A/D conversion losses [dB]
% 				%   Novatel: L = -4.0 	
% 				%   Plessey: L = -1.5		
% As              = getOdtbxOptions(options, 'SystemLoss', 0 ); % dB, System losses, in front of LNA
CN0_acq         = getOdtbxOptions(options, 'RecAcqThresh', 32 ); % dB-Hz, Receiver acquisition threshold
CN0_lim         = getOdtbxOptions(options, 'RecTrackThresh', CN0_acq ); % dB-Hz, Receiver tracking threshold


% the physical parameter results:
% epoch = out.epoch;       % epoch of first time [1x1]
TX_az = out.TX_az*d2r;   % The transmitter azimuth angle (rad) [nn x GPS_SIZE]
TX_link.alpha = out.TX_el*d2r; % The transmitter elevation angle (rad) [nn x GPS_SIZE]
RX_link.alpha = out.RX_el*d2r; % The receiver elevation angle (rad)

% Note, the receiver angles from "out" are below, in the ANT loop
Hrange = out.range;      % [GPS_SIZE x nn]
% Hrrate = out.rrate;      % [GPS_SIZE x nn]
rgps_mag = out.rgps_mag; % [nn x GPS_SIZE]
health = out.health;     % the health indicator, [nn x GPS_SIZE]
% dtsv = out.dtsv;          %individual satellite clock bias which is added to the epoch time to reflect the GPS time of measurement (using AFO and AF1)

% Set alpha_t for non-existent/unhealthy SVs to pi rad
TX_link.alpha(~health) = pi;

% Set alpha_r for non-existent/unhealthy SVs to 180 deg
RX_link.alpha(~health) = pi;

% Compute angle subtended by Earth and Earth mask angles for each SV
r_mask          = EARTH_RADIUS + AtmMask;	% Atmosphere mask radius (km)
denom=rgps_mag;
denom(denom==0)=NaN;
gamma = asin(EARTH_RADIUS./denom);   % Angle subtended by Earth at SV (nn,GPS_SIZE)
gamma_mask = asin(r_mask./denom);    % Angle subtended by Earth plus altitude mask (nn,GPS_SIZE)
%  Set prns visible if not blocked by Earth
vis_earth = (TX_link.alpha > gamma) | (Hrange' <= rgps_mag.*cos(gamma));   % (nn,GPS_SIZE)
%  Set prns visible if not subject to atmosphere mask
vis_atm = (TX_link.alpha > gamma_mask) | (Hrange' <= rgps_mag.*cos(gamma_mask)); % (nn,GPS_SIZE)

Hvis_earth = vis_earth';
Hvis_atm = vis_atm';

% Set receiver antenna loop number
loop = max([1,num_ant]);
GPS_SIZE = size(out.range,1);

%% Compute antenna specific data

% % set link budget params in structs
% RX_link.Nf = Nf;
% RX_link.L = L;
% % RX_link.freq = freq;
% RX_link.Ts = Ts;
% RX_link.As = As;
% RX_link.Ae = Ae;
% TX_link.P_sv = P_sv;
% 
% % Transmitter and receiver antenna masks
% TX_link.beta = xmit_ant_mask;  % User input from options
% RX_link.beta = rcv_ant_mask;   % User input from options

% ----------------------------------------
%  Antenna calculation loop
% ----------------------------------------
AntLB = cell(loop,1); % cell array of structs to hold link budget data
                      % for each antenna


for ANT=1:loop
    
    AntLB{ANT} = struct('Halpha_r',[],'Hvis_beta',[],'Hvis_CN0',[],...
        'HCN0',[],'HAd',[],'HAr',[],'HAP',[],'HRP',[],'HAt',[]);
    
    AntLB_raw = struct('CN0',[],'Ad',[],'Ar',[],'AP',[],'RP',[],'At',[]);
    
%     CN0 = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
%     Ar = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
%     At = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
%     Ad = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
%     AP = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
%     RP = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    
%     % The receiver elevation angle (rad)
%     alpha_r = out.RX_el(:,:,ANT)*d2r;
    
%     % Set alpha_r for non-existent/unhealthy SVs to 180 deg
%     RX_link.alpha(~health) = pi;
    
   
    % Determine if the pattern is elevation only (1-D) or azimuth and
    % elevation (2-D) and compute the receiver gain
    for j = 1:GPS_SIZE
        % Encapsulate RX and TX data
        % Originally set to be 1D receive, 1D transmit patterns
        RX_antenna = struct('pattern', RX_link.pattern{ANT}, ...
            'el', RX_link.alpha(:,j));
        TX_antenna = struct('pattern', TX_link.pattern, ...
            'el', TX_link.alpha(:,j));
        
        % Change dimensions on transmit patterns from 1D to 2D, if required
        if size(RX_link.pattern{ANT},2) > 2
            % 2D receive antenna
            RX_antenna.az = out.RX_az(:,j)*d2r;
        end
        if size(TX_link.pattern,2) > 2
            % 2D transmit antenna
            TX_antenna.az = TX_az(:,j);
        end
        
         % Compute gain/attenuation of receiving antenna pattern
        [AntLB_raw.CN0(:,j), AntLB_raw.Ar(:,j), AntLB_raw.At(:,j), ...
            AntLB_raw.Ad(:,j), AntLB_raw.AP(:,j), AntLB_raw.RP(:,j)] = ...
            linkbudget(Hrange(j,:)', RX_link, TX_link, RX_antenna, TX_antenna);
    end

    % Apply the receiver gain penalty for the user-defined mask angle
    AntLB_raw = gainpenalty_mask(RX_antenna, AntLB_raw, RX_link.alpha, RX_link.beta);
    
    % Apply the transmit gain penalty for the user-defined mask angle
    AntLB_raw = gainpenalty_mask(TX_antenna, AntLB_raw, TX_link.alpha, TX_link.beta);
    
    %------------------------------------------------------------------------------
    % EVALUATION OF GEOMETRIC CONSTRAINTS

    %  Set prns visible if los within antenna mask angles
    %  So far, alpha_t was computed assuming GPS antenna is nadir pointing
    vis_beta_t = (TX_link.alpha <= TX_link.beta);
    vis_beta = vis_beta_t & (RX_link.alpha <= RX_link.beta);    % (nn,GPS_SIZE)

    %  Set prns visible if CN0 is above acquisition/tracking threshold
    vis_CN0 = AntLB_raw.CN0 >= CN0_lim;                               % [nn,GPS_SIZE]

    %  OUTPUT PARAMETERS
    AntLB{ANT}.Halpha_r = RX_link.alpha';             % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_beta = vis_beta';             % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_CN0 = vis_CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HCN0 = AntLB_raw.CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAd = AntLB_raw.Ad';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAr = AntLB_raw.Ar';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAP = AntLB_raw.AP';             % [GPS_SIZE,nn]
    AntLB{ANT}.HRP = AntLB_raw.RP';             % [GPS_SIZE,nn]
    % Note: this variable can be used to pass out other values as well
    AntLB{ANT}.HAt = AntLB_raw.At';             % [GPS_SIZE,nn]

    % Mask undefined values for SV dependent parameters using (Health,
    %  Earth blockage, and xmit antennna masks)
    VIS_sv = vis_beta_t' & Hvis_earth & health';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAd(~VIS_sv) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HRP(~VIS_sv) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HAP(~VIS_sv) = -300;             % [GPS_SIZE,nn]

    % Mask undefined values for Antenna dependent parameters using
    % (Health, Earth blockage, and both antennna masks)
    VIS_ant = vis_beta' & Hvis_earth & health';             % [GPS_SIZE,nn]
    AntLB{ANT}.HCN0(~VIS_ant) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HAt(~VIS_ant) = -300;             % [GPS_SIZE,nn]
    
end % ANT loop


%% Combine results across multiple antennas
HVIS = visibility_constraints(AntLB, options, health, Hvis_earth, Hvis_atm);
end

