function [AntLB] = antennaCalc()
    AntLB{ANT} = struct('Halpha_r',[],'Hvis_beta',[],'Hvis_CN0',[],...
    'HCN0',[],'HAd',[],'HAr',[],'HAP',[],'HRP',[],'HAt',[]);
    
    CN0 = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    Ar = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    At = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    Ad = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    AP = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    RP = zeros(nn,GPS_SIZE);   % [nn,GPS_SIZE]
    
    % The receiver elevation angle (rad)
    alpha_r = out.RX_el(:,:,ANT)*d2r;
    
    % Set alpha_r for non-existent/unhealthy SVs to 180 deg
    alpha_r(~health) = pi;

    % Compute gain/attenuation of receiving antenna pattern

    % Determine if the pattern is elevation only (1-D) or azimuth and
    % elevation (2-D) and compute the receiver gain
    if rec_pattern_dim == 2 
        % 2D receive antenna

        for j = 1:GPS_SIZE
            if xmit_pattern_dim == 1 
                % 1D transmitter
                [CN0(:,j), Ar(:,j), At(:,j), Ad(:,j), AP(:,j), RP(:,j)] = gpslinkbudget(Hrange(j,:)', ...
                    RX_link, TX_link, ...
                    RXpattern{ANT}, alpha_r(:,j), out.RX_az(:,j)*d2r, ...
                    TXpattern, alpha_t(:,j), []);
            else
                % 2D transmitter
                [CN0(:,j), Ar(:,j), At(:,j), Ad(:,j), AP(:,j), RP(:,j)] = gpslinkbudget(Hrange(j,:)', ...
                    RX_link, TX_link, ...
                    RXpattern{ANT}, alpha_r(:,j), out.RX_az(:,j)*d2r, ...
                    TXpattern, alpha_t(:,j), TX_az(:,j));
            end
        end

    else
        % 1D receive antenna

        for j = 1:GPS_SIZE
            if xmit_pattern_dim == 1 
                % 1D transmitter
                [CN0(:,j), Ar(:,j), At(:,j), Ad(:,j), AP(:,j), RP(:,j)] = gpslinkbudget(Hrange(j,:)', ...
                    RX_link, TX_link, ...
                    RXpattern{ANT}, alpha_r(:,j), [], ...
                    TXpattern, alpha_t(:,j), []);
            else
                % 2D transmitter
                [CN0(:,j), Ar(:,j), At(:,j), Ad(:,j), AP(:,j), RP(:,j)] = gpslinkbudget(Hrange(j,:)', ...
                    RX_link, TX_link, ...
                    RXpattern{ANT}, alpha_r(:,j), [], ...
                    TXpattern, alpha_t(:,j), TX_az(:,j));
            end
        end

    end

    % Apply the receiver gain penalty for the user-defined mask angle
    if beta_r < (max(RXpattern{ANT}(:,1)))*pi/180
        % Check and apply additional mask angle penalty to:
        % Ar, RP, CN0
        % Note, the gpslinkbudget already applies a penalty for those 
        % angles outside the pattern so don't doubly-apply a penalty.
        Ar_mask_ind = (alpha_r > beta_r) & (Ar ~= -100.0);
        if sum(sum((Ar_mask_ind))) > 0
            Ar(Ar_mask_ind) = -100; % change Ar
            % alter dependent values:
            RP(Ar_mask_ind) = RP(Ar_mask_ind) - 100;
            CN0(Ar_mask_ind) = CN0(Ar_mask_ind) - 100;
        end
    end
    
    % Apply the transmit gain penalty for the user-defined
    % mask angle
    if beta_t < (max(TXpattern(:,1)))*pi/180
        % Check and apply additional mask angle penalty to: 
        % At, AP, RP, CN0
        % Note, the gpslinkbudget already applies a penalty for those 
        % angles outside the pattern so don't doubly-apply a penalty.
        At_mask_ind = (tx_angle_vector > beta_t)  & (At ~= -100.0);
        if sum(sum((At_mask_ind))) > 0
            At(At_mask_ind) = -100; % change At
            % alter dependent values:
            AP(At_mask_ind) = AP(At_mask_ind) - 100;
            RP(At_mask_ind) = RP(At_mask_ind) - 100;
            CN0(At_mask_ind) = CN0(At_mask_ind) - 100;
        end
    end
    
    %------------------------------------------------------------------------------
    % EVALUATION OF GEOMETRIC CONSTRAINTS

    %  Set prns visible if los within antenna mask angles
    %  So far, alpha_t was computed assuming GPS antenna is nadir pointing
    vis_beta_t = (alpha_t <= beta_t);
    vis_beta = vis_beta_t & (alpha_r <= beta_r);    % (nn,GPS_SIZE)

    %  Set prns visible if CN0 is above acquisition/tracking threshold
    vis_CN0 = CN0 >= CN0_lim;                               % [nn,GPS_SIZE]

    %  OUTPUT PARAMETERS
    AntLB{ANT}.Halpha_r = alpha_r';             % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_beta = vis_beta';             % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_CN0 = vis_CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HCN0 = CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAd = Ad';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAr = Ar';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAP = AP';             % [GPS_SIZE,nn]
    AntLB{ANT}.HRP = RP';             % [GPS_SIZE,nn]
    % Note: this variable can be used to pass out other values as well
    AntLB{ANT}.HAt = At';             % [GPS_SIZE,nn]

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
    
    
end