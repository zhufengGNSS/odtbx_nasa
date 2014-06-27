%
% gpsmeas_test
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

% Brent Wm. Barbee
%
% This function uses code from the following two scripts by Allen Brown:
%
% checkcalc_gps_posvel.m
% plot_checkcalc_gps_posvel.m

% Begin function.
function [failed] = gpsmeas_test(truth_file, TOL, yuma_file, ant_point, ...
                                 ant_pat, dump_results)
%
% Runs a GPS satellite simultaion test case.  Will return whether the
% simulation results match the results in a truth file.
%
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%   INPUTS
%    truth_file     string  the name of the file containing what GPS
%                           simulation results should be.  This also
%                           contains simulation setup parameters.
%    TOL            1x1     tolerance when comparing simulation results and
%                           expected results before flagging a failure.
%                           When comparing ranges, TOL is interpreted as in
%                           meters.  When comparing range rate, interpreted
%                           as m/s.
%    yuma_file      string  the name of the GPS almanac file used to
%                           determine GPS satellite positions
%    ant_point  1 x num_ant Specify attitude profile for each antenna
%                           (1) zenith pointing or (-1) nadir pointing wrt
%                               geocentric LVLH
%                           (2) parallel or (-2) antiparallel to the
%                               Earth-Sun vector
%                           (3) ecliptic north or south (-3)
%                           (4) fixed with respect to apogee zenith, or
%                               (-4) nadir vector apogee selected as the
%                               point of highest altitude (ephemeris must
%                               include apogee)
%                           (5) body fore and (-5) aft directions relative
%                               to geocentric LVLH
%                           (6) body port and (-6) starboard directions
%                               relative to geocentric LVLH  
%    ant_pat    1 x num_ant Specify antenna pattern for each antenna, e.g.                                              sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
%                           omni.txt - zero dB gain,  
%                                      180 degree half beamwidth
%                           trimblepatch_ant.txt - hemi antenna, 
%                                      4.5 dB gain, 90 deg half beamwidth
%    dump_results   1x1     Optional.  If 1, will store the simulation
%                           results in a .mat file which can later be used
%                           as a truth file for the test case.  The created
%                           file will have the same name as the truth file
%                           with "_new" appended.
%                           If not specified, will not dump results.
%   OUTPUTS
%    failed         1x1     0 if all test cases passed.  
%                           1 if one of the test cases failed.
    if (nargin < 6) 
        dump_results = 0;
    end
   
    % Default to 'pass'.
    failed = 0;

    % Load the truth data file.
    truth = load(truth_file);

    %% set up time
    % gpstools time is GPS seconds
    % gpsmeas time is seconds from epoch, with an epoch in UTC
    epoch = datenum(truth.sim_start_utc);
    % create the time vector using the same time intervals, referenced from the
    % first time
    t = truth.time - truth.time(1); 

    %% set up satellite state
    % gpstools state is in m (ecef for leo case) 
    % gpsmeas state is in km, default to eci unless a rotation function is
    % supplied
    xecef(1:3,:) = truth.sat_pos/1000;
    xecef(4:6,:) = truth.sat_vel/1000;

    %% Convert true position and velocity from ECEF to ECI
    ecef2eci  = jatDCM('ecef2eci', (t/86400)+epoch); 
    w         = [0; 0; JATConstant('wEarth')];
    % change name, added tmp_
    tmp_sat_pos = zeros(3, length(t));
    tmp_sat_vel = zeros(3, length(t));
    for n=1:length(t)
        tmp_sat_pos(1:3,n) = ecef2eci(:,:,n)*xecef(1:3,n);
        tmp_sat_vel(1:3,n) = ecef2eci(:,:,n)*(cross(w,xecef(1:3,n))+xecef(4:6,n));
    end
    x = zeros(10,length(t));
    x(1:3,:) = tmp_sat_pos; % km, eci
    x(4:6,:) = tmp_sat_vel; % km/s, eci
    % (end of state conversion)

    %% Compute attitude quaternion assuming zenith pointing
    % Note that the antenna z-axis should be pointing zenith.  The body to 
    % antenna rotation is done through the measOptions belows with the
    % parameter AntennaOrientation.
    % The following is used only if the receiver antenna is 2-D.
    R = dcm('ric', x(1:3,:),x(4:6,:)); % Rotation from ECI to RIC
    qatt = dcm2q(R);
        
    % gpsmeas uses atmosphere height above earth for 'AtmosphereMask' below, 
    % gpstools uses a total distance from the earth's center.  This is how
    % gpsmeas defines the earth radius and adds it to 'AtmosphereMask'.
    EARTH_RADIUS = JATConstant('rEarth','WGS84') / 1000;  % km Equatorial radius of Earth

    measOptions = odtbxOptions('measurement');
    measOptions = setOdtbxOptions(measOptions, ...
                  'epoch', epoch, ...
                  'useRange', true, ...
                  'useRangeRate', true, ...
                  'useDoppler', false, ...
                  'YumaFile', yuma_file, ...
                  'AntennaPointing', ant_point, ... % override of truth.pointing_ref
                  'useIonosphere', false, ...
                  'useTroposphere', false, ...
                  'PrecnNutnExpire', 0.1, ... % days
                  'AntennaOrientation', dcm('ax2', pi/2));
                         
    %rSigma (use default, no analogue in gpsdef)
    %Rotation2ECI (use default, no analogue in gpsdef)

    link_budget.GPSBand           = 'L1';                             % see truth.freq
    link_budget.AntennaPattern    = eval(ant_pat);
    link_budget.RXAntennaMask     = truth.rcv_ant_mask;
    link_budget.AtmosphereMask    = truth.r_mask/1000 - EARTH_RADIUS; % km
    link_budget.NoiseTemp         = truth.Ts;                         % K
    link_budget.AtmAttenuation    = truth.Ae;                         % dB
    link_budget.TransPowerLevel   = truth.sv_power;                   % enum
    link_budget.TransPowerOffset  = truth.xmit_power_offset;          % truth units? gpsmeas: dB
    link_budget.GPSBlock          = truth.sv_block;                   % enum
    link_budget.TXAntennaMask     = truth.xmit_ant_mask;              % rad
    link_budget.ReceiverNoise     = truth.Nf;                         % dB
    link_budget.RecConversionLoss = truth.L;                          % dB
    link_budget.SystemLoss        = truth.As;                         % dB
    link_budget.LNAGain           = truth.Ga;                         % dB
    link_budget.CableLoss         = truth.Ac;                         % dB
    link_budget.RecAcqThresh      = truth.CN0_lim;                    % dB-Hz
    link_budget.RecTrackThresh    = truth.CN0_lim;                    % dB-Hz
    link_budget.DynamicTrackRange = truth.dyn_range;                  % dB
    link_budget.TX_AntennaPointing= -1; % 1 for zenith pointing, -1 for nadir pointing
    measOptions = setOdtbxOptions(measOptions, 'linkbudget', link_budget);

    
    % Call gpsmeas.
    % Note that even though the spacecraft attitude is specified as the
    % fourth argument below, it is used only if the receiver antenna
    % pattern is 2-D.  Otherwise, the receiver antenna pointing is defined
    % by the AntennaPointing parameter in the measOptions structure.
    [y,H,~,AntLB] = gpsmeas(t, x, measOptions, qatt); 

    % Check the results:

    numsats = size(truth.Hrange, 1);
    if(size(y, 1) == numsats)
        % gpsmeas y output doesn't contain range rate data
        dorrate = 0; 
        yrange = y;
    else
        % gpsmeas y output does contain range rate data
        dorrate = 1; 
        yrange = y(1:2:size(y,1), :) * 1000; %m
        yrrate = y(2:2:size(y,1), :) * 1000; %m/s
    end
    yvis = ~isnan(yrange);
    
    if (dump_results) 
        strend = length(truth_file);
        suffix = truth_file(strend-3:strend);
        if (strcmp(suffix,'.mat'))
            new_file_name = truth_file(1:strend-4);
        else
            new_file_name = truth_file;
        end
        new_file_name = [new_file_name '_new.mat'];
        results.Hrange = yrange;
        results.Hrrate = yrrate;
        results.vis = yvis; 
        results.H = H;
        results.AntLB = AntLB; %#ok<STRNU>
        save(new_file_name, '-struct', 'truth');
        save(new_file_name, '-struct', 'results', '-append'); 
        display(['Saved results to ' new_file_name]);
        return
    end

    dispString = '';
    for s=1:numsats

        % Compute differences.
        dR = yrange(s,:) - truth.Hrange(s,:); % m
        if(dorrate)
            dV = yrrate(s,:) - truth.Hrrate(s,:); % m/s
        end
    
        % Perform visibility checks.
        vis_mismatch = xor(truth.vis(s,:), yvis(s,:));
        mismatch_ind = find(vis_mismatch);

        e1 = max(abs(dR));
        if(dorrate)
            e2 = max(abs(dV));
        end
        e3 = length(mismatch_ind);

        % Print statistics to the screen.
        if(dorrate)
            statString = sprintf('\nSat: %d: max range err %f (m), max rangerate err %f (m/s), vis mismatch: %d',...
                s, e1, e2, e3);
        else
            statString = sprintf('\nSat: %d: max range err %f (m), vis mismatch: %d\n',...
                s, e1, e3);
        end

        if(dorrate)
            if((e1 > TOL) || (e2 > TOL) || (e3 ~= 0))
                failed = 1;
                statString = upper(statString);
            end
        else
            if((e1 > TOL) || (e3 ~= 0))
                failed = 1;
                statString = upper(statString);
            end
        end
        dispString = strcat(dispString, statString);
        
    end
    
    if (failed) 
        disp(dispString);
    end
    
    if max(max(max(abs(H-truth.H)))) > TOL
        failed = 1;
        disp('Failed: Measurement partials do not match within tolerance.')
    end
    
    num_ant = length(truth.AntLB);
    for ANT = 1:num_ant
        if max(max(abs(truth.AntLB{ANT}.HCN0 - AntLB{ANT}.HCN0))) > TOL
            failed = 1;
            disp(['Failed: CNO for antennna',num2str(ANT),' do not match']);
        end
    end
    
    % quick sanity check:
    % check the y result visibility against the antenna visibility:
    num_ant = length(AntLB);
    max_CN0 = zeros(size(yvis,1),size(yvis,2));
    ck_antvis = logical(max_CN0);
    
    for ANT = 1:num_ant
        ck_antvis = ck_antvis | (AntLB{ANT}.Hvis_CN0 & AntLB{ANT}.Hvis_beta);
        max_CN0 = max(max_CN0,AntLB{ANT}.HCN0);
    end
        
    % y should never be a value when no receive antenna sees a transmit SV
    % or when the CN0 is really low (assume less than zero)
    sanity_check_fail = yvis & (~ck_antvis) & (~(max_CN0 > 0));
    if any(any(sanity_check_fail))
            failed = 1;
            disp('Failed sanity check between measurement y output, antenna visibility, and max CN0.');
    end
    

% End of function.
    end

