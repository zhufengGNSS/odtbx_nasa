function [CN0, Ar, At, Ad, AP, RP] = linkbudget(los_mag, RX_link, TX_link, RX_antenna, TX_antenna, CN0) 
%
% Calculates the GPS transmitter-receiver link budget.
%
% This calculates the full link budget between a GPS SV transmitter and a
% receiver system.  It includes the GPS SV transmitter power, the GPS
% antenna gain, the path loss, atmospheric attenuation, the receiver 
% antenna gain, the receiver system noise temperature, receiver system
% noise figures and losses, and received carrier to noise ratio.  It uses a
% set of system gain and loss constants with a set of distances and angles
% relative to the transmitter and receiver antennas to describe the link 
% budget geometry.
%
% Several methods of calculation can be used depending on the input
% specification:
% 1) Calculate CN0 (and other outputs) from both the transmit and receive
%    antenna patterns, the antenna angles, and the LOS magnitude.
% 2) Calculate Ar (and other outputs) from the transmit antenna pattern,
%    the transmit antenna angle(s), the LOS magnitude, and the CN0.
% 3) Calculate At (and other outputs) from the receive antenna pattern,
%    the receive antenna angle(s), the LOS magnitude, and the CN0.
%
% Either antenna pattern can be 1D, using only the antenna elevation
% angles, or 2D, using the antenna azimuth and elevation angles.  1D 
% patterns will ignore any supplied azimuth angle data.  Both types of 
% antenna angles must be supplied for a 2D antenna.  The transmit 
% antenna and receive antenna patterns can be different dimensions (i.e.
% one can be 1D while the other is 2D).
%
% If the receive pattern and receiver angles arguments are all empty then
% an ideal omnidirectional receive antenna with zero dB gain is assumed.
%
% Transmit or receive gains, At or Ar, are set to -100 dB if the angles are
% outside the prescribed pattern.  (This is not applied to the receive 
% gains with the omni receive option.)  All dependent outputs (CN0, RP, or 
% AP, as applicable) are biased as well.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   los_mag         double  Nx1     line of sight range magnitude (km)
%   RX_link         struct  1       struct of receiver link budget params
%                                   parameters with the following fields:
%       .Ts         double  1       System noise temp [K]
%       .Ae         double  1       attenuation due to atmosphere (should
%                                   be negative) [dB]
%       .Nf         double  1       dB, Noise figure of receiver/LNA
%       .L          double  1       Receiver implementation, A/D conversion losses [dB]
%       .As         double  1       dB, System losses, in front of LNA
%       .freq       double  1       carrier frequency [Hz]
%
%   TX_link         struct  1       struct of transmitter link budget params
%       .P_sv       double  1       spacecraft transmit power [dB]
%
%   RX_antenna      struct  1       struct of receiver antenna information
%                                   with the following fields:
%       .pattern double  AxG,[]  receive antenna gain pattern (deg & dB),
%                                   NOTE: can be 1D or 2D, optional: can be
%                                   empty if given CN0 with tx_pattern, or
%                                   can be empty if using an omni antenna
%       .el      double  Nx1     receiver antenna boresight
%                                   elevation angle (rad), optional: can be
%                                   empty if using an omni antenna
%       .az      double  Nx1,[]  receiver antenna boresight
%                                   azimuth angle (rad), optional: use with
%                                   a 2D rx_pattern
%
%   TX_antenna      struct  1       struct of receiver antenna information
%                                   with the following fields:
%       .tx_pattern double  AxG,[]  transmit antenna gain pattern (deg & dB),
%                                   NOTE: can be 1D or 2D, optional: can be
%                                   empty if given CN0 with rx_pattern
%       .tx_el      double  Nx1     transmitter antenna boresight
%                                   elevation angle (rad)
%       .tx_az      double  Nx1,[]  transmitter antenna boresight
%                                   azimuth angle (rad), optional: use with
%                                   a 2D tx_pattern
%   CN0             double  Nx1,[]  Signal carrier to noise ratio,
%                                   optional: only required when either
%                                   antenna pattern is not supplied
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   CN0             double  Nx1     Signal carrier to noise ratio, this is
%                                   either calculated or copied from the
%                                   given input
%   Ar              double  Nx1     receive antenna gain (dB)
%   At              double  Nx1     transmit anteanna gain (dB)
%   Ad              double  Nx1     Attenuation from R^2 losses (dB)
%   AP              double  Nx1     budget gain before receiver antenna
%                                   (dB)
%   RP              double  Nx1     budget gain before receiver amplifiers
%                                   and conversion (dB)
%
% The following are the link budget equations:
%   Ad = 20.*log10((C/freq)./(4*pi.*los_mag));
%   AP = P_sv + At + Ad + Ae;
%   RP = AP + Ar + As;
%   CN0 = RP - (10*log10(Ts)) + 228.6 + Nf + L;
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



% GPS SV (transmitter) antenna file definitions:
%  The format for 1D (elevation) antenna patttern files is as follows:
% 
%  % Header lines
%  % .
%  % .
%  % .
%  % Header lines
%  % angle [deg] gain [db]
%  0.0  gain1
%  .    .
%  .    .
%  angN gainN
% 
%  The comment character for antenna pattern files is the
%  "%" character and there can be as many header lines as
%  desired. The antenna pattern is specified as two columns
%  of numbers where the first column contains the angle
%  off the antenna boresight in units of degrees and the
%  second column contains the gain value in units of dB.
%  The angular resolution at which these data are
%  specified is up to the user, but the range is usually from 0 deg to 
%  90 deg or less.
%
%  The format for 2D antenna patttern files is as follows:
% 
%  % Header lines
%  % .
%  % .
%  % .
%  % Header lines
%  % One row with a placeholder scalar and then the applicable azimuth 
%  % angles [deg]
%  % Then on subsequent rows: the elevation angle [deg]
%  % deg or less and the gain value [db] for each applicable azimuth angle
%
%  -50  0.0     az2     ...     azM
%  0.0  gain1,1 gain1,2 ...     gain1,M
%  .    .
%  .    .
%  elN  gainN,1 gainN,2 ...     gainN,M
% 
%  The comment character for antenna pattern files is the
%  "%" character and there can be as many header lines as
%  desired. The antenna pattern is specified as a table
%  of numbers where the first column contains the angle
%  off the antenna boresight (elevation) in units of degrees and the
%  subsequent columns contains the gain value in units of dB.
%  The angular resolution at which these data are
%  specified is up to the user, but the elevation range is usually from 0 
%  deg to 90 deg or less.  The azimuth angles should range from 0 to < 360
%  providing complete azimuth coverage.

%% argument checks:
% data types that are provided:
have_CN0 = 0;
have_tx = 0;
have_rx = 0;
have_rx_omni = 0; % special rx case

if exist('CN0','var') && ~isempty(CN0) 
    have_CN0 = 1;
else
    CN0 = [];
end

if (exist('RX_antenna','var') && isfield(RX_antenna,'pattern') && ~isempty(RX_antenna.pattern))
    if (isfield(RX_antenna,'el') && ~isempty(RX_antenna.el))
        have_rx = 1; % we have a pattern and elevation angles
    end
    
    % I have questions about the logic in this section. The original
    % version had ~isempty(pattern). If we're looking to see if it's
    % empty (no inputs), wouldn't we want isempty(pattern)?
elseif (~exist('RX_antenna','var')) || (~isstruct(RX_antenna)) || ...
        (~isfield(RX_antenna,'pattern') || isempty(RX_antenna.pattern)) && ...
        (~isfield(RX_antenna,'el') || isempty(RX_antenna.el)) && ...
        (~isfield(RX_antenna,'az') || isempty(RX_antenna.az))
    have_rx = 1;
    have_rx_omni = 1; % no inputs means omni
end

if (exist('TX_antenna','var') && isstruct(TX_antenna) && isfield(TX_antenna,'tx_pattern') && ~isempty(TX_antenna.tx_pattern))
    if (isfield(TX_antenna,'tx_el') && ~isempty(TX_antenna.tx_el))
        have_tx = 1; % we have a pattern and elevation angles
    end
end

inputcheck = have_CN0 + have_rx + have_tx;
if (inputcheck) > 2
    error('Improper number of inputs, too many inputs specified.');
elseif (inputcheck) < 2
    error('Improper number of inputs, too few inputs specified.');
end

%% Calculate transmitter gain
if have_tx
    
    % Determine if the pattern is elevation only (1-D) or azimuth and
    % elevation (2-D)
    if(size(TX_antenna.tx_pattern,2) > 2)
        
        % arg check:
        if (~exist('TX_antenna','var') || ~isstruct(TX_antenna) || ~isfield(TX_antenna,'tx_az') || isempty(TX_antenna.tx_az))
            error('Presented 2D transmit model without azimuth angles - aborting.');
        end
        
        % use both transmitter azimuth and elevation to compute gain from a
        % 2D antenna model
        At = interp2(TX_antenna.tx_pattern(1,2:end)*pi/180,TX_antenna.tx_pattern(2:end,1)*pi/180,TX_antenna.tx_pattern(2:end,2:end),TX_antenna.tx_az,TX_antenna.tx_el,'spline');

        % The anglemask is the max el angle for which gain will be evaluated.  It is 
        % the minimum of the max defined angle for the antenna pattern.
        anglemask = (max(TX_antenna.tx_pattern(2:end,1)))*pi/180;  % [1,1]
    else
        % use only transmitter elevation angle to compute gain from a 1D
        % antenna model
        At = interp1((TX_antenna.tx_pattern(:,1))*pi/180,TX_antenna.tx_pattern(:,2),TX_antenna.tx_el,'spline');
        
        % The anglemask is the max angle for which gain will be evaluated.  It is 
        % the minimum of the max defined angle for the antenna pattern.
        anglemask = (max(TX_antenna.tx_pattern(:,1)))*pi/180;  % [1,1]
    end
    
    % Remove gains computed for invalid angles
    At(TX_antenna.tx_el>anglemask) = -100;

end

%% Calculate receiver gain
if have_rx
    if have_rx_omni
        Ar = 0;
    else
        % Determine if the pattern is elevation only (1-D) or azimuth and
        % elevation (2-D)
        if(size(RX_antenna.pattern,2) > 2)

            % arg check:
            if (~exist('RX_antenna.az','var') || isempty(RX_antenna.az))
                error('Presented 2D receiver model without azimuth angles - aborting.');
            end

            % use both receiver azimuth and elevation to compute gain from a
            % 2D antenna model
            Ar = interp2(RX_antenna.pattern(1,2:end)*pi/180,RX_antenna.pattern(2:end,1)*pi/180,RX_antenna.pattern(2:end,2:end),RX_antenna.az,RX_antenna.el,'spline');

        else
            % use only receiver elevation angle to compute gain from a 1D
            % antenna model
            Ar = interp1((RX_antenna.pattern(:,1))*pi/180,RX_antenna.pattern(:,2),RX_antenna.el,'spline');

        end

        % The anglemask is the max angle for which gain will be evaluated.  It
        % is  the minimum of the max defined angle for the antenna pattern.
        anglemask = (max(RX_antenna.pattern(:,1)))*pi/180;  % [1,1]

        % Remove gains computed for angles outside the pattern
        Ar(RX_antenna.el>anglemask) = -100;
    end
end

%% Link budget calculations
% the scalars that are independent of the unknonwn
scalars = TX_link.P_sv + RX_link.Ae + RX_link.As - (10*log10(RX_link.Ts)) + 228.6 + RX_link.Nf + RX_link.L;

% path loss
C = JATConstant('c') / 1000;  % km/s Speed of light
Ad = 20.*log10((C/RX_link.freq)./(4*pi.*los_mag));

if have_rx && have_tx
    CN0 = At + Ar + Ad + scalars;
elseif have_rx && have_CN0
    At = CN0 - Ar - Ad - scalars;
elseif have_tx && have_CN0
    Ar = CN0 - At - Ad - scalars;
else
    error('Logic error.'); % shouldn't hit this
end

AP = TX_link.P_sv + At + Ad + RX_link.Ae;
RP = AP + Ar + RX_link.As;

%% validation checks:
bias = CN0-Ar-At-Ad;
if (max(bias) - min(bias)) > 1e-6 % arbitrary dB tolerance check
    warning('Internal gpslinkbudget check failed.'); %#ok<WNTAG>
end
