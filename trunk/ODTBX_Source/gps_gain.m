function gain = gps_gain(gps_meas, phys_param, RX_link, TX_link, pattern, ant)
%
% Calculates antenna gains from GPS physical parameters, CN0, and link budget data.
%
% Description: gps_gain calculates antenna gain as a function of boresight 
% azimuth and elvation angles and time for either the transmit antenna or 
% the receive antenna when provided with the remaining parts of the
% link budget.  The provided antenna pattern and ant flag specify which
% antenna gain is to be calculated.  The measured carrier to noise ratio 
% from the gps_meas struct is used along with the RX_link and TX_link structs
% and pattern to complete the link budget.
%
% The resulting gain data is only calculated if the provided data passes
% the following checks.  An empty gain array will be returned if no data
% can be processed.  gps_gain checks the phys_param.meta.TX_ID.GPS_PRN 
% against the values in gps_meas.GPS_PRN to be sure they are for the same 
% PRN value.  If there is no match then no calculations occur because the
% data is not for the same transmitter.  This tool checks the epoch 
% timestamps between the phys_param.epoch, and associated 
% gps_meas.PRN_data{}.epoch.  Only epochs that exactly match are processed.
% For a description of the link budget equations, see gpslinkbudget.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   gps_meas        struct  1       gps measurements struct containing receiver
%                                   data, see makeGpsData
%           gps_meas fields:
%       .RX_meta    struct  1       A struct to hold the metadata
%                                   associated with this gps_data, see
%                                   below for its fields
%       .GPS_PRN    double  1xM     The PRN numbers for the data in this
%                                   struct.  Each GPS_PRN matches once cell
%                                   in the PRN_data.  -1 means "no data"
%       .PRN_data   cell    1xM     Measurement data for a PRN in a struct,
%           PRN_data fields:
%       .epoch      double  1xT     measurement time, using the ODTBX 
%                                   convention of datenum with convertTime 
%                                   GPS epoch
%       .raw_snr    double  1xT     raw signal strength measurement
%       (other fields in .PRN_data{} are ignored)
%
%   phys_param      struct  1       GPS physical parameter struct
%                                   containing data for one GPS
%                                   PRN/SV-receiver combination.
%           phys_param fields:
%       .epoch      double	1xZ	    time, UTC epoch of the measurement, as
%                                   a MATLAB serial date, see ConvertTime
%       .TX_az      double	1xZ     transmitter antenna azimuth angle (rad)
%       .TX_el      double	1xZ     transmitter antenna boresight elevation
%                                   angle (rad)
%       .RX_az      double	1xZ     receiver antenna azimuth angle (rad)
%       .RX_el      double	1xZ     receiver antenna boresight elevation angle
%                                   (rad)
%       .range      double	1xZ     absolute value of the transmitter-receiver
%                                   range (km)
%       .range_rate	double	1xZ the transmitter-receiver relative
%                                   velocity projected along the LOS (km/s)
%       .GPS_yaw	double	1xZ     the GPS yaw model yaw angle of the
%                                   satellite (deg)
%       .meta       struct	1       metadata for this phys_param
%           meta fields:
%       .TX_ID      struct	1	    GPS SV/PRN identifier
%           (other fields in phys_param.meta are ignored)
%
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
%                                   fields:
%       .P_sv       double  1       spacecraft transmit power [dB]
%
%   pattern         double  AxG     provided antenna gain pattern (deg &
%                                   dB) for receive or transmit antenna
%   ant             double  1       provided antenna pattern is: 1=receive 
%                                   pattern to calculate transmit gain, 
%                                   2=transmitpattern to calculate receive 
%                                   gain
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   gain            double  4xN     Empty, or:
%                                   (1,:) time, UTC epoch of the measurement,
%                                       as a MATLAB serial date, see
%                                       ConvertTime
%                                   (2,:) the calculated gain [dB]
%                                   (3,:) the antenna azimuth angle (rad)
%                                   (4,:) the antenna elevation angle (rad)
%
% See also: makeGpsData, gpslinkbudget, gps_phys_params, makeGpsTXID
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

%% check inputs
if isempty(gps_meas) || isempty(phys_param) || isempty(RX_link) || ...
        isempty(TX_link) || isempty(pattern) || isempty(ant)
    error('gps_gain: One or more required arguments were empty, aborting.');
end

if size(pattern,1) < 3 || size(pattern,2) < 3
    error('gps_gain: Not provided a 2D pattern, aborting.');
end

if (ant < 1) || (ant > 2)
    error('gps_gain: invalid antenna enum value (%d), aborting.',ant);
end

%% processing
gain = []; % for early exits

% get the PRN from the TX_ID
PRN = phys_param.meta.TX_ID.GPS_PRN;
prnind = find(gps_meas.GPS_PRN == PRN); % index of data in gps_meas.PRN_data

% cross-check args
if length(prnind) < 1
    % warn here because this could occur with upstream filtering
    warning('gps_gain: GPS PRN mismatch between the GPS measurement data and physical params (PRN=%d), (no processing).\n',PRN); %#ok<WNTAG>
    return;
elseif length(prnind) > 1
    % this really shouldn't happen - the must be an error in how gps_meas was constructed.
    error('gps_gain: Too many matching GPS PRNs in a (badly constructed) GPS measurement data struct (PRN=%d).',PRN);
end

% find the common times between the gps_meas and the phys_param
% ipp = the data to use from phys_param.*(ipp)
% igm = the data to use from gps_meas.PRN_data{prnind}.*(igm)
% Use a 0.1 sec tolerance for the conversion
ipp = compValsTol(phys_param.epoch, convertTime('UTC','GPS',gps_meas.PRN_data{prnind}.epoch), 0.1/86400);
igm = compValsTol(gps_meas.PRN_data{prnind}.epoch, convertTime('GPS','UTC',phys_param.epoch), 0.1/86400);

% check: no common times
if isempty(ipp) || isempty(igm)
    % warn here because this could occur with upstream filtering
    warning('gps_gain: No matching epoch timestamps found between the provided gps_meas and phys_param.'); %#ok<WNTAG>
    return;
end

CN0 = gps_meas.PRN_data{prnind}.raw_SNR(igm)';
range = phys_param.range';

gain = zeros(4,length(ipp));
gain(1,:) = phys_param.epoch(ipp);

%% Calculate CNO from the physical parameters and link budget parameters
if ant == 1
    % use receive pattern to calculate transmit gains
    [~, ~, gain(2,:)] = gpslinkbudget(range, RX_link, TX_link, pattern, phys_param.RX_el(ipp)', phys_param.RX_az(ipp)', [], phys_param.TX_el(ipp)', phys_param.TX_az(ipp)', CN0);
    gain(3,:) = phys_param.TX_az(ipp);
    gain(4,:) = phys_param.TX_el(ipp);
elseif ant == 2
    % use transmit pattern to calculate receive gains
    [~, gain(2,:)] = gpslinkbudget(range, RX_link, TX_link, [], phys_param.RX_el(ipp)', phys_param.RX_az(ipp)', pattern, phys_param.TX_el(ipp)', phys_param.TX_az(ipp)', CN0);
    gain(3,:) = phys_param.RX_az(ipp);
    gain(4,:) = phys_param.RX_el(ipp);
end

