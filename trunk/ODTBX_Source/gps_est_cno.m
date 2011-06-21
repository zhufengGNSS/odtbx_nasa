function [CN0, fil_phys_param] = gps_est_cno(phys_param, RX_link, TX_link, rx_pattern, tx_pattern, filters)
%
% Calculates the estimated carrier to noise ratio from physical parameters and link budget data.
%
% Description: gps_est_cno estimates the carrier to noise ratio based on
% the provided physical parameters (including azimuth and elvation angles,
% and the antenna-to-antenna line of sight range), the 2-dimensional
% antenna patterns the transmit antenna and receive antenna with the
% provided remaining parts of the link budget.
%
% The physical parameters can optionally be filetred before processing using
% one or more provided FilterPpData instances.  If multiple filters are
% supplied then they are applied in the given order.  The filtered phys_param
% results are returned as an optional output argument.  This argument can 
% be left empty.
%
% An empty CN0 will be returned if no data can be processed.
%
% For a description of the link budget equations, see gpslinkbudget.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
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
%       .RX_meta	struct	1       Measurement metadata from the gps_meas 
%                                   input, see makeGpsData
%       .RX_state_source	string  1   see inputs
%       .TX_ID      struct	1	GPS SV/PRN identifier
%       .TX_state_source	string  1   The GPS almanac file used for
%                                   processing.
%       .gen_date double	1       MATLAB datenum date/time of when this
%                                   data was processed
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
%   rx_pattern      double  AxG     receive antenna gain pattern (deg & dB)
%   tx_pattern      double  AxG     transmit antenna gain pattern (deg &
%                                   dB)
%   filters         cell    1xF     (optional) A cell array of one or more
%                                   FilterPpData filters to be applied to
%                                   the phys_param data before processing.
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   CN0             double  2xN     Empty, or 
%                                   (1,:) time, UTC epoch of the measurement,
%                                   as a MATLAB serial date, see ConvertTime
%                                   (2,:) Signal carrier to noise ratio [dB]
%   fil_phys_param  struct  1       GPS physical parameter struct
%                                   containing the filtered input
%                                   argument phys_param.
%                                   
% See also: gpslinkbudget, gps_phys_params, FilterPpData
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
if isempty(phys_param) || isempty(RX_link) || ...
        isempty(TX_link) || isempty(rx_pattern) || isempty(tx_pattern)
    error('gps_est_cno: One or more required arguments were empty, aborting.');
end

if size(rx_pattern,1) < 3 || size(rx_pattern,2) < 3
    error('gps_est_cno: Not provided a 2D receive pattern, aborting.');
end

if size(tx_pattern,1) < 3 || size(tx_pattern,2) < 3
    error('gps_est_cno: Not provided a 2D transmit pattern, aborting.');
end

if ~exist('filters','var') || isempty(filters)
    dofilter = 0;
else
    dofilter = 1;
end

if dofilter  && ~iscell(filters)
    if isa(filters,'FilterPpData')
        % place it inside a cell array for below
        tmp = filters;
        filters = cell(1);
        filters{1} = tmp;
    else
        error('gps_est_cno: Invalid filter argument type (must be a cell array), aborting.');
    end
end

%% filter the phys_param
fil_phys_param = phys_param;
if dofilter
    for i = 1:length(filters)
        if ~isa(filters{i},'FilterPpData')
            error('gps_est_cno: Invalid filter argument type (must be a cell array), aborting.');
        end
        f = filters{i};
        fil_phys_param = f.filter(phys_param);
    end
end

%% post-filtering checks
% check for any remaining data, if none then exit early
if isempty(fil_phys_param.epoch)
    CN0 = [];
    return;
end

%% Calculate CNO from the physical parameters and link budget parameters
CN0 = zeros(2,size(fil_phys_param.epoch,2));
CN0(1,:) = fil_phys_param.epoch;
CN0(2,:) = gpslinkbudget(fil_phys_param.range', RX_link, TX_link, rx_pattern, fil_phys_param.RX_el', fil_phys_param.RX_az', tx_pattern, fil_phys_param.TX_el', fil_phys_param.TX_az');

