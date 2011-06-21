function phys_param = makePpData(gps_data)
%
% Creates a blank data structure for holding GPS physical parameters.
%
% Generally this function is used by gps_phys_params or FilterPpData and
% not by an end user.  If the gps_data is provided the RX_meta metadata is
% copied into the output phys_param struct.  If no argument is provided 
% then the phys_param metadata is left blank/default.
%
% INPUTS:
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   gps_data        struct  1       (Optional) gps_meas struct from
%                                   makeGpsData
%
% Returns:
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   phys_param      struct  1       Struct of GPS physical parameter data.
%           each phys_param fields:
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
% See also: makeGpsData, gps_phys_params, FilterPpData, makeGpsTXID
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

phys_param = struct('epoch',[],'TX_az',[],'TX_el',[],'RX_az',[],'RX_el',[],...
    'range',[],'range_rate',[],'GPS_yaw',[],'meta',[]);
phys_param.meta = struct('RX_meta',[],'RX_state_source',[],...
    'TX_ID',[],'TX_state_source',[],'gen_date',[]);

if nargin > 0
    % copy the RX metadata
    if ~isstruct(gps_data) || ~isfield(gps_data,'RX_meta')
        error('makePpData given invalid GPS data argument');
    end
    if ~isfield(gps_data.RX_meta,'RX_ID') || ...
            ~isfield(gps_data.RX_meta,'meas_file') || ...
            ~isfield(gps_data.RX_meta,'obs_metadata')
        error('makePpData given GPS data argument that has invalid metadata');
    end
    phys_param.meta.RX_meta.RX_ID = gps_data.RX_meta.RX_ID;
    phys_param.meta.RX_meta.meas_file = gps_data.RX_meta.meas_file;
    phys_param.meta.RX_meta.obs_metadata = [];
    for i = length(gps_data.RX_meta.obs_metadata)
        phys_param.meta.RX_meta.obs_metadata{i} = gps_data.RX_meta.obs_metadata{i};
    end
else
    phys_param.meta.RX_meta = struct('RX_ID',-1,'meas_file',[],'obs_metadata',[]);
end
