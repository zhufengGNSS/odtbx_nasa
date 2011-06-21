function dest = copyGpsMetaData(source)
%
% Copies the metadata from a GPS data structure into a blank GPS data structure.
%
% INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   source          struct  1       A struct of populated GPS measurements
%                                   & metadata
%
% OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   dest            struct  1       A struct of GPS measurements
%                                   & metadata.  The data is blank but the
%                                   metadata is copied from source.
%
% GPS data description:
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   source/dest     struct  1       A struct of raw GPS measurements &
%                                   metadata
%       gps_data fields:
%   .RX_meta        struct  1       A struct to hold the metadata
%                                   associated with this gps_data, see
%                                   below for its fields
%   .GPS_PRN        double  1xM     The PRN numbers for the data in this
%                                   struct.  Each GPS_PRN matches once cell
%                                   in the PRN_data.  -1 means "no data"
%   .PRN_data       cell    1xM     Measurement data for a PRN
%       PRN_data fields:
%   .epoch          double  1xT     measurement time, using the ODTBX 
%                                   convention of datenum with convertTime 
%                                   GPS epoch
%   .raw_snr        double  1xT     raw signal strength measurement
%   .pseudorange    double  1xT     raw pseudorange measurement
%   .doppler        double  1xT     raw doppler measurement
%   .phase          double  1xT     raw phase measurement
%
% Note, M is the number of different PRNs in this dataset with
% measurements.  T is the number of measurement times, unique to each PRN.
%
% Fields of the gps_data.RX_meta:
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   .RX_ID          double  1       Receiver identifier (optional),
%                                   defaults to -1 which means not
%                                   specified
%   .meas_file      char    1xN     full file path and name of measurement
%                                   data (optional).
%   .obs_metadata   varies  1       Metadata from the measurement source,
%                                   type and content can vary by source,
%                                   (optional).
%
% Special handling for GPS times:
% Note that GPS time inputs and outputs are based on the GPS epoch of
% Jan 6, 1980 and not the standard datenum epoch.  The following correction
% is used:
% datenum('Jan 6 1980') = 723186 (exact).
%
% See also: convertTime, makeGpsData
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

dest = makeGpsData; % blank struct

% copy over the metadata
dest.RX_meta.RX_ID = source.RX_meta.RX_ID;
dest.RX_meta.meas_file = source.RX_meta.meas_file;

if iscell(source.RX_meta.obs_metadata)
    
    dest.RX_meta.obs_metadata = cell(size(source.RX_meta.obs_metadata));
    
    % copy the cell array contents
    for i = 1:length(source.RX_meta.obs_metadata)
        dest.RX_meta.obs_metadata{i} = source.RX_meta.obs_metadata{i};
    end
else
    % not a cell? straight copy
    dest.RX_meta.obs_metadata = source.RX_meta.obs_metadata;
end
