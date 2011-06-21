function TX_ID = makeGpsTXID(GPS_ID, GPS_PRN, block_type)
%
% Creates a blank data structure for describing a GPS transmitter.
%
% A GPS Transmitter ID is a way to uniquely identify a particular GPS
% transmitter.  It relates the GPS vehicle block type and the current 
% PRN assignment for that transmitter to a user-supplied GPS ID.  The PRN
% assignment for a particular GPS transmitter may change over time but the
% assigned GPS ID and block type should not change.
%
% Returns:
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   TX_ID           struct  1       A struct to hold raw GPS measurements
%       its fields:
%   .GPS_ID         double  1       A user-supplied GPS satellite vehicle 
%                                   identifier (note, this is not the GPS PRN)
%                                   (defaulted to -1)
%   .GPS_PRN        double  1       The assigned GPS PRN number for this
%                                   transmitter (defaulted to -1)
%   .block_type     double  1       The GPS Satellite Block type:
%                                               1=II/IIA (default)
%                                               2=IIR
%                                               3=IIR-M
%                                               4=IIF
%                                               5=III
%                                   (This numerical value should be
%                                   consistent with the gpsmeas GPSBlock 
%                                   option.)
%
% See also: gpsmeas
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

if isempty(GPS_ID)
    GPS_ID = -1;
end

if isempty(GPS_PRN)
    GPS_PRN = -1;
end

if isempty(block_type)
    block_type = 1;
end

TX_ID = struct('GPS_ID',GPS_ID,'GPS_PRN',GPS_PRN,'block_type',block_type);

end
