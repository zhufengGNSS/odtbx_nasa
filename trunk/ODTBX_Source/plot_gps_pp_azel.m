function [h] = plot_gps_pp_azel(phys_param, ant)
% Plots antenna az and el angle physical parameters from gps_phys_params vs time.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   phys_param      struct  1       Struct of GPS physical parameters
%                                   containing data for one GPS PRN/SV
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
%   ant             double  1       provided antenna data is: 1=receive, 
%                                   2=transmit
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   h               handle  1       figure handle
%
% See also: gps_phys_params
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

r2d=180/pi; % rad to deg

if ant == 1
    antname = 'Receive';
    az = phys_param.RX_az*r2d;
    el = phys_param.RX_el*r2d;
else
    antname = 'Transmit';
    az = phys_param.TX_az*r2d;
    el = phys_param.TX_el*r2d;
end

h = figure;
subplot(2,1,1);
plot((phys_param.epoch-phys_param.epoch(1))*24,el,'.');
ylabel('Elevation (deg)');

ts = cell(2,1);
ts{1} = sprintf('%s Angles for GPS SV PRN %d',antname,phys_param.meta.TX_ID.GPS_PRN);
ts{2} = ['Start date: ' datestr(phys_param.epoch(1,1),31) ' UTC'];
title(ts);

subplot(2,1,2);
plot((phys_param.epoch-phys_param.epoch(1))*24,az,'.');
ylabel('Azimuth (deg)');

xlabel('Time (hours)');
