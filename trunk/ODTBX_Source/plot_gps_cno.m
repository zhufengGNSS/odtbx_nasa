function [h] = plot_gps_cno(gps_data, PRNs)
% Plots the raw GPS signal to noise ratio vs time.
%
% Plots the GPS signal to noise ratio for each GPS PRN specified in the
% argument.  Multiple PRNs are plotted on subplots, but four or less PRN 
% subplots are recommended for readability.  If there is no data for that 
% PRN then there will be an empty subplot.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   gps_data        struct  1       A struct of raw GPS measurements
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
%   .raw_SNR        double  1xT     raw signal strength measurement
%   .pseudorange    double  1xT     raw pseudorange measurement
%   .doppler        double  1xT     raw doppler measurement
%   .phase          double  1xT     raw phase measurement
%
%   PRNs            double  1xN     An array of PRN numbers to
%                                   plot/subplot.
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   h               handle  1       figure handle
%
% See also: makeGpsData, rinexo2gpsdata, FilterGpsData
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

sp = length(PRNs);

h = figure;

% PRNs in the data
[prns, datind, argind] = intersect(gps_data.GPS_PRN, PRNs);
missingprns = setdiff(PRNs, gps_data.GPS_PRN);

% find the earliest common time
epochs = zeros(1,length(datind));
for i = 1:length(datind)
    epochs(i) = gps_data.PRN_data{datind(i)}.epoch(1); % note, GPS datenum
end
baseepoch = min(epochs);

yl = cell(1,2);

for i = 1:sp
    if sp > 1
        subplot(sp,1,i)
    end
    if i <= length(prns)
        plot((gps_data.PRN_data{datind(i)}.epoch-baseepoch)*24,gps_data.PRN_data{datind(i)}.raw_SNR,'.');
    end
    yl{1} = 'C/N_0';
    if i <= length(prns)
        yl{2} = ['PRN:' num2str(gps_data.GPS_PRN(datind(i)))];
    else
        yl{2} = ['PRN:' num2str(missingprns(i-length(prns))) ];
    end
    ylabel(yl);
    
    if i == 1 && ~isempty(prns)
        ts = cell(2,1);
        if sp < 2
            ts{1} = sprintf('Raw C/N_0 for GPS SV PRN(s): %s',num2str(sort(PRNs)));
        else
            ts{1} = sprintf('Raw C/N_0 for GPS SV PRN(s): [%s]',num2str(sort(PRNs)));
        end
        % Note, correct for the GPS data's datenum bias, see convertTime
        ts{2} = ['Start date: ' datestr(baseepoch+723186,31) ' GPS'];
        title(ts);
    end
    
    if i == sp
        xlabel('Time (hours)');
    end
end