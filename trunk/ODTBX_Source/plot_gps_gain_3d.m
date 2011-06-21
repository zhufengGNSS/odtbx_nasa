function [h] = plot_gps_gain_3d(ant, gain, PRN, pattern)
% Plots the gain from gps_gain and/or antenna gain pattern vs az and el angles (3D plots).
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   ant             double  1       provided antenna data is: 1=receive, 
%                                   2=transmit
%   gain            double  4xN     (Optional, see below)
%                                   (1,:) time, UTC epoch of the measurement,
%                                       as a MATLAB serial date, see
%                                       ConvertTime
%                                   (2,:) the calculated gain [dB]
%                                   (3,:) the antenna azimuth angle (rad)
%                                   (4,:) the antenna elevation angle (rad)
%   PRN             double  1       (Optional) The PRN number for this data
%   pattern         double  AxG     (Optional, see below) Co-plot the 2D 
%                                   antenna gain pattern (deg & dB)
%
% Note, either gain, pattern, or both must be supplied.
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   h               handle  1       figure handle
%
% See also: gps_gain
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
else
    antname = 'Transmit';
end

dosurf = 0;
if exist('pattern','var')
    dosurf = 1;
end

dodata = 0;
if exist('gain','var')
    dodata = 1;
end

if ~dosurf && ~dodata
    error('plot_gps_gain_3d: Either gain data or an antenna pattern must be supplied as arguments.');
end

h = figure;
if dodata
    plot3(gain(3,:)*r2d, gain(4,:)*r2d, gain(2,:),'.');
end
if dosurf
    hold on;
    surf(pattern(1,2:end),pattern(2:end,1),pattern(2:end,2:end),'LineStyle','none');
    hold off;
end
view(-113,34);

xlabel('Azimuth Angle (deg)');
ylabel('Elevation Angle (deg)');
zlabel('Gain (dB)');

ts = cell(2,1);
if exist('PRN','var')
    ts{1} = sprintf('%s Gain for GPS SV PRN %d',antname, PRN);
else
    ts{1} = [antname ' Gain'];
end
ts{2} = ['Data dates: ' datestr(gain(1,1),31) ' to ' datestr(gain(1,end),31) ' UTC'];
title(ts);
