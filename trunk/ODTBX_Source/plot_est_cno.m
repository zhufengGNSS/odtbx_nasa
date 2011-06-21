function [h] = plot_est_cno(CN0, PRN)
% Plots the estimated Carrier to Noise ratio from gps_est_cno vs time.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   CN0             double  2xN     (1,:) time, UTC epoch of the measurement,
%                                   as a MATLAB serial date, see ConvertTime
%                                   (2,:) Signal carrier to noise ratio
%                                   [dB]
%   PRN             double  1       (Optional) The PRN number for this data
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   h               handle  1       figure handle
%
% See also: gps_est_cno
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

h = figure;
plot((CN0(1,:)-CN0(1,1))*24,CN0(2,:),'.');
ylabel('C/N_0');

xlabel('Time (hours)');

ts = cell(2,1);
if exist('PRN','var')
    ts{1} = sprintf('Estimated C/N_0 for GPS SV PRN %d',PRN);
else
    ts{1} = 'Estimated C/N_0';
end
ts{2} = ['Start date: ' datestr(CN0(1,1),31) ' UTC'];
title(ts);