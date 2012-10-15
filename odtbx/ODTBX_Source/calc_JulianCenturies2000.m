function [T_JC] = calc_JulianCenturies2000(JD)
%
% Compute the time in Julian Centuries since the ephemeris epoch of 2000
%
% USAGE:   [T_JC] = calc_JulianCenturies2000(JD)
%
% PURPOSE: Compute the time in Julian Centuries since the ephemeris
%          epoch of 2000.
%
% INPUTS:  JD   - Time in Julian Date format.
%
% OUTPUTS: T_JC - Time in Julian Centuries since the ephemeris
%                 epoch of 2000.
%
% References:
%
% [1] http://en.wikipedia.org/wiki/Obliquity_of_the_ecliptic#Values
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

% Brent Wm. Barbee
% brent.barbee@emergentspace.com
%
% September 2009

% Begin function.

    % Compute the time in Julian Centuries since the ephemeris
    % epoch of 2000 (JD = 2451545.0).
    T_JC = (JD - 2451545.0)/36525;

% End of function.
end
