function [ob] = calc_obliquity(T_JC)
%
% Compute Earth obliquity of the ecliptic
%
% PURPOSE: Compute the angle known as the obliquity of the ecliptic,
%          which is the angle between the Earth's axis of rotation
%          the normal to the ecliptic plane.
%
% USAGE:   [ob] = calc_obliquity(T_JC)
%
% INPUTS:  T_JC - Time, measured in Julian Centuries, since the
%                 ephemeris epoch of 2000.
%
% OUTPUTS: ob   - Obliquity of the ecliptic, measured in radians.
%
% REFERENCES:
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

    % Compute the obliquity of the ecliptic in units of radians.
    ob = 0.409092802283074 - (0.000226966106587847*T_JC) - ((2.8623399732707e-09)*(T_JC^2)) + ((8.79645943005142e-09)*(T_JC^3));

% End of function.
end
