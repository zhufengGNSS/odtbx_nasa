function x = LLA2ecef(lat, lon, alt)

% LLA2ECEF  Returns the ECEF position from latitude, longitude and altitude
%
%  x = LLA2ecef(lat, lon, alt) returns the ECEF positions of the input latitude,
%  longitude and altitude positions.
%
% INPUTS
%   VARIABLE    SIZE        DESCRIPTION (Optional/Default)
%     lat       (1xN)       Latitudes in degs
%     lon       (1xN)       Longitudes in degs
%     alt       (1xN)       Altitudes in km
%
% OUTPUTS
%     x         (3xN)       ECEF Positions in km
%
% VALIDATION/REGRESSION TEST
%
%   These tests were moved to LLA2ecef_test.m to conform to the new
%   regression testing framework.
%
%   keywords: Coordinate Transformations, JAT Adapter
%   See also ECEF2LLA
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
 
%   REVISION HISTORY
%   Author                    Date        Comment
%                         (MM/DD/YYYY)
%   Keith Speckman         05/19/2008     Original
%   Ravi Mathur              08/29/2012       Extracted regression test

lat = lat/180*pi;   % rads
lon = lon/180*pi;   % rads
alt = alt*1000;     % m

n = size(lat,2);
x = zeros(3,n);

for k=1:n
    station = jat.groundstations.GroundStation('tmp',lat(k),lon(k),alt(k));
    x(:,k) = station.getECEFPosition();
end


x = x/1000;         % back to km

end