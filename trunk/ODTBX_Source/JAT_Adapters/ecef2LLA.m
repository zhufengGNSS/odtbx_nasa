function [lat,lon,alt] = ecef2LLA(x,ErrTol)

% ECEF2LLA  Returns the latitude, longitude and altitude of an ECEF position.
%
%   [lat,lon,alt] = ecef2LLA(x) returns the latitude, longitude and
%   altitude of input ECEF positions. Latitude and longitude both range
%   from -180 to 180 (deg).
%
%   INPUTS
%   VARIABLE    SIZE        DESCRIPTION (Optional/Default)
%     x         (3xN)       ECEF Positions in km
%     ErrTol    (1x1)       (optional) error (km) above which a warning is
%                           produced; default = 1e-9 km
%
%   OUTPUTS
%     lat       (1xN)       Latitudes in degs
%     lon       (1xN)       Longitudes in degs
%     alt       (1xN)       Altitudes in km
%
% VALIDATION/REGRESSION TEST
%
%   These tests were moved to ecef2LLA_test.m to conform to the new
%   regression testing framework.
%
%   keyword: Coordinate Transformations, JAT Adapter
%   See also JATDCM
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
%   Author      		Date         	Comment
%                           (MM/DD/YYYY)
%   Derek Surka              09/18/2007         Original
%   Keith Speckman           05/30/2008       Corrected error in JAT call
%                                             Modified JAT algorithm
%                                             Added regression and validation
%   Allen Brown              3/30/2009        Fixed regression check logic,
%                                             fixed regression data error
%                                             (now using
%                                             ecef2LLA_RegressionData6_08b)
%   Ravi Mathur              08/29/2012       Extracted regression test


if nargin == 1
        ErrTol = 1e-9; % km
end

x = x*1000; % m

n   = size(x,2);
lat = zeros(1,n);
lon = zeros(1,n);
alt = zeros(1,n);

for k=1:n
    LLA = jat.groundstations.GroundStation.getLLA(x(:,k));
    lat(1,k) = LLA(1)*180/pi;
    lon(1,k) = LLA(2)*180/pi;
    alt(1,k) = LLA(3)/1000; %km
    ecef_error = LLA(4:6); %m

    % If the error excedes the tolerance, produce a warning
    % This may be due to the JAT iteration scheme not converging resulting in a small error
    if norm(ecef_error) > ErrTol*1e3
	    warning(['Error exceeds tolerance.  Error = ',num2str(norm(ecef_error/1e3)),' km'])
    end
end

end