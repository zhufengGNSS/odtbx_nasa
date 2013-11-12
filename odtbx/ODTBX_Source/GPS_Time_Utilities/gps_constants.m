% CONSTANTS       - Script loads constants used in GPS/orbital sims
%
% Last modified: 5/2002

L1 = 1575.42e6;    % Hz
L2 = 1227.6e6;
L5 = 1176.45e6;
C = 2.99792458e8;  % m/s       % Speed of light from GPS ICD 1991

chipRateL1CA=1.023e6;
L1CAperiod=1024/chipRateL1CA;
chipRateL1PY=10.23e6;


% Reference: 
% WGS-84: "World Geodetic System 1984"; the ECEF spatial coordinate system
% used by GPS since 22 Jan 1987 (see, for example, "Session 2: World Geodetic
% System 1984" pp. 67-134, in Proceedings of the Fourth International Geodetic
% Symposium on Satellite Positioning, Apr. 28 - May 2, 1986, Volume 1);
MU = 3.986005e14;              % Earth's grav. parameter [m^3/s^2] from GPS ICD 1991
EARTH_RATE = 7.2921151467e-5;  % WGS-84 Earth's rotation rate [rad/s] from GPS ICD 1991
EARTH_RADIUS = 6378137.0;      % WGS-84 Earth's sma [m] 
EARTH_B = 6356752.3142;        % WGS-84 Earths semiminor axis [m]

EARTH_GM = 3.986005e14;        % Earth's grav. parameter [m^3/s^2] from GPS ICD 1991

SECONDS_IN_WEEK = 7*24*3600;
SECONDS_IN_DAY = 24*3600;

% References: NIMA WGS Definitions, wgs84fin.pdf, 2 Jan 2000, pp3-1 to 3-7
% A more recent definition of the WGS ECEF coordinate frame?
if 0
EARTH_RATE = 7.292115e-5;      % Earth's rotation rate [rad/s] from NIMA 2000 spec
EARTH_GM = 3986004.418e8;      % Earths gravitational parameter [m^3/s^2] from NIMA 2000 spec
EARTH_RADIUS = 6378137.0;      % Earth's sma [m] from NIMA 2000 spec
EARTH_B = 6356752.3142;        % Earths semiminor axis [m] from NIMA 2000 spec
EARTH_ECC = 8.1819190842622e-2;  % Earths first eccentricity from NIMA 2000 spec
end

% defined in the ICD_GPS-200 as: 
% *       origin is the center of mass of the Earth 
% *       z-axis is parallel to the direction of the Conventional
%    International Origin (CIO) as defined by the Bureau International de'l Heure
%    (BIH), and passes through instantaneous pole of epoch 1984.0 
% *       x-axis is the intersection of the reference meridian plane and the
%    plane of the mean astronomic equator, with the reference meridian being
%    parallel to the zero meridian defined by the BIH 
% *       y-axis completes the system as a right-handed rectangular system 
% Some of the major constants in WGS-84, which models the Earth as an
% ellipsoid of revolution, are: 
% *       semimajor axis (origin to equator on x-y plane) = 6378.137 km 
% *       semiminor axis (origin to either pole) = 6356.7523142 km 
% *       flattening = 1/298.257223563 
% *       angular velocity of Earth = 7.292115e-5 radians/s 
% *       angular velocity of Earth (untruncated) = 7.2921151467e-5 radians/s 
%         (for precise satellite applications) 
% *       G x mass of Earth = 3.986005e14 m**3/s**2 
