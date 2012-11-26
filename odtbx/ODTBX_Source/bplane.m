function [Bt, Br, LTOF] = bplane(X,N,mew)
% Function to calculate B-Plane parameters for a hyperbolic orbit
%
%    [Bt Br LTOF] = bplane(X,N,mew) returns the B-Plane coordinates Bt, Br,
%    and Linearized Time of Flight (LTOF) for the given spacecraft state, X.
%    The reference direction is defined by the unit vector N.  The central
%    body gravitational parameter is defined in the variable mew.
%
%      INPUTS                 DESCRIPTION
%       X                      Position and Velocity (hyperbolic)
%       N                      Reference Direction {[0 0 1]'}
%       mew                    Gravitational Parameter {3.986004415e+5}
%
%keyword: Coordinate Transformations, Interplanetary, Targetting
% See also: DCM
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

% Kenneth M. Getzandanner
% Navigation & Mission Design Branch
% NASA Goddard Space Flight Center

if isempty(N)
    N = [0 0 1]';
end

if isempty(mew)
    mew = 3.986004415000000e+5;
end

% Ensure N is a unit vector
N = N/norm(N);

% Get position and velocity
P = X(1:3,1);
V = X(4:6,1);

% Position and velocity magnitude
p = norm(P);
v = norm(V);

% Angular momentum vector
h = cross(P,V)/norm(cross(P,V));

% Eccentricity
E = (v^2/mew-1/p)*P - V*dot(P,V)/mew;
e = norm(E);

% Semi-major axis
a = 1/(2/p-v^2/mew);

% B-vector magnitude
b = -a*sqrt(e^2-1);

% Beta angle
beta = pi/2-asin(1/e);

% True anomaly (position and B-vector)
theta = acos(-a*(e^2-1)/(p*e)-1/e);
thetab = pi-(pi/2+beta);

% Hyperbolic excess velocity vector
Vinf = cos(beta)*E/e+sin(beta)*cross(h,E/e)/norm(cross(h,E/e));

% Construct B-Plane coordinate frame
S = Vinf;

T = cross(S,N)/norm(cross(S,N));

R = cross(S,T)/norm(cross(S,T));

B = b*cross(S,h)/norm(cross(S,N));

% Calculate BdotT and BdotR
Bt = dot(B,T);
Br = dot(B,R);

% Hyperbolic anomaly (position and B-vector)
Fr = acosh((e+cos(theta))/(1+e*cos(theta)));
Fb = acosh((e+cos(thetab))/(1+e*cos(thetab)));

% Time since periapsis (position and B-vector)
tr = (e*sinh(Fr)-Fr)/sqrt(mew/(-a)^3);
tp = (e*sinh(Fb)-Fb)/sqrt(mew/(-a)^3);

% Calculate Linearized Time of Flight
LTOF = tr-tp;

end