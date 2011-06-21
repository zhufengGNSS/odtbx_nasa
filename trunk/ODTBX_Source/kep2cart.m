function X = kep2cart(KOE,GM)
% KEP2CART  Convert Keplerian orbital elements to cartesian states.
%   X = KEP2CART(KOE) converts the two-body Keplerian orbital elements in
%   the structure KOE to Cartesian position/velocity vectors, X = [R;V],
%   using the EGM-96 value of the earth's gravitational constant, 
%   GM=3.986004415e+14 m^3/sec^2.  KOE is defined as in KEPEL, as follows:
%      KOE.sma = semi-major axis
%      KOE.ecc = eccentricity
%      KOE.incl = inclination
%      KOE.raan = right ascension of the ascending node
%      KOE.argp = argument of periapse
%      KOE.tran = true anomaly
%   All angles are in radians, and SMA is in the length units of GM (meters
%   by default).  The output state, X, will be expressed in the same
%   planet-centered inertial frame to which the elements are referenced.
%   The output will be a vector if KOE is a "scalar" structure;" otherwise
%   it will be a matrix where each column is the position/velocity state of
%   the corresponding set of osculating elements in the fields of KOE.
%
%   X = KEP2CART(KOE,GM) uses the input GM instead of the EGM-96 value.
%
% keyword: Utilities, 
%
% See also
%      Converting from Cartesian states to Keplerian elements: KEPEL
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

% Bob Merriam (el2xyz)
% NASA JSC

% Russell Carpenter
% NASA JSC and NASA GSFC

if nargin == 1,
    GM=3.986004415e+14;
end

% el2xyz code (vectorized by RC):
a = [KOE(:).sma]         ;% semi-major axis
e = [KOE(:).ecc]         ;% eccentricity
i = [KOE(:).incl]        ;% inclination (radians)
p = [KOE(:).argp]        ;% argument of perigee (radians)
n = [KOE(:).raan]        ;% ascending node (radians)
w = [KOE(:).tran]        ;% true anomaly (radians)

sp = sin(p);   cp = cos(p);
sn = sin(n);   cn = cos(n);
si = sin(i);   ci = cos(i);
sw = sin(w);   cw = cos(w);

P = [ cn.*cp-sn.*sp.*ci           ;% compute position of periapsis
      sn.*cp+cn.*sp.*ci           ;
      sp.*si             ]        ;

Q = [-cn.*sp-sn.*cp.*ci           ;% compute the vector at a right angle
     -sn.*sp+cn.*cp.*ci           ;% to P in the direction of motion
      cp.*si             ]        ;

slr = a.*(1 - e.^2);              % semi-latus rectum
rmag = slr./(1 + e.*cw);          % position vector magnitude

% Position and velocity vectors:
r(1:3,:) = [ P*diag(cw) + Q*diag(sw)     ] * diag(rmag);
v(1:3,:) = [-P*diag(sw) + Q*diag(e + cw) ] * diag(sqrt(GM./slr));
X = [r;v];
