function varargout = kepel(r,v,GM)
% KEPEL  Compute Keplerian orbital elements from Cartesian states.
%   KOE = KEPEL(R,V) computes the two-body Keplerian orbital elements
%   from input position vectors R and V, using the EGM-96 value of the
%   earth's gravitational constant, GM = 3.986004415e+14 m^3/sec^2.  
%   r and v units are meters.
%
%   KOE = KEPEL(R,V,GM) uses the input GM instead of the EGM-96 value.
%   r and v units must match GM.
%
%   KOE = KEPEL(X,[],GM) assumes that X=[R;V].
%   r and v units must match GM.
%
%   The elements are returned in the data structure KOE which has the
%   following fields, referenced to the basis in which R & V reside:
%      KOE.sma = semi-major axis
%      KOE.ecc = eccentricity
%      KOE.incl = inclination
%      KOE.raan = right ascension of the ascending node
%      KOE.argp = argument of periapse
%      KOE.tran = true anomaly
%   All angles are in radians, and SMA is in the units of R and GM.
%
%   If called with six output arguments instead of one, the KOE's are
%   returned in the order listed above.
%
% keyword: Utilities
%
% See also
%      Converting back to Cartesian states: KEP2CART
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

% Steven Hughes
% NASA GSFC

% Russell Carpenter
% NASA GSFC

% Parse inputs and check for errors:
if nargin < 3,
    GM = 3.986004415e+14;
end
[nr,nc] = size(r);
if ~xor(nr==1,nc==1), % NOT(Either but not both)
    error('KEPEL:inp1NotVec', 'First input must be a vector.')
end
el = length(r);
if nargin == 1 || isempty(v),
    if el ~= 6,
        error('KEPEL:xNot6d', 'X must be 6x1 or 1x6.')
    end
    v = r(4:6);
    r = r(1:3);
else
    [nrv,ncv] = size(v);
    if ~xor(nrv==1,ncv==1) || length(v)~=3,
        error('KEPEL:vNot3d', 'V must be 3x1 or 1x3.')
    end
    if el ~= 3,
        error('KEPEL:rNot3d', 'R must be 3x1 or 1x3.')
    end
    
    % now that the lengths are checked, ensure r and v have the same
    % orientation:
    if nr ~= nrv
        v = v'; % force v orientation to match r
    end
end

% Compute elements (cart2oe.m code):
k = [0;0;1];
h = cross(r,v);
n = cross(k,h);
N = norm(n);
H2 = dot(h,h);
V2 = dot(v,v);
R = norm(r);
e = ( (V2-GM/R)*r - (dot(r,v))*v )/GM;
p  = H2/GM;

KOE.ecc = norm(e);
KOE.sma = p/(1 - KOE.ecc^2);
KOE.incl = acos(h(3)/sqrt(H2));
KOE.raan = acos(n(1)/N);
if n(2) < -eps % Fix quadrant.
   KOE.raan = 2*pi - KOE.raan;
end
KOE.argp = acos(dot(n',e)/N/KOE.ecc);
if e(3) < -eps % Fix quadrant.
   KOE.argp = 2*pi - KOE.argp;
end
KOE.tran = acos(dot(e',r)/KOE.ecc/R);
if dot(r,v) < -eps	% Fix quadrant.
   KOE.tran = 2*pi - KOE.tran;
end
KOE = orderfields(KOE,[2 1 3:6]);
if nargout == 1
    varargout{1} = KOE;
else
    varargout{1} = KOE.sma;
    varargout{2} = KOE.ecc;
    varargout{3} = KOE.incl;
    varargout{4} = KOE.raan;
    varargout{5} = KOE.argp;
    varargout{6} = KOE.tran;
end 