function M = rotransf(w,D)
% ROTRANSF  Rotational Transformation of Position, Velocity, & Acceleration.
%   M = ROTRANSF(W,D) where W is the constant angular velocity
%   of a frame "B" relative to another frame "A," where D is
%   the direction cosine matrix transforming positions from "A" to "B," and
%   where W is expressed in B's coordinates,  returns a 9x9 matrix, M,
%   that may be used to rotate the coordinates of position,
%   velocity, and acceleration from "A" to "B," and simultaneously to
%   "correct" the velocity and acceleration for the rotation of "B" as
%   shown below.
%
%           A            B
%            d            d               B A
%       v = -- r  ;  v = -- r  =  D * v -  w  x  D * r
%        A  dt  A     B  dt  B   B A   A    B   B A   A
%
%           A 2          B 2
%            d            d                B A
%       a = --- r  ; a = --- r  =  D * v -  w  x  D * v
%        A    2  A    B    2  B   B A   A    B   B A   A
%           dt            dt
%
%   If W and D are 3xN and 3x3xN respectively, then the corresponding
%   output will be 9x9xN.  Note that INV(M) ~= M', but an easy way to get
%   INV(M) is with ROTRANSF(-W,D').
%
%   keyword: Coordinate Transformations, Attitude, Utilities,
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

% Russell Carpenter
% NASA GSFC
% Date: 2005-02-22 12:07:47 -0500 (Tue, 22 Feb 2005)

% Error Checking:
[nw,pw] = size(w);
[nD,mD,pD] = size(D);
if pw~=pD,
    error('ROTRANSF:inpLen','Input lengths not consistent.')
end
if nw~=3,
    error('ROTRANSF:wSize', 'Input angular velocity must be 3x1xN.')
end
if nD~=3 | mD~=3,
    error('ROTRANSF:dcmSize','Input DCM must be 3x3xN.')
end

% Create M:
wX = xmat(w);
M = zeros(9,9,pw);
for k = pw:-1:1,
    wXD(:,:,k) = wX(:,:,k)*D(:,:,k);
    wXwXD(:,:,k) = wX(:,:,k)*wXD(:,:,k);
end
M(1:3,1:3,:) = D;
M(4:6,1:3,:) = -wXD;
M(7:9,1:3,:) =  wXwXD; % Needs to be "+" (prev. vers. had error)
M(4:6,4:6,:) = D;
M(7:9,4:6,:) = -2*wXD;
M(7:9,7:9,:) = D;
