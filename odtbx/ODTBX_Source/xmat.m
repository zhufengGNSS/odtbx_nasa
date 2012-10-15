function x = xmat(v)
% XMAT  Cross-product matrix for arrays of vectors.
%   X = XMAT(V) where V is a 3x1 vector returns the skew-symmetric
%   "cross-product matrix" X such that X*W = CROSS(V,W):  
%      X = [  0   -V(3)  V(2)
%            V(3)   0   -V(1)
%           -V(2)  V(1)   0  ];
%   If V is 3xN, then N 3x3 cross-product matrices, one for each column of
%   V, are returned.
%
%   keyword: Utilities, 
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
% NASA Johnson Space Center, ~1990

% Russell Carpenter
% NASA Goddard Space Flight Center
% Date: 2005-02-22 12:07:47 -0500 (Tue, 22 Feb 2005)

[m,n] = size(v);
if m == 1 & n == 3,
    el = 1;
    v = v';
elseif m == 3,
    el = n;
else,
    error('XMAT:tooBig', 'Input must be 3xN, 3x1, or 1x3.')
end
x = zeros(3,3,el);
x(1,2,:) = -v(3,:);
x(2,1,:) =  v(3,:);
x(1,3,:) =  v(2,:);
x(3,1,:) = -v(2,:);
x(2,3,:) = -v(1,:);
x(3,2,:) =  v(1,:);
