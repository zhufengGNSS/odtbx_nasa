function R=q2dcm(q)
%
% Converts from a quaternion [vector;scalar] to a direction cosine matrix
%
% Usage:  R=q2dcm(q)
% 
% INPUT:
%   q    4xN      Quaternion of the form: [vector scalar]
%
% OUTPUT: 
%   R    3x3xN    Orthogonal direction cosine matrix
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

%
% Written by Sun Hur-Diaz
%
[m,N] = size(q);
if m ~= 4
    error('Quaternion input dimensions need to be 4xN')
end

R = NaN(3,3,N);
for i = 1:N
    qi = q(:,i)/norm(q(:,i));
    R(:,:,i) = (qi(4)^2-qi(1:3)'*qi(1:3))*eye(3)+2*qi(1:3)*qi(1:3)'-2*qi(4)*cr(qi(1:3));
end

function xc = cr(x)
xc = [ 0   -x(3),   x(2)
       x(3)  0,    -x(1)
      -x(2)  x(1)    0 ];