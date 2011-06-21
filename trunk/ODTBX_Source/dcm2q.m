function q = dcm2q(R)
% 
% Converts the direction cosine matrix to quaternion
%
% Usage:  q = dcm2q(R)
% 
% INPUT: 
%   R      3x3xN      Orthogonal direction cosine matrix
%
% OUTPUT:
%   q      4xN        Quaternion of the form: [vector;scalar]
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

% Written by S. Hur-Diaz

N = size(R,3);
q = NaN(4,N);
for j = 1:N
    q4 = 0.5*sqrt(1+trace(R(:,:,j)));
    q13 = [ R(2,3,j)-R(3,2,j); R(3,1,j)-R(1,3,j); R(1,2,j)-R(2,1,j)]/4/q4;
    qq = [q13;q4];
    q(:,j) = qq/norm(qq);   % make sure it's unit length
end
