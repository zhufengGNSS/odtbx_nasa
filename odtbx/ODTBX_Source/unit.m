function [u,l] = unit(v)
% UNIT  Efficiently compute unit vectors and norms for vector arrays.
%   U = UNIT(V) where V is a vector returns the unit vector in the
%   direction of V.  The output U will have the same row or column
%   orientation as the input V.  If the input V is a matrix, UNIT
%   interprets each column of the matrix to be a separate vector, and
%   returns a matrix in U that has the unit vectors of the columns of V as
%   its columns.
%
%   [U,L] = UNIT(V) also returns the length (2-norm) L of the vector V.
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
% NASA GSFC, 08/30/04
%   (for other changes, see the svn repository)

[m,n] = size(v);
if m==1 && n==1, % Input is a scalar
    warning('UNIT:scalar', 'Input is a scalar.')
    l = v;
    u = 1;
elseif xor(m==1,n==1), % One or the other but not both, i.e. vector input
    l = norm(v);
    u = v/l;
else % Matrix input; treat as collection of column vectors
    l = sqrt(sum(v.^2));
    % u = v.*repmat(1./l,m,1);
    u = bsxfun(@times,v,1./l); % Faster for large n
end
