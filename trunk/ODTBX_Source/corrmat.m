function C = corrmat(P)
% CORRMAT  Correlation Matrix.
%   C = CORRMAT(P) with P = a positive semi-definite, symmetric covariance
%   matrix P, returns a matrix C where
%   (1) the diagonals of C are the square roots of the diagonals of P, i.e.
%       the diagonal of C contains the standard deviations;
%   (2) the upper triangle of C contains the correlation coefficients; and
%   (3) the lower triangle of C is a copy of the lower triangle of P.
%
%keyword: Utilities,
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

% Error Checking
s2 = diag(P);
if any(s2<0),
    error('CORRMAT:notPosSemiDef', 'Input is not positive semidefinite.')
end
[n,m] = size(P);
if n~=m,
    error('CORRMAT:notSquare', 'Input is not square.')
end
if any(abs(P-P')>eps),
    error('CORRMAT:notSymmetric', 'Input is not symmetric to within eps.')
end

% Place sigmas along main diagonal of C:
C = diag(sqrt(s2));

% Compute correlation coefficients for upper triangle of C:
n = length(P);
for i = 1:n-1, % Loop over rows
    for j = i+1:n, % Loop over columns
        CiiCjj = C(i,i)*C(j,j);
        if CiiCjj, % i.e. if ~= 0
            C(i,j) = P(i,j)/CiiCjj;
            if C(i,j) < -1 | C(i,j) > 1,
                error('CORRMAT:corrTooBig', ...
                    'Non-physical correlation; maybe input not >= 0')
            end
        else,
            C(i,j) = 0;
        end
    end
end

% Copy lower triangle of P into lower triangle of C
C = C + tril(P,-1);
