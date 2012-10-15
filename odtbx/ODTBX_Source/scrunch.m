function x = scrunch(A)
% SCRUNCH  "Vectorize" a square symmetric matrix.
%   X = SCRUNCH(A) takes the upper triangular part of a square symmetric
%   matrix, A, and puts it in vector form, starting with the first column.
%   It takes the elements in a column from the top to the diagonal element
%   and continually appends it to the vector.
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

% Scott Paynter
% University of Texas at Austin, ~1991

% Russell Carpenter (added input warnings and generalized for 3-D arrays)
% NASA Goddard Space Flight Center

[n,m,p]=size(A);
if n~=m,
    ns=num2str(n);
    warning('SCRUNCH:notSquare', ...
        ['Input not square; upper left ',ns,'x',ns,' will be SCRUNCHED.'])
    A=A(1:n,1:n,:);
end
dA = abs(A-permute(A,[2 1 3]));
if any(dA>eps),
    warning('SCRUNCH:notSymmetric', ...
        ['Input not symmetric; upper triangular part will be SCRUNCHED.\n',...
        '         max(A-A'') = ', num2str(max(max(dA)))])
end
for i=p:-1:1
    y=A(1,1,i);
    if n>1
        for k=2:n
            y=[y;A(1:k,k,i)];
        end
    end
    x(:,i) = y;
end
