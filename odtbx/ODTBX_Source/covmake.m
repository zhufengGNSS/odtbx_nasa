function P = covmake(C)
% COVMAKE  Covariance maker.
%   P = COVMAKE prompts the user to input standard deviations and
%   correlation coefficients, and creates the corresponding covariance
%   matrix.
%
%   P = COVMAKE(C) creates the covariance from the correlation matrix C,
%   where the main diagonal of C contains the standard deviations, and the
%   upper triangle of C contains the correlation coefficients.  The lower
%   triangle of C is ignored.
% 
% keyword: Utilities,
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

% CREATED:
% Russell Carpenter 
% NASA Johnson Space Center, ~1990
% UPDATED:
% Russell Carpenter 
% NASA Goddard Space Flight Center ~2005

if nargin == 0,
    n = input('Enter dimension: ');
    for k = 1:n,
        sig(k) = input(['Enter sig_', num2str(k), ': ']); %#ok<AGROW>
    end
    for j = 1:(n-1),
        for k = (j+1):n,
            rho(j,k) = input(['Enter rho_', num2str(j), num2str(k), ': ']); %#ok<AGROW>
        end
    end
elseif nargin == 1,
    [m,n] = size(C);
    if m~=n,
        error('COVMAKE:notSquare', 'Input is not square.')
    end
    rho = C;
    sig = diag(C);
else
     error('COVMAKE:wrongInput', 'COVMAKE must have zero or one input.')
end

n = length(sig);
for k = n:-1:1,
    P(k,k) = sig(k)^2;
end
for j = 1:(n-1),
    for k = (j+1):n,
        if rho(j,k) >= 1 || rho(j,k) <= -1,
            error('COVMAKE:rhoBig', ['-1<rho_', num2str(j), num2str(k), ...
                '<+1 is not satisified.'])
        else
            P(j,k) = rho(j,k)*sig(j)*sig(k); %#ok<AGROW>
            P(k,j) = P(j,k); %#ok<AGROW>
        end
    end
end
e = eig(P);
if any(e < 0),
    error('COVMAKE:nonPos',...
        'Covariance non-positive; try different correlations.');
end