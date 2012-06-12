function x = covsmpl(P,varargin)
% COVSMPL Covariance sampling function.
%    COVSMPL(P) returns 1 sample of a zero-mean Gaussian distribution with
%    covariance P.  P must be real and positive definite. 
% 
%    COVSMPL(P,k) returns k samples, with each sample in a column.
%
%    COVSMPL(P,k,s) returns k samples and sets the random number state to s
%    (see the documentation or help on RANDN for more info on how Matlab
%    handles random number seeds). NOTE: Passing in s = NaN will not
%    affect the random number state, i.e., is equivalent to not passing in
%    s at all.
%
%    COVSMPL(P,'principalAxis') returns "plus one sigma" principal axis
%    perturbation vectors derived from an eigenvalue decomposition of P.
%    (To get +/-n sigma pertubations, multiply the output by +/-n.)
%
%    NOTE: If P is a 3-D array, with each "slice" representing a different
%    covariance matrix, COVSMPL will return a set of samples for each
%    covariance.  The output will be squeezed to remove any singleton
%    dimensions, so that in such cases, the output of COVSMPL(P) will be
%    2-D, but the output of COVSMPL(P,...) will be 3-D.
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

% Russell Carpenter
% NASA Goddard Space Flight Center, 12/24/2003.

% Kevin Berry   12/05/2008  Added check for infinite covariance
% R Carpenter   05/05/2009  Fixed principal axis option & updates for 7.7

         
% Parse inputs.
rndsmp = true; % Default is to return random samples.
if nargin > 2,
    if(~isnan(varargin{2}))
        warning('COVSMPL:seedReset', 'Resetting random number seed')
        try % new syntax for ML 7.7+
%             RandStream('mcg16807', 'Seed', varargin{2})
            RandStream.setGlobalStream(RandStream('shr3cong','Seed',varargin{2}));
        catch %#ok<CTCH>
            randn('state', varargin{2}) %#ok<RAND>
        end
    end
end
if nargin > 1,
    if ischar(varargin{1}),
        rndsmp = false;
        if strcmpi(varargin{1}, 'principalaxis')
            warning('COVSMPL:princAx', ...
                'Returning principal axis perturbations');
            k = 1;
        end
    else
        k = varargin{1};
    end
else
    k = 1;
end
if nargin > 3 || nargin == 0,
    error('COVSMPL:numInputs', 'Unsupported number of inputs')
end
% Compute samples.
[nr,nc,ns] = size(P);
x = zeros(nr,k,ns);
for j = ns:-1:1,
    fin = ~isinf(diag(P(:,:,j))); %finite elements along the diagonal
    Pj  = P(fin,fin,j); %finite section of the covariance matrix 
    nj  = sum(fin); %dimension of Pj  
    if rndsmp,
        % Use Cholesky decomposition to get the "matrix square root" of P.
        % Need to use the transpose because of Matlab's convention that
        % R'*R = P, where R = chol(P).  If P has any zeros eigenvalues,
        % this won't work, so use the eigenvalue decomposition instead, and
        % force any negative eigenvalues to be zero.
        try
            randarray = RandStream.getGlobalStream.randn(nj,k);
        catch
            randarray = randn(nj,k);
        end
        try
            x(fin,:,j) = chol(Pj)'*randarray; 
        catch %#ok<CTCH>
            [V,D] = eig(Pj);
            D(D<0) = 0;
            x(fin,:,j) = V*sqrt(D)*randarray; 
        end
    else
        % Find "plus one sigma" principal axis perturbations.
        [V,D] = eig(Pj);
        x = V*sqrt(D); 
    end
end
if rndsmp,
    x = squeeze(x);
    if nr == 1,
        x = x';
    end
end
