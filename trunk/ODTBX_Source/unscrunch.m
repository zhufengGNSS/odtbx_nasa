function A = unscrunch(x,sigflag)
% UNSCRUNCH  Recreate a matrix that has been vectorized with SCRUNCH.
%    A = UNSCRUNCH(X) reconstructs a symetric matrix which is in
%    the vector form, X, consistent with "scrunch.m".
%
%    SIG = UNSCRUNCH(X,SIGFLAG) with SIGFLAG = 1 returns only the square
%    roots of the diagonal elements of the original matrix, i.e. the
%    standard deviations, if the original matrix was a covariance.
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

% Russell Carpenter (added SIGFLAG and generalized for 3-D arrays)
% NASA GSFC 
%   (for other changes, see the svn repository)

if nargin==1,
    sigflag=0;
end
n=size(x,1);
m=n;
for k=1:n
    m=m-k;
    if m==0
        n=k;
    end
end
for el = size(x,2):-1:1
    A(1,1,el)=x(1,el);
    if n>1
        range = 2:n;
        if sigflag,
            i=1;
            for j=range
                i=i+j;
                A(j,el)=x(i,el);
            end
        else
            for k=range
                m=sum(1:k-1);
                A(1:k,k,el)=x(m+1:m+k,el);
            end
        end
    end
    if sigflag
        A(:,:,el)=sqrt(A(:,:,el));
    else
        A(:,:,el)=A(:,:,el)+A(:,:,el)'-diag(diag(A(:,:,el)));
    end
end