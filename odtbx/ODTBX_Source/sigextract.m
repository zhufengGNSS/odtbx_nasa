function sig = sigextract(Pscrnch)
% SIGEXTRACT  Extract standard deviations from "scrunched" covariances.
%   SIGEXTRACT(P) where P is a rectangular array whose columns are the 
%   unique upper triangular elements of one or more positive semidefinite 
%   matrices as created for example by SCRUNCH.M, will return a rectangular
%   array in which each column contains the square roots of the diagonals
%   of each matrix "scrunched" into P It is most useful for extracting the
%   standard deviations from a sequence of "scrunched" covariances.
%
%   keyword: Utilities,
%
%   See also
%      covariance storage:     SCRUNCH, UNSCRUNCH
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

% Based on Scott Paynter's "unscrunch.m"
% Russell Carpenter

n=size(Pscrnch,1);
m=n;
for k=1:n
   m=m-k;
   if m==0
      n=k;
   end
end
sig(n,:,:)=sqrt(Pscrnch(end,:,:));
sig(1,:,:)=sqrt(Pscrnch(1,:,:));
if n>1
   i=1;
   for j=2:n-1,
      i=i+j;
      sig(j,:,:)=sqrt(Pscrnch(i,:,:));
   end
end
%sig=squeeze(sig);