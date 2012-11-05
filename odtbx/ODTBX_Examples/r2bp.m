function [xd,A,Q] = r2bp(t,x,mu)
% R2BP  Restricted Two-body problem dynamics equations.
%   xd = r2bp(t,x,mu) returns the derivatives for the restricted
% two-body problem.  If the gravity coefficient, mu = G*M is ommited, the
% default (EGM-96) value (in km^3/sec^2) is used.
%   [xd,A] = r2bp(t,x,mu) also returns the Jacobian matrix, A.
%   [xd,A,Q] = r2bp(t,x,mu) also returns a process noise power spectral
% density based on a random walk process noise input.  At present the
% magnitude of the PSD is hardcoded to be 1e-9, where the units will depend
% on the value of mu (ie 1e-9 km^2/sec if default mu is used).
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

% Version History:
% Original    Kathryn Gregory, based on Russell Carpenter's pr2bp
% Revision 1: Russell Carpenter - Cleanup and corrections
% (for additional changes, see the svn repository)

m = length(t);
if nargin == 2 || isempty(mu),
    mu = 3.986004415000000e+5; % EGM-96, km^3/sec^2
end
R = x(1:3,:);
[Ur,r] = unit(R);
g = mu./r.^2;
V = x(4:6,:);
xd = [V; repmat(-g,3,1).*Ur];
if nargout > 1,
   I = repmat(eye(3),[1 1 m]);
   UUT = bsxfun(@times,reshape(Ur,3,1,m),reshape(Ur,1,3,[]));
   G = bsxfun(@times,reshape(-g./r,1,1,m),(I - 3*UUT));
   A(4:6,1:3,:) = G;
   A(1:3,4:6,:) = I;
   %Above does the following, but is often faster.
   %I = eye(3);
   %for k = m:-1:1,
   %    A(4:6,1:3,k) = -g(k)./r(k)*(I - 3*Ur(:,k)*Ur(:,k)');
   %    A(1:3,4:6,k) = I;
   %end
end
if nargout == 3,
    Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2),[1 1 m]);
end
