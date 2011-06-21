function [xd,A,Q] = pr2bp(t,x,mu)
% PR2BP  Planar restricted two-body problem dynamics equations.
%   [xd,A,Q] = pr2bp(t,x,mu)
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

if nargin == 2 || isempty(mu),
    mu = 3.986004415000000e+5; % EGM-96, km^3/sec^2
end
R = x(1:2,:);
[Ur,r] = unit(R);
g = mu./r.^2;
V = x(3:4,:);
xd = [V; repmat(-g,2,1).*Ur];
if nargout > 1,
    I = repmat(eye(2),[1 1 length(t)]);
    UUT = repmat(permute(shiftdim(Ur,-1),[2 1 3]),[1 2]) .* ...
        repmat(shiftdim(Ur,-1),[2 1]);
    G = repmat(shiftdim(-g./r,-1),[2 2]).*(I - 3*UUT);
    A(3:4,1:2,:) = G;
    A(1:2,3:4,:) = I;
    %Above does the following, but is about 3x faster.
    %I = eye(2);
    %for k = length(t):-1:1,
    %    A(3:4,1:2,k) = -g(k)./r(k)*(I - 3*Ur(:,k)*Ur(:,k)');
    %    A(1:2,3:4,k) = I;
    %end
end
if nargout == 3,
    Q = repmat(diag([0 0 1e-9 1e-9].^2),[1 1 length(t)]);
end