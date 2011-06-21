function [y,H,R] = rrdot2D(t,x,sig)
% RRDOT2D  Range and range-rate measurement model for 2-D inertial state
% space.
%   [y,H,R] = rrdot2D(t,x,sig)
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

if nargin == 2 || isempty(sig),
    sig = diag([1e-3 1e-6]); % 1 m, 1 mm/sec
end
rp = 6378;
omega = 2*pi/86400;
theta = omega*t(:)';
rs = rp*[cos(theta);sin(theta)];
vs = omega*rp*[-sin(theta);cos(theta)];
dr = x(1:2,:)-rs;
dv = x(3:4,:)-vs;
rho = sqrt(sum(dr.^2));
y = [rho; dot(dr,dv)./rho];
udr = unit(dr);
if nargout > 1,
    Hvv = shiftdim(udr,-1);
    Hrr = shiftdim(udr,-1);
    Hvr = shiftdim(dv.*repmat(1./rho,2,1) ...
        - dr.*repmat(dot(dr,dv)./rho.^3,2,1),-1);
    Hrv = zeros(size(Hvr));
    H = [Hrr, Hrv; Hvr, Hvv];
    %Above does this:
    %for k = length(t):-1:1,
    %    H(2,4:6,k) = udr(:,k)';
    %    H(2,1:3,k) = dv(:,k)'/rho(k) - dv(:,k)'*dr(:,k)*dr(:,k)'/rho(k)^3;
    %    H(1,1:3,k) = udr(:,k)';
    %end
end
if nargout == 3,
    R = repmat(sig.^2,[1,1,length(t)]);
end