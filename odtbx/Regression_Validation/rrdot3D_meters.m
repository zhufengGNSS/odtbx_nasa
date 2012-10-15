function [y,H,R] = rrdot3D_meters(t,x,sig)
% RRDOT3D_METERS  Range and range-rate measurement model for 3-D inertial state
% space.
%   [y,H,R] = rrdot3D(t,x,sig)
% This measurement file was used for the ODTBX batch validation against ODEAS. The omega listed is based on the
% Rotation Rate used in ODEAS. The station is initially located at [rp,0,0], where rp is the radius 
% of the Earth. This case uses a single polar orbit spacecraft which is initially directly over 
% the station. The measurement step is 10 seconds.
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

% Revision History:
% Original   Kathryn Gregory, based on Russell Carpenter's rrdot2D
% Revision 1: Kathryn Gregory & Kevin Berry  
% Date: 2-12-07  Updated function for batch validation example

x = x/1000;  %DMS convert incoming metres to km to be consistent with other code

if nargin == 2 || isempty(sig),
    sig = diag([1e-3 1e-6]); % 1 m, 1 mm/sec
end
rp = 6378.1363;
%omega = 2*pi/86164;
omega = 7.2921158494e-05; %Rotation Rate used in ODEAS
theta = (omega*t(:)');
yc = zeros(1,length(theta)); % Need this because need only first row otherwise too large and screws up dimensions of rs 
rs = rp*[cos(theta);sin(theta); yc(1,:)];

vs = omega*rp*[-sin(theta);cos(theta); yc(1,:)];
dr = x(1:3,:)-rs;
dv = x(4:6,:)-vs;
rho = sqrt(sum(dr.^2));
y = NaN(2,length(t));
k = find(dot(dr,rs)>=0); % This is the zero-deg elev mask
for(i=1:length(k))
    y(:,k(i)) = [rho(k(i)); dot(dr(:,k(i)),dv(:,k(i)))./rho(k(i))];
end
%y(k) = sqrt(sum(dr(:,k).^2)); % Only pos elev gets a value
%y(:,k) = [rho(k); dot(dr(:,k),dv(:,k))./rho(k)];

udr = unit(dr);
if nargout > 1,
    Hvv = shiftdim(udr,-1);
    Hrr = shiftdim(udr,-1);
    Hvr = shiftdim(dv.*repmat(1./rho,3,1) ...
        - dr.*repmat(dot(dr,dv)./rho.^3,3,1),-1);
    Hrv = zeros(size(Hvr));
    H = 2*[Hrr, Hrv; Hvr, Hvv];
    %Above does this:
%     for k = length(t):-1:1,
%        H(2,4:6,k) = udr(:,k)';
%        H(2,1:3,k) = dv(:,k)'/rho(k) - dv(:,k)'*dr(:,k)*dr(:,k)'/rho(k)^3;
%        H(1,1:3,k) = udr(:,k)';
%     end
end
if nargout == 3,
    R = repmat(sig.^2,[1,1,length(t)]); %Original Line
  % R = repmat(sig.^2,[1,1,length(y)]);
end
