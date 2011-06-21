function [y,H,R] = rrdot3D_1way(t,x,sig)
%RRDOT3D_1way Range and range-rate measurement model for 3-D inertial state
%   [y,H] = rrdot3D_1way(t,x,sig) returns the measurement, y, and 
% measurement partials, H, for 1-way range and range-rate from an 
% artificial ground station located at [rp,0,0], where rp, the radius of 
% the earth, and the earth rotation rate (omega in the code) are based on 
% the parameters used in the NASA GSFC Orbit Determination Error Analysis 
% System. If the measurement noise standard deviation matrix, sig, is 
% ommited, a default value of sig = diag([1e-3 1e-6]) is used.  
% NOTE: The value of rp, and the default values of sig, are specified in km
% and km/sec.  The state is position and velocity in km and km/s.
%   [y,H,R] = rrdot3D_1way(t,x,sig) also outputs measurement noise 
% covariance matrix.
%
% See Also:  RRDOT3D
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
% Original by Kathryn Gregory, based on Russell Carpenter's rrdot2D.
% Revised by Kathryn Gregory & Kevin Berry - Updated function for batch
% validation example.
% Revised by Russell Carpenter - Cleanup and corrections.
% Sun Hur-Diaz   7/30/09    Removed the 2-multiplier for 1-way

lent = length(t);
if nargin < 3 || isempty(sig),
    sig = diag([1e-3 1e-6]); % 1 m, 1 mm/sec
end
rp = 6378.1363; % km
omega = 7.2921158494e-05; % Rotation Rate used in ODEAS
theta = omega*t(:)';
o = zeros(size(theta));
rs = rp*[cos(theta); sin(theta); o];
vs = omega*rp*[-sin(theta); cos(theta); o];
dr = x(1:3,:) - rs;
dv = x(4:6,:) - vs;
rho = sqrt(sum(dr.^2));
% The next 3 lines implement a zero-deg elev mask
y = NaN(2,lent);
k = (dot(dr,rs)>=0);
y(:,k) = [rho(k); dot(dr(:,k),dv(:,k))./rho(:,k)];
udr = unit(dr);
if nargout > 1,
    Hvv = shiftdim(udr,-1);
    Hrr = shiftdim(udr,-1);
    Hvr = shiftdim(dv.*repmat(1./rho,3,1) ...
        - dr.*repmat(dot(dr,dv)./rho.^3,3,1),-1);
    Hrv = zeros(size(Hvr));
    H = [Hrr, Hrv; Hvr, Hvv];
    % Above does this:
    %for k = lent:-1:1,
    %   H(2,4:6,k) = 2*udr(:,k)';
    %   H(2,1:3,k) = 2*dv(:,k)'/rho(k) - dv(:,k)'*dr(:,k)*dr(:,k)'/rho(k)^3;
    %   H(1,1:3,k) = 2*udr(:,k)';
    %end
end
R = repmat(sig.^2,[1,1,lent]);
