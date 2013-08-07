function [y H R] = dualIADat(t,x,options)
% Get measurement values for the chaser_target tutorial.
% This function is a wrapper for 2 measurement functions: a GPS measurement
% (via the ODTBX gpsmeas function) and a 'pose' measurement (using range).
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

sig = options.sig; % Extract sigma from input options

lent = length(t);
I = eye(3,3);
O = zeros(3,3);

xc = x(1:3,:); % Chaser state
xt = x(7:9,:); % Target state
yp = xt - xc;  % Compute the 'pose' range measurement

% Compute the chaser's state via GPS
q = repmat([0; 0; 0; 1],1,lent);
[yg Hg Rg] = gpsmeas(t,x(1:6,:),options,q);

y = [yp;yg]; % The output measurement is a concatenation of the pose and GPS measurements

% Compute the measurement partials and noise covariance
pad = zeros(size(yg,1), 6);
for i=lent:-1:1
    H(:,:,i) = [-I O I O; Hg(:,:,i) pad];       % Measurement partials
    R(:,:,i) = blkdiag(diag(sig.^2),Rg(:,:,i)); % Measurement noise covariance
end

end