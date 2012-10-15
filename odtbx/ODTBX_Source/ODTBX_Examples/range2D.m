function [y,H,R] = range2D(t,x,sig)
% RANGE2D  Range measurement model for 2-D inertial state space.
%   [y,H,R] = range2D(t,x,sig)
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
    sig = 1e-3; % 1 m
end
theta = 2*pi/86400*t(:)';
rs = 6378*[cos(theta);sin(theta)];
dr = x(1:2,:)-rs;
k = find(dot(dr,rs)>0); % This is the zero-deg elev mask
y = NaN(1,length(t));
y(k) = sqrt(sum(dr(:,k).^2)); % Only pos elev gets a value
if nargout > 1,
    Hrr = shiftdim(unit(dr),-1);
    Hrv = zeros(size(Hrr));
    H = [Hrr, Hrv];
    %Above does this:
    %for k = length(t):-1:1,
    %    H(:,:,k) = [unit(dr(:,k))', zeros(1,2)];
    %end
end
if nargout == 3,
    R = repmat(sig^2,[1,1,length(t)]);
end
