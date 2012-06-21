function Po = initialize_cov(r,v,sig_pos,sig_vel,sig_rad,sig_spd,mu)
% INITIALIZE_COV Creates an initial covariance matrix
% 
% Po = init_cov(r,v,sig_r,sig_v,sig_e,sig_h,mu)
% 
% INITIALIZE_COV was designed to include pseudo-measurements of certain
% orbital quantities in the initial covariance matrix. A spherical initial
% covariance matrix (sigma values along the diagonal) is unrealistic and
% INITIALIZE_COV is needed to make the covariance more realistic. See
% Richard Battin's book "An Introduction to the Mathematics and Methods
% of Astrodynamics" pgs 678-679 for more detailed information.
% 
% INPUTS
% VARIABLE      SIZE    DESCRIPTION
%    r          3x1     Inertial starting position of spacecraft (km)
%    v          3x1     Inertial starting velocity of spacecraft (km/s)
%    sig_pos    1       Standard deviation of positions (km)
%    sig_vel    1       Standard deviation of velocities (km/s)
%    sig_rad    1       Standard deviation of total radius (km)
%    sig_spd    1       Standard deviation of total speed (km/s)
%    mu         1       Gravitational constant of central body (km^3/s^2)
% 
% OUTPUTS 
%    Po                 Initial covariance matrix
% 
% 
% keyword: initialize, initial, covaraince 
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
% 
%   REVISION HISTORY
%    Author               Date         	    Comment
%    Benjamin Asher       06/22/2012        Original
% 

%% Check inputs

if nargin < 6
    err = MException('CovarianceInitializer:NotEnoughInputs', ...
        'Not enough input arguments.');
    throw(err)
end
if nargin == 6
    mu = 3.986004418e5;                             % km^3/s^2
end

%% Calculate initial covariance matrix

R = norm(r);
V = norm(v);
h = cross(r,v);
H = norm(h);
be = [r*(mu/R^3);v];
bh = (1/H)*[cross(v,h);cross(h,r)];
Prv = diag([sig_pos^2*ones(1,3) sig_vel^2*ones(1,3)]);
sig_e = mu/R^2*sig_rad+V*sig_spd;
sig_h = (V*sig_rad+R*sig_spd);
W = (diag([sig_e sig_h]))\eye(2);
Ho = [be';bh'];
Po = (Prv\eye(6)+Ho'*W*Ho)\eye(6);

end