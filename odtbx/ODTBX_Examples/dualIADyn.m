function [xdot A Q] = dualIADyn(t,x,~)
% Get dynamics & estimation values for the chaser_target tutorial.
% This function is a wrapper for the individual dynamics of the chaser and
% target spacecraft.
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

lent = length(t);

O = zeros(6,6);

% Extract states from input
xc = x(1:6,:);  % Chaser state
xt = x(7:12,:); % Target state

% Compute accelations, Jacobian, and Process Noise for chaser and target
[xdotC Ac Qc] = r2bp(t,xc); % Chaser state equations of motion
[xdotT At Qt] = r2bp(t,xt); % Target state equations of motion

xdot = [xdotC; xdotT]; % The output dynamics is a concatenation of the chaser and target dynamics

% Concatenate the chaser and target Jacobian and Process Noise
Amat = [Ac O; O At];
Qmat = blkdiag(Qc, Qt);
A = repmat(Amat, [1 1 lent]);
Q = repmat(Qmat, [1 1 lent]);
    
end