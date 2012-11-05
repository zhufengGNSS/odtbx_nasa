function [xdot A Q] = dualIADyn(t,x,~)
% dualIADyn Get dynamics & estimation values for the pancake_demo example.

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

I = eye(6,6);
O = zeros(6,6);

xc = x(1:6,:);
xt = x(7:12,:);

[xdotC Ac Qc] = r2bp(t,xc);
[xdotT At Qt] = r2bp(t,xt);

xdot = [xdotC;
        xdotT];

A = nan(12,12,lent);
Q = nan(12,12,lent);
for i=1:lent
    A(:,:,i) = [Ac O;
                O At];
    Q(:,:,i) = blkdiag(Qc,Qt);
end
    
end