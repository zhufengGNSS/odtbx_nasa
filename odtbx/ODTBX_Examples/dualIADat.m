function [y H R] = dualIADat(t,x,options)
% dualIADat Data file for chaser_target_demo example.

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

sig = options.sig;

lent = length(t);
I = eye(3,3);
O = zeros(3,3);

xc = x(1:3,:);
xt = x(7:9,:);

q = repmat([0 0 0 1]',1,lent);
[yg Hg Rg] = gpsmeas(t,x(1:6,:),options,q);

n = size(yg,1);

yp = xt-xc;

y = [yp;yg];

for i=lent:-1:1
    H(:,:,i) = [-I O I O; Hg(:,:,i) zeros(n,6)];
    R(:,:,i) = blkdiag(diag(sig.^2),Rg(:,:,i));
end

end