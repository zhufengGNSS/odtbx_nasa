function x = chi2inv_odtbx(p,v)
% ODTBX version of chi2inv
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

% Sun Hur-Diaz
% July 19, 2010

x = 7;
maxcount = 200;
i = 1;
dxmin = 1e-6;

while i <= maxcount
    [po,dpdx] = chi2cdf_odtbx(x,v);
    dx = (p-po)/dpdx;
    dx = (.2/v)*dx;
    if x + dx < 0
        x = x/2;
    else
        x = x + dx;
    end
    if norm(dx) < dxmin
        break;
    end
    i = i + 1;
end
if i>maxcount
    disp('Max iteration in chi2inv_odtbx reached.')
end
