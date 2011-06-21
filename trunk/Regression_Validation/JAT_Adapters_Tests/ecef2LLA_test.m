function failed = ecef2LLA_test()
% Regression Test Case
% Function(s) ecef2LLA
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

failed  = 0;
tol     = 1e-10;

x = [   0     0   6500    0
        0    6500   0     0
        6500  0     0   -6500 ];
    
latTarget = [   0                       90       0         0 ];
lonTarget = [   90                       0       0       -90 ];
altTarget = [  1.432476857548207e+002 121.863 121.863 1.432476857548207e+002 ]; 
    
[lat,lon,alt] = ecef2LLA(x);

if( (any(abs(latTarget-lat))>tol) || (any(abs(lonTarget-lon))>tol) || (any(abs(altTarget-alt))>tol) )
    failed = 1;
end