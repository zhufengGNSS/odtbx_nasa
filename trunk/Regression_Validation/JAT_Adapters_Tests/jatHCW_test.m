function failed = jatHCW_test()
% Regression Test Case
% Function(s) jatHCW
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
orbRate = 1.1e-3; % rad/sec
dT      = 10; %step size in sec
t       = 0; % sec
x       = [10;1;20;0.3;0.4;0.5]; % km

targetXDot = [  3.000000000000000e-001
                4.000000000000000e-001
                5.000000000000000e-001
                9.163000000000000e-004
               -6.600000000000000e-004
               -2.420000000000000e-005]; % km, sec


targetA    = [  1.000181498169882e+000    0                        0    9.999798334553415e+000  1.099988908377773e-001                         0
               -1.330991947469384e-006    1                        0   -1.099988908377773e-001  9.999193338213656e+000                         0
                         0                0   9.999395006100392e-001                         0                       0    9.999798334553415e+000
                3.629926795442890e-005    0                        0    9.999395006100392e-001  2.199955633601751e-002                         0
               -3.992959737411317e-007    0                        0   -2.199955633601751e-002  9.997580024401569e-001                         0
                         0                0  -1.209975598480963e-005                         0                       0    9.999395006100392e-001 ];

[xDot, a] = jatHCW(t,x,{orbRate,dT});

if( any( abs(xDot-targetXDot) > tol ) )
    failed = 1;
end
if( any( any( abs(a-targetA) > tol ) ) )
    failed = 1;
end

