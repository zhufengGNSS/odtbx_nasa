function fail = clk_all_test
% Regression tests for clkdyn2, clkdyn3, clkprop2, and clkprop3
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
% Modification History
% ---------------------
%   Author                  Date            Comments
%   Unknown                 ?               Created
%   Ravi Mathur             08/27/2012      Rename to conform to new
%                                           regression test format

% Test data
t = [0 30];
x = [...
    1e-0 2e-0; 
    1e-3 2e-3; 
    1e-6 2e-6];
opt.q1 = 1.2e-22;
opt.q2 = 1.6e-26;
opt.q3 = 1.0e-35;

% Test cases
[x2d,A2,Q2] = clkdyn2(t,x(1:2,:),opt);
[x3d,A3,Q3] = clkdyn3(t,x,opt);
[x2,rtQd2] = clkprop2(t,x(1:2,2),opt);
[x3,rtQd3] = clkprop3(t,x(:,2),opt);

% Create regression test archive; usually this should be commented out
% x2d_a = x2d; A2_a = A2; Q2_a = Q2;
% x3d_a = x3d; A3_a = A3; Q3_a = Q3;
% x2_a = x2; rtQd2_a = rtQd2;
% x3_a = x3; rtQd3_a = rtQd3;
% save clktests x2d_a A2_a Q2_a x3d_a A3_a Q3_a x2_a rtQd2_a x3_a rtQd3_a

% Load regression archive and test against current values
load clktests.mat
xdtol = eps;
xdfail = any(abs(x2d - x2d_a)>xdtol) + any(abs(x3d - x3d_a)>xdtol);
Atol = eps;
Afail = any(any(abs(A2 - A2_a)>Atol)) + any(any(abs(A3 - A3_a)>Atol));
Qtol = eps;
Qfail = any(any(abs(Q2 - Q2_a)>Qtol)) + any(any(abs(Q3 - Q3_a)>Qtol));
xtol = eps;
xfail = any(abs(x2 - x2_a)>xtol) + any(abs(x3 - x3_a)>xtol);
rtQdtol = eps;
rtQdfail = any(any(abs(rtQd2 - rtQd2_a)>rtQdtol)) ...
    + any(any(abs(rtQd3 - rtQd3_a)>rtQdtol));

fail = any(xdfail + shiftdim(Afail,1) + shiftdim(Qfail,1) + xfail ...
    + shiftdim(rtQdfail,1) > 0);