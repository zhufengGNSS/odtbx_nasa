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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Ravi Mathur        11/20/2013      Original (adapted from FORTRAN code)

%% Start regression test for AutoDX
function [failed] = AutoDX_test

% Test function and true derivative
F = @(t, X) X^3.5;
dFdX = @(t, X) 3.5*X^2.5;

% AutoDX inputs
t0 = 0;             % Reference t
X0 = 9.1;           % Reference X
dX_max = abs(X0)+1; % Max step size for X
iX = 1;             % Element of X for which to evaluate derivative
order = 2;          % Order of finite-difference derivative method

[dFdX_opt, dX_opt, dXmax, err, fcnerr] = ...
    AutoDX(F, t0, X0, dX_max, iX, order, [])
dFdX_true = dFdX(t0, X0)
dFdX_err = abs(dFdX_opt - dFdX_true)/dFdX_true

tol = 1.0e-10;

failed = (dFdX_err > tol);

end