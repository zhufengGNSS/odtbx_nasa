function AutoDX_example()
%AUTODX_EXAMPLE An example of using AutoDX to compute optimal step sizes for finite-difference derivatives.
% sizes to minimize error in a finite-difference approximation. This
% example uses a vector function F(t, X) of a vector input X.
%
% See also AutoDX
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Ravi Mathur        10/20/2014       Original


% Set up AutoDX inputs
t0 = 1;               % Reference t
X0 = [30; 45]*pi/180; % Reference X
dX_max = [0.26, 0.26];          % Auto-compute initial large step size for X
dFdX_known = [];      % Assume all gradient elements are unknown

% Display true derivative
F0 = testfcn(t0, X0);
dFdX_true = testfcnder(t0, X0); % True derivative
format long;
disp('True dFdX:');
disp(dFdX_true);

% Use AutoDX to compute optimal step size, Jacobian, etc...
adx = AutoDX;
adx.order = 2; % Use 2nd-order central difference approximation
[dFdX_adx, dX_adx, dXmax_adx, err_adx, fcnerr, nfcalls_adx] = ...
    adx.GetOptimalGradient(@testfcn, t0, X0, F0, dX_max, dFdX_known);

% Compute the true relative error for comparison
err_adx_true = abs((dFdX_adx - dFdX_true)/dFdX_true);

disp('AutoDX Results:');
disp('dFdX:');
disp(dFdX_adx);
disp('optimal dX:');
disp(dX_adx);
disp('max safe dX:');
disp(dXmax_adx);
disp('estimated relative error:');
disp(err_adx);
disp('true relative error:');
disp(err_adx_true);
disp('relative condition error:');
disp(fcnerr);
fprintf('# fcn calls: %i\n\n', nfcalls_adx);

end

% Multivariate example function of scalar t and vector X.
% Note that the input t is treated as a parameter, so it is ignored for the
% purposes of this example.
function F = testfcn(t, X)
F = [sin(X(1))*X(2)^3; ...
    cos(X(2))*X(1)^3];
end

% Analytical derivative, purely for comparison purposes
function dFdX = testfcnder(t, X)
dFdX = [cos(X(1))*X(2)^3, 3*sin(X(1))*X(2)^2; ...
        3*cos(X(2))*X(1)^2, -sin(X(2))*X(1)^3];
end
