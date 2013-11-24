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
%   Ravi Mathur        11/20/2013      Original

%% Start regression test for AutoDX
function [failed] = AutoDX_test

% Test function and analytical derivative
F = @(t, X) (exp(X)/sqrt(sin(X^3) + cos(X^3)));
dFdX = @(t, X) (exp(X)*((3*X^2+2)*sin(X^3) + (2-3*X^2)*cos(X^3))) /...
               (2*(sin(X^3)+cos(X^3))^(3/2));
           
% AutoDX inputs
t0 = 0;             % Reference t
X0 = 1.33;          % Reference X (near a singularity)
dX_max = [];        % Auto-compute initial large step size for X
order = 2;          % 2nd-order central difference method

F0 = F(t0, X0);
dFdX_true = dFdX(t0, X0); % True derivative
fprintf('True dFdX             = %d\n\n', dFdX_true);

% Use AutoDX to compute optimal step size and associated derivative
adx = AutoDX;
adx.order = order;
[dFdX_adx, dX_adx, dXmax_adx, err_adx, fcnerr] = ...
    adx.GetOptimalDX(F, t0, X0, F0, dX_max, []);
err_adx_true = abs((dFdX_adx - dFdX_true)/dFdX_true);
fprintf('AutoDX dFdX           = %d\n', dFdX_adx);
fprintf('AutoDX Optimal dX     = %d\n', dX_adx);
fprintf('AutoDX Max Safe dX    = %d\n', dXmax_adx);
fprintf('AutoDX Function Error = %d\n', fcnerr);
fprintf('AutoDX estimated err  = %d\n', err_adx);
fprintf('AutoDX true err       = %d\n\n', err_adx_true);

% Use numjac to compute derivative
[dFdX_nj, fac] = numjac(F, t0, X0, F0, 10, [], false);
err_nj_true = abs((dFdX_nj - dFdX_true)/dFdX_true);
fprintf('numjac dFdX           = %d\n', dFdX_nj);
fprintf('numjac true err       = %d\n', err_nj_true);

tol = 1.0e-10;
failed = (err_adx_true > tol);

end