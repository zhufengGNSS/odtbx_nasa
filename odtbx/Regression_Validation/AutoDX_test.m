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
%   Ravi Mathur        11/24/2013      Original

%% Start regression test for AutoDX
function [failed] = AutoDX_test(order)

% Test function and analytical derivative
F = @(t, X) (exp(X)/sqrt(sin(X^3) + cos(X^3)));
dFdX = @(t, X) (exp(X)*((3*X^2+2)*sin(X^3) + (2-3*X^2)*cos(X^3))) /...
               (2*(sin(X^3)+cos(X^3))^(3/2));
           
% AutoDX inputs
t0 = 0;             % Reference t
X0 = 1.33;          % Reference X (near a singularity)
dX_max = [];        % Auto-compute initial large step size for X
dFdX_known = [];    % Assume all gradient elements are unknown

if(nargin == 0)
    order = 1; % Default to 1st-order forward difference method
end

F0 = F(t0, X0);
dFdX_true = dFdX(t0, X0); % True derivative
fprintf('True dFdX             = %d\n\n', dFdX_true);

% Use AutoDX to compute optimal step size and associated derivative
adx = AutoDX;
adx.order = order;
[dFdX_adx, dX_adx, dXmax_adx, err_adx, fcnerr, nfcalls_adx] = ...
    adx.GetOptimalGradient(F, t0, X0, F0, dX_max, dFdX_known);
err_adx_true = abs((dFdX_adx - dFdX_true)/dFdX_true);
fprintf('AutoDX dFdX           = %d\n', dFdX_adx);
fprintf('AutoDX optimal dX     = %d\n', dX_adx);
fprintf('AutoDX max safe dX    = %d\n', dXmax_adx);
fprintf('AutoDX rel cond err   = %d\n', fcnerr);
fprintf('AutoDX rel est err    = %d\n', err_adx);
fprintf('AutoDX true err       = %d\n', err_adx_true);
fprintf('AutoDX # fcn calls    = %i\n\n', nfcalls_adx);

% Use numjac to compute derivative
[dFdX_nj, ~, ~, ~, nfcalls_nj] = numjac(F, t0, X0, F0, 0);
err_nj_true = abs((dFdX_nj - dFdX_true)/dFdX_true);
fprintf('numjac dFdX           = %d\n', dFdX_nj);
fprintf('numjac true err       = %d\n', err_nj_true);
fprintf('numjac # fcn calls    = %i\n', nfcalls_nj);

tol = 1.0e-6;
failed = (err_adx_true > tol);

end