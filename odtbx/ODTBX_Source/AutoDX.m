function [dFdX_out, dX_out, dXmax_out, err_out, fcnerr_out] = ...
    AutoDX(F, t0, X0, dX_max, iX, order, varargin)
%AutoDX Numerically compute the optimal step size dX to minimize the error
%in the finite-difference derivative dF/dX of the function F(t,X). The
%derivative is computed with respect to a specified element Xi of X.
%   [DFDX,DX] = AUTODX(F, T, X, DX_MAX, IX, ORDER)
%   F: The function (column vector) whose derivative should be analyzed
%   t: Scalar independent variable, passed to F
%   X: Column vector of dependent variables for differentiation
%   dX_max: Column vector of maximum perturbations for each element of X
%   iX: Element of X to be analyzed. 1 <= i <= length(X)
%   ORDER: Desired finite-difference truncation order. Valid values:
%      ORDER = 1: 1st-order Forward Difference
%      ORDER = 2: 2nd-order Central Difference
%      ORDER = 4: 4th-order Central Difference
%      ORDER = 6: 6th-order Central Difference
%   Additional input arguments are passed to function as F(t,X,varargin)
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
%   Ravi Mathur        11/20/2013      Original (adapted from FORTRAN code)

numdX_exp = 60; % Number of dX exponents to test
npoints = 2; % Number of points (gradients) used to compute estimated truncation error
np = npoints + 1; % Number of saved points (extra one for the "previous" gradient)
valid_orders = [1, 2, 4, 6]; % Possible values for FD method order

if(~any(order == valid_orders))
    error('AutoDX:InvalidOrder', 'Invalid input order');
end

args = varargin;
numX = length(X0);
F0 = feval(F, t0, X0, args{:}); % Function at input X0
numF = length(F0);

slope_errt_starttol = 0.1; % Truncation error slope valid below this value
slope_errt_endtol = 1.0; % Truncation error slope invalid above this value
valid_errt_req = 7; % Number of required consecutive valid truncation error points

% Compute maximum perturbation size based on current parameter's size. We
% do this by seeing which power of 2 approximately equals the parameter's
% maximum value. log2(X) = ln(X)/ln(2) NOTE that for higher order CD
% methods, we will be computing order*dx, so make sure that the maximum
% perturbation size does not exceed dx_max/order.
if(order <= 1)
    X_size = floor(log(dX_max(iX))/log(2.0)) + 1;
else
    X_size = floor(log(dX_max(iX)/(order/2))/log(2.0)) + 1;
end

numdX = 0; % Number of computed gradients
i_dX = np; % This will be incremented into the range [1,np]
err_t_prev = 0.0*F0; % Previous dX value's truncation error (for comparison)
dX_vec = inf(np,1); % Tested dX values
dFdX = zeros(numF, np); % Matrix of gradient vectors for each tested dX value
true_cn = 0.0*F0; % Corrected value of Cn coefficient
slope_errt_best = inf(numF,1); % Best value of truncation error slope
num_valid_errt = 0*F0; % Number of potentially valid truncation errors
num_zero_errt = 0*F0;  % Number of potentially zero truncation errors
errt_valid = false(numF,1); % Indicate that truncation error is invalid
max_valid_dx = 0.0*F0;

% Perturbed function vectors for each dX
% Submatrices contain F(X+-dX), F( X+- 2dX), F(X +- 3dX), stored as column
% vectors for each tested dX.
Fpert = zeros(numF, np, abs(order)); % Perturbed function vectors for each tested dX

for j = 1:numdX_exp
    dX_test = (2.0)^(X_size - j); % Step size for this iteration
    
    % Make sure perturbation is larger than machine precision
    Xplus = X0(iX) + dX_test;
    Xminus = X0(iX) - dX_test;
    if(((Xplus - Xminus) ~= 2.0*dX_test) && (dX_test < 16.0*eps*abs(X0(iX))))
        warning('AutoDX:SkipDX', 'Skipping dX because it is too small to reliably affect X');
        continue
    end
    
    i_dX = advance(i_dX, 1, 1, np); % Increment index in range [1, np]
    dX_vec(i_dX) = dX_test; % Store current dX value
    
    % Get gradient using the appropriate FD method
    [dFdX(:,i_dX), Fvals(:,i_dX,:)] = ...
        getGradient_Order(F, t0, X0, iX, dX_test, order, F0, args{:});
    
    
    
end

dFdX_out = dFdX(:,i_dX);

end % function AutoDX()

%%
function i = advance(i, di, imin, imax)
% Advance integer i by di, wrapping it to the range [imin, imax]. The
% direction of wrapping is consistent with the sign of di. It is assumed
% that imin < imax, and that initially imin <= i <= imax.

i = i + mod(di, imax-imin+1);

if(i < imin)
    i = imax - (imin - i - 1);
elseif(i > imax)
    i = imin + (i - imax - 1);
end

end % function advance

%%
function [dFdX, Fpert] ...
    = getGradient_Order(F, t, X, iX, dx, order, FX, varargin)
% Compute the gradient dF/dX with a given order finite-difference method

Xcurr = X; % Initialize perturbed parameter vector

args = varargin;

if(order == 1) % First-order forward difference method, O(dx)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    dFdX = (Fplus - FX)/dX;
    Fpert(:,1,1) = Fplus;

elseif(order == 2) % Second-order central difference method, O(dx^2)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dx;
    Fminus = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    dFdX = (Fplus - Fminus)/(2.0*dx);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    
elseif(order == 4) % Fourth-order CD method, O(dx^4)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dx;
    Fminus = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 2.0*dx;
    Fplus2 = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 2.0*dx;
    Fminus2 = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    dFdX = (Fminus2 - Fplus2 + 8.0*(Fplus - Fminus))/(12.0*dx);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    Fpert(:,1,3) = Fplus2;
    Fpert(:,1,4) = Fminus2;
    
elseif(order == 6) % Sixth-order CD method, O(dx^6)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dx;
    Fminus = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 2.0*dx;
    Fplus2 = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 2.0*dx;
    Fminus2 = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 3.0*dx;
    Fplus3 = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 3.0*dx;
    Fminus3 = feval(F, t, Xcurr, args{:}); % Backward perturbed 
    
    dFdX = (Fplus3 - Fminus3 ...
            + 9.0*(Fminus2 - Fplus2) ...
            + 45.0*(Fplus - Fminus)) / (60.0*dx);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    Fpert(:,1,3) = Fplus2;
    Fpert(:,1,4) = Fminus2;
    Fpert(:,1,5) = Fplus3;
    Fpert(:,1,6) = Fminus3;
else
    error('order = %d is not supported', order);
end

end % function getGradient_Order

%%
function [err_cond, err_sub] = getErrorCoeff(order, Fcurr, Fpert)
% Computes the coefficients of condition and subtractive cancellation
% errors in dFdX, given perturbations in the current value of F.

error in F given perturbed values

if(order == 1) % O(dx) FD
    Fplus = Fpert(:,1,1);
    err_cond = abs(Fplus) + abs(Fcurr);
    err_sub = max([abs(Fplus), abs(Fcurr)]);

elseif(order == 2) % O(dx^2) CD
    Fplus = Fpert(:,1,1);
    Fminus = Fpert(:,1,2);
    err_cond = (abs(Fplus) + abs(Fminus))/2.0;
    err_sub = max([abs(Fplus), abs(Fminus)])/2.0;

elseif(order == 4) % O(dx^4) CD
    Fplus = Fpert(:,1,1);
    Fminus = Fpert(:,1,2); 
    Fplus2 = Fpert(:,1,3);
    Fminus2 = Fpert(:,1,4);
    err_cond = (abs(Fplus2) + abs(Fminus2) ...
                 + 8.0*(abs(Fplus) + abs(Fminus)))/12.0;
    err_sub = (max([abs(Fplus2), abs(Fminus2)]) ...
                 + 8.0*max([abs(Fplus), abs(Fminus)]))/12.0;

elseif(order == 6) % O(dx^6) CD
    Fplus = Fpert(:,1,1);
    Fminus = Fpert(:,1,2); 
    Fplus2 = Fpert(:,1,3);
    Fminus2 = Fpert(:,1,4);
    Fplus3 = Fpert(:,1,5);
    Fminus3 = Fpert(:,1,6);
    err_cond = (abs(Fplus3) + abs(Fminus3) ...
                 + 9.0*(abs(Fplus2) + abs(Fminus2)) ...
                 + 45.0*(abs(Fplus) + abs(Fminus)))/60.0;
    err_sub = (max([abs(Fplus3), abs(Fminus3)]) ...
                 + 9.0*max([abs(Fplus2), abs(Fminus2)]) ...
                 + 45.0*max([abs(Fplus), abs(Fminus)]))/60.0;

else
    error('order = %d is not supported', order);
end

end % function getErrorCoeff