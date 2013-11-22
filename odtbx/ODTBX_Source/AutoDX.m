function [dFdX_out, dX_out, dXmax_out, err_out, fcnerr_out] = ...
    AutoDX(F, t0, X0, dX_max, order, dFdX_known, varargin)
%AutoDX Numerically compute the optimal step size dX to minimize the error
%in the finite-difference derivative dF/dX of the function F(t,X). The
%derivative is computed with respect to a specified element Xi of X.
%   [DFDX,DX] = AUTODX(F, T, X, DX_MAX, IX, ORDER)
%   F: The function (column vector) whose derivative should be analyzed
%   t: Scalar independent variable, passed to F
%   X: Column vector of dependent variables for differentiation
%   dX_max: Column vector of maximum perturbations for each element of X
%   dFdX_known: Boolean matrix (length(F)-by-length(X)) 
%   ORDER: Desired finite-difference truncation order. Valid values:
%      ORDER = 1: 1st-order Forward Difference
%      ORDER = 2: 2nd-order Central Difference
%      ORDER = 4: 4th-order Central Difference
%      ORDER = 6: 6th-order Central Difference
%   Additional input arguments are passed through as F(t,X,varargin)
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
%   Ravi Mathur      Aug-Nov 2013      Original adapted from Fortran code
%                                      and Ravi's dissertation

numX = length(X0);
args = varargin;

for iX = numX:1
    [dFdX_out(:,iX), dX_out(:,iX), dXmax_out(:,iX), err_out(:,iX), fcnerr_out(:,iX)] = ...
        AutoDX_iX(F, t0, X0, dX_max, iX, order, dFdX_known, args{:});
end

end

%%
% AutoDX_iX: Performs the actual stepsize-search algorithm on a given
% element of the input X vector (specified with iX)
function [dFdX_out, dX_out, dXmax_out, err_out, fcnerr_out] = ...
    AutoDX_iX(F, t0, X0, dX_max, iX, order, dFdX_known, varargin)

args = varargin;
numX = length(X0);
F0 = F(t0, X0, args{:}); % Function at input X0
numF = length(F0);

numdX_exp = 60; % Number of dX exponents to test
npoints = 2; % Number of points (gradients) used to compute estimated truncation error
np = npoints + 1; % Number of saved points (extra one for the "previous" gradient)

valid_orders = [1, 2, 4, 6]; % Possible values for FD method order
if(~any(order == valid_orders))
    error('AutoDX:InvalidOrder', 'order=%d is not supported', order);
end

dFdX_computed = false(numF,1);
if(~isempty(dFdX_known))
    dFdX_computed(:) = dFdX_known(:,iX);
end

numdX = 0; % Number of computed gradients
i_dX = np; % This will be incremented into the range [1,np]
err_t_prev = 0.0*F0; % Previous dX value's truncation error (for comparison)
dX_vec = inf(np,1); % Tested dX values
dFdX = zeros(numF, np); % Matrix of gradient vectors for each tested dX value
Cn_true = 0.0*F0; % Corrected value of Cn coefficient
slope_errt_best = inf(numF,1); % Best value of truncation error slope
num_valid_errt = 0*F0; % Number of potentially valid truncation errors
num_zero_errt = 0*F0;  % Number of potentially zero truncation errors
errt_valid = false(numF,1); % Indicate that truncation error is invalid
max_valid_dx = 0.0*F0;

slope_errt_starttol = 0.1; % Tolerance for valid truncation error slope
slope_errt_endtol = 1.0; % Tolerance for invalid truncation error slope
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

% Perturbed function vectors for each dX
% Submatrices contain F(X+-dX), F(X+- 2dX), F(X +- 3dX), stored as column
% vectors for each tested dX.
Fpert = zeros(numF, np, abs(order)); % Perturbed function vectors for each tested dX

for j = 1:numdX_exp
    if(all(dFdX_computed))
        break;
    end
    
    dX_test = (2.0)^(X_size - j); % Step size for this iteration
    
    % Make sure perturbation is larger than machine precision
    Xplus = X0(iX) + dX_test;
    Xminus = X0(iX) - dX_test;
    if((Xplus - Xminus) ~= 2.0*dX_test)
        warning('AutoDX:SkipDX', 'Skipping dX=%d because it is too small to reliably affect X', dX_test);
        continue
    end
    
    i_dX = advance(i_dX, 1, 1, np); % Increment index in range [1, np]
    dX_vec(i_dX) = dX_test; % Store current dX value
    
    % Get gradient using the appropriate FD method
    [dFdX(:,i_dX), Fpert(:,i_dX,:)] = ...
        getGradient_Order(F, t0, X0, iX, dX_test, order, F0, args{:});
    
    if(~all(isfinite(dFdX(:,i_dX))))
        warning('AutoDX:NotFinite', 'Skipping dX because some elements are not finite');
        i_dX = advance(i_dX, -1, 1, np); % Reset index
        continue
    end
    
    numdX = numdX + 1; % Increment number of valid gradient points
    
    % Make sure that there are enough points to start approximating
    % truncation error
    if(numdX < npoints) 
        continue
    end
    
    % Shift the 'current' point back (will be reset later). This allows us
    % to work with the "current point, next point" wording instead of
    % "current point, previous point".  It's mainly just to be consistent
    % with the wording of the derived equations.
    i_dX = advance(i_dX, -(npoints-1), 1, np);
    i_dX_prev = advance(i_dX, -1, 1, np);
    i_dX_next = advance(i_dX, 1, 1, np);
    
    % Compute error estimates independently for each element of F(t,X)
    % Note that this (probably) cannot be vectorized since each element of
    % F (each subfunction) may have a different optimal step size.
    for currF = numF:1
        if(dFdX_computed(currF))
            continue;
        end
        
        f1 = dFdX(currF, i_dX);      % Derivative from current point
        f2 = dFdX(currF, i_dX_next); % Derivative from next point
        dx1 = dX_vec(i_dX);          % Perturbation from current point
        dx2 = dX_vec(i_dX_next);     % Purturbation from next point
        
        % Compute truncation error for the current point and load into the
        % error vector. Here we have a special case: if dF/dX is not
        % changing wrt dX, then the Cn value should be zero. BUT because of
        % numerical errors, we need to consider "zero" as some max cutoff.
        Cn = f2 - f1;
        err_t = 10.0*eps*max(abs(f1), abs(f2)); % Cutoff for Cn
        if(abs(Cn) <= err_t)
            Cn = 0; % Change in dFdX is negligible, so Cn should be zero
        else
            Cn = Cn/(dx1^order - dx2^order); % Otherwise, compute normally
        end
        err_t = Cn * dx1^order; % Estimated truncation error
        
        % Check if truncation error is actually O(dx^order). To do so, we
        % test the slope of log(|err_t|) wrt log(dX), and test whether it
        % is near the analytically predicted value.
        if((err_t == 0.0) || (err_t_prev(currF) == 0.0))
            slope_errt = inf;
        else
            % slope = [log(|err_t_prev|) - log(|err_t|)]/[log(dX_prev) - log(dX)]
            %       = log(|err_t_prev/err_t|) / log(dX_prev/dX)
            slope_errt = log(abs(err_t_prev(currF)/err_t))/log(dX_vec(i_dX_prev)/dx1);
            slope_errt = slope_errt - order; % Slope error wrt expected value
        end
        err_t_prev(currF) = err_t; % Update truncation error history
        
        % Test for a valid nonzero truncation error
        if(abs(slope_errt) <= order*slope_errt_starttol)
            num_zero_errt(currF) = 0; % No more nonzero truncation errors
            
            % For the first pair of valid dX points, both are valid
            if(num_valid_errt(currF) == 0)
                num_valid_errt(currF) = 2;
            else
                num_valid_errt(currF) = num_valid_errt(currF) + 1;
            end
            
            % The first dX whose truncation error slope is nearly correct
            % is the maximum safe dX
            if(max_valid_dx(currF) == 0.0)
                max_valid_dx(currF) = dX_vec(i_dX_prev);
            end
                                                    
            % If the current truncation error slope is closer to ideal than
            % the previous slope, then use the current Cn value as the
            % correct one. Note that Cn is proportional to the (order+1)-th
            % derivative of the function.
            if(abs(slope_errt) < abs(slope_errt_best(currF)))
                Cn_true(currF) = Cn;
                slope_errt_best(currF) = slope_errt; % Update best slope with current one
            end
            
            % If we have a sufficient number of consecutive valid points,
            % then we are on the way to computing an optimal dX value. Note
            % that "sufficent number" is currently empirically obtained,
            % and must be at least 2.
            if(num_valid_errt(currF) == max(2,valid_errt_req - floor(order/2)))
                errt_valid(currF) = true;
            end
            
        % Even if the current truncation error is not valid from a log-log
        % slope point of view, it may actually be zero. A zero truncation
        % error corresponds to the function behaving like a polynomial of
        % degree O(order).
        elseif(~errt_valid(currF))
            num_valid_errt(currF) = 0; % No more valid nonzero truncation errors
            slope_errt_best(currF) = inf;
            
            % Check for zero truncation error
            if(err_t == 0.0)
                num_zero_errt(currF) = num_zero_errt(currF) + 1;
                
                % First zero truncation error
                if(num_zero_errt(currF) == 1) 
                    max_valid_dx(currF) = dx1;

                % Truncation error is zero
                elseif(num_zero_errt(currF) == (valid_errt_req - floor(order/2)))  
                    % If the computed dF/dX = 0, then there is no way to
                    % tell if the derivative is just zero at the current X,
                    % or if it is zero at every X (i.e. constant function).
                    % So the max valid dX is set to zero to tell the user
                    % to redo this analysis when X changes.
                    if(abs(f1) < eps)
                        max_valid_dx(currF) = 0.0;
                    end
                    
                    errt_valid(currF) = true;
                    Cn_true(currF) = 0.0;
                end

            % Otherwise truncation error estimates are just plain invalid,
            % so reset the counters to prepare for possible upcoming valid
            % truncation errors.
            else
                max_valid_dx(currF) = 0.0;
                num_zero_errt(currF) = 0;
                Cn_true(currF) = 0.0;
            end
        end
        
        % If the truncation error is valid, but its slope goes beyond an
        % allowable tolerance, then we assume that roundoff error has crept
        % in enough to consider this the 'optimal' dX value.
        if(errt_valid(currF) && (abs(slope_errt) > order*slope_errt_endtol))
            
            % Estimate the current function's condition error assuming the
            % current dX value minimizes total truncation+roundoff errors.
            
            i_dX_opt = i_dX_prev; % Indicate optimal dX index
            
            % Compute the corrected optimal step size
            t_ratio = dx2/dx1;
            t_ratio = ((1.0 - t_ratio^order)/(1.0 + 1.0/t_ratio))^(1.0/(order+1));
            hopt_corrected = t_ratio*dx1;
            %dX_vec(i_dX_opt) = hopt_corrected;
            
            % Estimate maximum subtractive cancellation error
            [err_cond, err_sub] = getErrorCoeff(order, F0(currF), Fpert);
            err_sub = eps*err_sub;
            
            % Estimate maximum condition error for this function
            if(Cn_true(currF) == 0.0)
                err_a = 0.0; % Cn = 0 means that dFdX is not changing wrt dX, so assume that function has no error
            else
                err_a = order*abs(Cn_true(currF))*dX_vec(i_dX_opt)^(order+1) - err_sub;
            end
            
            % Store the estimated condition error for this function
            if(err_cond == 0.0)
                fcnerr_out(currF) = 1.0; % Indicate a large relative error
            else
                fcnerr_out(currF) = err_a/err_cond; % Convert to relative error
            end
            
            % Store the current dX value as the optimal dX value
            dX_out(currF) = dX_vec(i_dX_opt);
            
            % Store the gradient from the current dX value as the optimal gradient
            dFdX_out(currF) = dFdX(currF,i_dX_opt); % Uncorrected gradient
            
            % Store the max dX value
            dXmax_out(currF) = max_valid_dx(currF);
            
            % Store the estimated roundoff + truncation error
            err_out(currF) = (err_a + err_sub)/dX_vec(i_dX_opt) + abs(Cn_true(currF))*dX_vec(i_dX_opt)^order;
            
            if(dFdX(currF,i_dX) ~= 0.0) % Relative error if possible
                err_out(currF) = err_out(currF)/(dFdX(currF,i_dX_opt) + err_out(currF));
            end
            
            dFdX_computed(currF) = true; % Done analyzing this element of dFdX
            errt_valid(currF) = false; % Truncation error is no longer valid
        end
                
    end
    
    i_dX = advance(i_dX, npoints-1, 1, np); % Reset the 'current' point
    
end

end % function AutoDX_iX()

%%
function iout = advance(i, di, imin, imax)
% Advance integer i by di, wrapping it to the range [imin, imax]. The
% direction of wrapping is consistent with the sign of di. It is assumed
% that imin < imax, and that initially imin <= i <= imax.

iout = i + mod(di, imax-imin+1);

if(iout < imin)
    iout = imax - (imin - iout - 1);
elseif(iout > imax)
    iout = imin + (iout - imax - 1);
end

end % function advance

%%
function [dFdX, Fpert] ...
    = getGradient_Order(F, t, X, iX, dX, order, FX, varargin)
% Compute the gradient dF/dX with a given order finite-difference method

Xcurr = X; % Initialize perturbed parameter vector

args = varargin;

if(order == 1) % First-order forward difference method, O(dx)
    Xcurr(iX) = X(iX) + dX;
    Fplus = F(t, Xcurr, args{:}); % Forward perturbed
    
    dFdX = (Fplus - FX)/dX;
    Fpert(:,1,1) = Fplus;

elseif(order == 2) % Second-order central difference method, O(dx^2)
    Xcurr(iX) = X(iX) + dX;
    Fplus = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dX;
    Fminus = F(t, Xcurr, args{:}); % Backward perturbed
    
    dFdX = (Fplus - Fminus)/(2.0*dX);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    
elseif(order == 4) % Fourth-order CD method, O(dx^4)
    Xcurr(iX) = X(iX) + dX;
    Fplus = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dX;
    Fminus = F(t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 2.0*dX;
    Fplus2 = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 2.0*dX;
    Fminus2 = F(t, Xcurr, args{:}); % Backward perturbed
    
    dFdX = (Fminus2 - Fplus2 + 8.0*(Fplus - Fminus))/(12.0*dX);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    Fpert(:,1,3) = Fplus2;
    Fpert(:,1,4) = Fminus2;
    
elseif(order == 6) % Sixth-order CD method, O(dx^6)
    Xcurr(iX) = X(iX) + dX;
    Fplus = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dX;
    Fminus = F(t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 2.0*dX;
    Fplus2 = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 2.0*dX;
    Fminus2 = F(t, Xcurr, args{:}); % Backward perturbed
    
    Xcurr(iX) = X(iX) + 3.0*dX;
    Fplus3 = F(t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - 3.0*dX;
    Fminus3 = F(t, Xcurr, args{:}); % Backward perturbed 
    
    dFdX = (Fplus3 - Fminus3 ...
            + 9.0*(Fminus2 - Fplus2) ...
            + 45.0*(Fplus - Fminus)) / (60.0*dX);
    Fpert(:,1,1) = Fplus;
    Fpert(:,1,2) = Fminus;
    Fpert(:,1,3) = Fplus2;
    Fpert(:,1,4) = Fminus2;
    Fpert(:,1,5) = Fplus3;
    Fpert(:,1,6) = Fminus3;
else
    error('AutoDX:getGradient_Order:InvalidOrder', 'order = %d is not supported', order);
end

end % function getGradient_Order

%%
function [err_cond, err_sub] = getErrorCoeff(order, Fcurr, Fpert)
% Computes the coefficients of condition and subtractive cancellation
% errors in dFdX, given perturbations in the current value of F.

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
    error('AutoDX:getErrorCoeff:InvalidOrder', 'order = %d is not supported', order);
end

end % function getErrorCoeff