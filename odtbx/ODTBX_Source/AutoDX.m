classdef AutoDX < handle
    %AutoDX Numerically compute the optimal step size dX to minimize the error
    %in the finite-difference derivative dF/dX of the function F(t,X).
    %   [DFDX,DX] = AutoDX.GetOptimalGradient(F, T0, X0, F0, DX_MAX, DFDX_KNOWN)
    %   dFdX: Jacobian matrix (length(F)-by-length(X)) 
    %   dX: Matrix of optimal dX values. Element (i,j) corresponds to the
    %       optimal dX value for subfunction Fi and variable Xj. If the
    %       values down a particular column vary significantly, then some
    %       subfunctions of F may be ill-conditioned wrt that Xj variable.
    %   F: The vector function whose derivative should be analyzed.
    %   t0: Scalar independent variable, passed to F
    %   X0: Vector of dependent variables for differentiation
    %   F0: Vector F(t0, X0)
    %   dX_max: Vector of maximum perturbations for each element of X
    %           Set to [] to use dX_max = 1 + abs(X) (a safe default)
    %   dFdX_known: Boolean matrix (length(F)-by-length(X)). Can be []
    %   ORDER: Desired finite-difference truncation order. Valid values:
    %      ORDER = 1: 1st-order Forward Difference
    %      ORDER = 2: 2nd-order Central Difference
    %      ORDER = 4: 4th-order Central Difference
    %      ORDER = 6: 6th-order Central Difference
    %   Additional input arguments are passed through as F(t,X,varargin)
    %
    %   [DFDX,DX,DXMAX,ERR,FCNERR,NFCALLS] = AutoDX.GetOptimalGradient(...)
    %   dXmax: Matrix of maximum safe step sizes. Element (i,j) corresponds
    %          to the maximum dX value for subfunction Fi and variable Xj.
    %          If the values down a particular column vary significantly,
    %          then some subfunctions of F may be stiff wrt that Xj.
    %   err: Matrix of estimated errors in corresponding dFdX elements.
    %        Errors are upper-bound magnitudes relative to dFdX.
    %   fcnerr: Matrix of estimated condition errors in corresponding
    %           subfunctions. Element (i,j) corresponds to the condition
    %           error in function Fi wrt variable Xj. Condition errors are
    %           relative, and imply the number of digits in Fi that are
    %           reliable given a small change in Xj. If some values down a
    %           particular column are large, then those subfunctions of F
    %           may be losing precision wrt that Xj.
    %           e.g. fcnerr(2,3)=1e-9 means that subfunction F2 has ~9
    %           digits of accuracy wrt a small change in variable X3, and
    %           has lost the remaining ~6 digits of accuracy due to errors
    %           building up within its implementation.
    %   nfcalls: Total number calls made to F(t,X).
    %
    %   Example:
    %      adx = AutoDX;  % Initialize AutoDX object
    %      adx.order = 2; % Specify 2nd-order Central Differences
    %      [...] = adx.GetOptimalGradient(...); % Optimize
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
    
    properties
        order; % Order of finite-difference method (see constant valid_orders)
    end
    
    % Internal working arrays
    properties (SetAccess = private, GetAccess = private)
        numF;            % Number of elements in function output
        err_t_prev;      % Previous dX value's truncation error (for comparison)
        dX_vec;          % Tested dX values
        dFdX;            % Matrix of gradient vectors for each tested dX value
        Cn_true;         % Corrected value of Cn coefficient
        slope_errt_best; % Best value of truncation error slope
        num_valid_errt;  % Number of potentially valid truncation errors
        num_zero_errt;   % Number of potentially zero truncation errors
        errt_valid;      % Indicate that truncation error is invalid
        max_valid_dx;    % Maximum valid step size
        dFdX_computed;   % Keep track of successfully computed derivatives
        
        % Fpert: Perturbed function vectors for each tested dX. Submatrix i
        % contains F(X+-i*dXj), stored as column vectors for each tested
        % dXj where 1 <= j <= np (see Constants for np).
        Fpert;           

    end
    
    properties (Constant)
        valid_orders = [1, 2, 4, 6]; % Valid finite-difference method orders
        numdX_exp = 60;            % Number of dX exponents to test
        npoints = 2;               % Number of points (gradients) used to compute estimated truncation error
        np = AutoDX.npoints + 1;   % Number of saved points (extra one for the "previous" gradient)
        slope_errt_starttol = 0.1; % Tolerance for valid truncation error slope
        slope_errt_endtol = 1.0;   % Tolerance for invalid truncation error slope
        valid_errt_req = 7;        % Number of required consecutive valid truncation error points
    end
    
    methods
        function obj = AutoDX()
            obj.order = 2; % Default to O(h^2) central differences
            obj.numF = nan;
        end
        
        %% set.order 
        % Check that desired order is valid when it is set.
        % This function is called by MATLAB whenever 'order' is set
        % explicitely through a 'obj.order = value' assignment.
        function obj = set.order(obj, value)
            if(~any(value == obj.valid_orders))
                warning('AutoDX:InvalidOrder', 'order=%d is not supported, using %d instead', value, obj.order);
            else
                obj.order = value;
            end
        end
        
        %% GetOptimalGradient
        function [dFdX_out, dX_out, dXmax_out, err_out, fcnerr_out, numFCalls] = ...
                GetOptimalGradient(obj, F, t0, X0, F0, dX_max, dFdX_known, varargin)
            %GetOptimalGradient Compute the optimal step size for each
            %element of X
            % See also: AutoDX
                    
            numX = length(X0);
            if(isempty(dX_max))
                dX_max = abs(X0) + 1.0;
            end
            
            numFCalls = 0;
            
            for iX = numX:1
                [dFdX_out(:,iX), dX_out(:,iX), dXmax_out(:,iX), err_out(:,iX), fcnerr_out(:,iX), nFc] = ...
                    obj.AutoDX_iX(F, t0, X0, F0, dX_max, iX, dFdX_known, varargin{:});
                numFCalls = numFCalls + nFc;
            end
            
        end
        
        %% GetGradient
        function dFdX_out = GetGradient(obj, F, t0, X0, F0, dX, varargin)
            %GetGradient Compute the gradient dF/dX of a function F(t,X)
            % DFDX = GETGRADIENT(F, T0, X0, F0, DX)
            % dFdX: Matrix of dF/dX elements
            % F:  Function handle for F(t, X)
            % t0: Scalar independent variable, passed to F
            % X0: Vector of dependent variables for differentiation
            % F0: F(t0, X0)
            % dX: Vector of perturbations for each element of X0

            for iX = 1:length(X0)
                dFdX_out(:,iX) = ...
                    obj.getGradient_iX(F, t0, X0, iX, dX(iX), F0, varargin{:});
            end
        end
        
    end
    
    methods (Access = private)
        
        %% AutoDX_iX
        % Performs the actual stepsize-search algorithm on a given element
        % of the input X vector (specified with iX)
        function [dFdX_out, dX_out, dXmax_out, err_out, fcnerr_out, numFCalls] = ...
                AutoDX_iX(obj, F, t0, X0, F0, dX_max, iX, dFdX_known, varargin)
            
            obj.InitArrays(length(F0)); % Initialize arrays as needed
            
            % Extract known elements of dFdX
            if(~isempty(dFdX_known))
                obj.dFdX_computed(:) = dFdX_known(:,iX);
            end
            
            numdX = 0; % Number of computed gradients
            numFCalls = 0; % Number of function calls
            i_dX = obj.np; % This will be incremented into the range [1,np]
            
            % Compute maximum perturbation size based on current
            % parameter's size. We do this by seeing which power of 2
            % approximately equals the parameter's maximum value. NOTE that
            % for higher-order CD methods, we will be computing order*dx,
            % so make sure that the largest perturbation size actually used
            % does not exceed dx_max/order.
            X_size = floor(log2(dX_max(iX)/obj.order));
                        
            for j = 1:obj.numdX_exp
                
                % Exit search if all subfunctions of F have been optimized
                if(all(obj.dFdX_computed))
                    break;
                end
                
                % Quickly compute power-of-2 step size for this iteration
                % The benefit of this is explained in the dissertation
                dX_test = pow2(X_size-j); % Step size for this iteration
                
                % Make sure perturbation is larger than machine precision
                Xplus = X0(iX) + dX_test;
                Xminus = X0(iX) - dX_test;
                if((Xplus - Xminus) ~= 2.0*dX_test)
                    warning('AutoDX:SkipDX', 'Skipping dX=%d because it is too small to reliably affect X', dX_test);
                    continue
                end
                
                i_dX = obj.advance(i_dX, 1); % Increment index
                obj.dX_vec(i_dX) = dX_test; % Store current dX value
                
                % Get gradient using the appropriate FD method
                [obj.dFdX(:,i_dX), obj.Fpert(:,i_dX,:)] = ...
                    obj.getGradient_iX(F, t0, X0, iX, dX_test, F0, varargin{:});
                numFCalls = numFCalls + obj.order;
                
                if(any(~isfinite(obj.dFdX(:,i_dX))))
                    warning('AutoDX:NotFinite', 'Skipping dX=%d because dFdX is not finite', dX_test);
                    i_dX = obj.advance(i_dX, -1); % Reset index
                    continue
                elseif(any(~isreal(obj.dFdX(:,i_dX))))
                    warning('AutoDX:NotReal', 'Skipping dX=%d because dFdX is complex', dX_test);
                    i_dX = obj.advance(i_dX, -1); % Reset index
                end
                
                numdX = numdX + 1; % Increment number of valid gradient points
                
                % Make sure that there are enough points to start approximating
                % truncation error
                if(numdX < obj.npoints)
                    continue
                end
                
                % Shift the 'current' point back (will be reset later). This allows us
                % to work with the "current point, next point" wording instead of
                % "current point, previous point".  It's mainly just to be consistent
                % with the wording of the derived equations.
                i_dX = obj.advance(i_dX, -(obj.npoints-1));
                i_dX_prev = obj.advance(i_dX, -1);
                i_dX_next = obj.advance(i_dX, 1);
                
                % Compute error estimates independently for each element of F(t,X)
                % Note that this (probably) cannot be vectorized since each element of
                % F (each subfunction) may have a different optimal step size.
                for currF = 1:obj.numF
                    if(obj.dFdX_computed(currF))
                        continue;
                    end
                    
                    f1 = obj.dFdX(currF, i_dX);      % Derivative from current point
                    f2 = obj.dFdX(currF, i_dX_next); % Derivative from next point
                    dx1 = obj.dX_vec(i_dX);          % Perturbation from current point
                    dx2 = obj.dX_vec(i_dX_next);     % Purturbation from next point
                    
                    % Compute truncation error for the current point and load into the
                    % error vector. Here we have a special case: if dF/dX is not
                    % changing wrt dX, then the Cn value should be zero. BUT because of
                    % numerical errors, we need to consider "zero" as some max cutoff.
                    Cn = f2 - f1;
                    err_t = 10.0*eps*max(abs(f1), abs(f2)); % Cutoff for Cn
                    if(abs(Cn) <= err_t)
                        Cn = 0; % Change in dFdX is negligible, so Cn should be zero
                    else
                        Cn = Cn/(dx1^obj.order - dx2^obj.order); % Otherwise, compute normally
                    end
                    err_t = Cn * dx1^obj.order; % Estimated truncation error
                    
                    % Check if truncation error is actually O(dx^order). To do so, we
                    % test the slope of log(|err_t|) wrt log(dX), and test whether it
                    % is near the analytically predicted value.
                    if((err_t == 0.0) || (obj.err_t_prev(currF) == 0.0))
                        slope_errt = inf;
                    else
                        % slope = [log(|obj.err_t_prev|) - log(|err_t|)]/[log(dX_prev) - log(dX)]
                        %       = log(|obj.err_t_prev/err_t|) / log(dX_prev/dX)
                        slope_errt = log(abs(obj.err_t_prev(currF)/err_t))/log(obj.dX_vec(i_dX_prev)/dx1);
                        slope_errt = slope_errt - obj.order; % Slope error wrt expected value
                    end
                    obj.err_t_prev(currF) = err_t; % Update truncation error history
                    
                    % Test for a valid nonzero truncation error
                    if(abs(slope_errt) <= obj.order*obj.slope_errt_starttol)
                        obj.num_zero_errt(currF) = 0; % No more nonzero truncation errors
                        
                        % For the first pair of valid dX points, both are valid
                        if(obj.num_valid_errt(currF) == 0)
                            obj.num_valid_errt(currF) = 2;
                        else
                            obj.num_valid_errt(currF) = obj.num_valid_errt(currF) + 1;
                        end
                        
                        % The first dX whose truncation error slope is nearly correct
                        % is the maximum safe dX
                        if(obj.max_valid_dx(currF) == 0.0)
                            obj.max_valid_dx(currF) = obj.dX_vec(i_dX_prev);
                        end
                        
                        % If the current truncation error slope is closer to ideal than
                        % the previous slope, then use the current Cn value as the
                        % correct one. Note that Cn is proportional to the (order+1)-th
                        % derivative of the function.
                        if(abs(slope_errt) < abs(obj.slope_errt_best(currF)))
                            obj.Cn_true(currF) = Cn;
                            obj.slope_errt_best(currF) = slope_errt; % Update best slope with current one
                        end
                        
                        % If we have a sufficient number of consecutive valid points,
                        % then we are on the way to computing an optimal dX value. Note
                        % that "sufficent number" is currently empirically obtained,
                        % and must be at least 2.
                        if(obj.num_valid_errt(currF) == max(2,obj.valid_errt_req - floor(obj.order/2)))
                            obj.errt_valid(currF) = true;
                        end
                        
                    elseif(~obj.errt_valid(currF))
                        % Even if the current truncation error is not valid from a log-log
                        % slope point of view, it may actually be zero. A zero truncation
                        % error corresponds to the function behaving like a polynomial of
                        % degree O(order).
                        obj.num_valid_errt(currF) = 0; % No more valid nonzero truncation errors
                        obj.slope_errt_best(currF) = inf;
                        
                        % Check for zero truncation error
                        if(err_t == 0.0)
                            obj.num_zero_errt(currF) = obj.num_zero_errt(currF) + 1;
                            
                            % First zero truncation error
                            if(obj.num_zero_errt(currF) == 1)
                                obj.max_valid_dx(currF) = dx1;
                                
                                % Truncation error is zero
                            elseif(obj.num_zero_errt(currF) == (obj.valid_errt_req - floor(obj.order/2)))
                                % If the computed dF/dX = 0, then there is no way to
                                % tell if the derivative is just zero at the current X,
                                % or if it is zero at every X (i.e. constant function).
                                % So the max valid dX is set to zero to tell the user
                                % to redo this analysis when X changes.
                                if(abs(f1) < eps)
                                    obj.max_valid_dx(currF) = 0.0;
                                end
                                
                                obj.errt_valid(currF) = true;
                                obj.Cn_true(currF) = 0.0;
                            end
                            
                            % Otherwise truncation error estimates are just plain invalid,
                            % so reset the counters to prepare for possible upcoming valid
                            % truncation errors.
                        else
                            obj.max_valid_dx(currF) = 0.0;
                            obj.num_zero_errt(currF) = 0;
                            obj.Cn_true(currF) = 0.0;
                        end
                    end
                    
                    % If the truncation error is valid, but its slope goes beyond an
                    % allowable tolerance, then we assume that roundoff error has crept
                    % in enough to consider this the 'optimal' dX value.
                    if(obj.errt_valid(currF) && (abs(slope_errt) > obj.order*obj.slope_errt_endtol))
                        
                        % Estimate the current function's condition error assuming the
                        % current dX value minimizes total truncation+roundoff errors.
                        
                        i_dX_opt = i_dX_prev; % Indicate optimal dX index
                        
                        % Compute the corrected optimal step size
                        t_ratio = dx2/dx1;
                        t_ratio = ((1.0 - t_ratio^obj.order)/(1.0 + 1.0/t_ratio))^(1.0/(obj.order+1));
                        hopt_corrected = t_ratio*dx1;
                        %obj.dX_vec(i_dX_opt) = hopt_corrected;
                        
                        % Estimate maximum subtractive cancellation error
                        [err_cond, err_sub] = obj.getErrorCoeff(F0(currF), obj.Fpert(currF,i_dX,:));
                        err_sub = eps*err_sub;
                        
                        % Estimate maximum condition error for this function
                        if(obj.Cn_true(currF) == 0.0)
                            err_c = 0.0; % Cn = 0 means that dFdX is not changing wrt dX, so assume that function has no error
                        else
                            err_c = obj.order*abs(obj.Cn_true(currF))*obj.dX_vec(i_dX_opt)^(obj.order+1) - err_sub;
                        end
                        
                        % Store the estimated condition error for this function
                        if(err_cond == 0.0)
                            fcnerr_out(currF) = 1.0; % Indicate a large relative error
                        else
                            fcnerr_out(currF) = err_c/err_cond; % Convert to relative error
                        end
                        
                        % Store the current dX value as the optimal dX value
                        dX_out(currF) = obj.dX_vec(i_dX_opt);
                        
                        % Store the gradient from the current dX value as the optimal gradient
                        dFdX_out(currF) = obj.dFdX(currF,i_dX_opt); % Uncorrected gradient
                        
                        % Store the max dX value
                        dXmax_out(currF) = obj.max_valid_dx(currF);
                        
                        % Store the estimated roundoff + truncation error
                        err_out(currF) = (err_c + err_sub)/obj.dX_vec(i_dX_opt) ...
                            + abs(obj.Cn_true(currF))*obj.dX_vec(i_dX_opt)^obj.order;
                        
                        if(obj.dFdX(currF,i_dX) ~= 0.0) % Relative error if possible
                            err_out(currF) = err_out(currF)/(obj.dFdX(currF,i_dX_opt) + err_out(currF));
                        end
                        
                        obj.dFdX_computed(currF) = true; % Done analyzing this element of dFdX
                        obj.errt_valid(currF) = false; % Truncation error is no longer valid
                    end
                    
                end
                
                i_dX = obj.advance(i_dX, obj.npoints-1); % Reset the 'current' point
                
            end
            
        end % function AutoDX_iX()
        
        %% advance
        function iout = advance(obj, i, di)
            % Advance integer i by di, wrapping it to the range [1, np].
            % The wrapping direction is consistent with the sign of di. It
            % is assumed that, initially, 1 <= i <= np.
            
            iout = i + mod(di, obj.np);
            
            if(iout > obj.np)
                iout = iout - obj.np;
            end
        end % function advance
        
        %% getGradient_iX
        function [dFdX, Fpert] ...
                = getGradient_iX(obj, F, t0, X0, iX, dX, F0, varargin)
            % Compute the gradient dF/dX with a given order finite-difference method
            
            Xcurr = X0; % Initialize perturbed parameter vector
                        
            if(obj.order == 1) % First-order forward difference method, O(dx)
                Xcurr(iX) = X0(iX) + dX;
                Fplus = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                dFdX = (Fplus - F0)/dX;
                
                if(nargout == 2)
                    Fpert(:,1,1) = Fplus;
                end
                
            elseif(obj.order == 2) % Second-order central difference method, O(dx^2)
                Xcurr(iX) = X0(iX) + dX;
                Fplus = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - dX;
                Fminus = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                dFdX = (Fplus - Fminus)/(2.0*dX);
                
                if(nargout == 2)
                    Fpert(:,1,1) = Fplus;
                    Fpert(:,1,2) = Fminus;
                end
                
            elseif(obj.order == 4) % Fourth-order CD method, O(dx^4)
                Xcurr(iX) = X0(iX) + dX;
                Fplus = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - dX;
                Fminus = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                Xcurr(iX) = X0(iX) + 2.0*dX;
                Fplus2 = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - 2.0*dX;
                Fminus2 = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                dFdX = (Fminus2 - Fplus2 + 8.0*(Fplus - Fminus))/(12.0*dX);
                
                if(nargout == 2)
                    Fpert(:,1,1) = Fplus;
                    Fpert(:,1,2) = Fminus;
                    Fpert(:,1,3) = Fplus2;
                    Fpert(:,1,4) = Fminus2;
                end
                
            elseif(obj.order == 6) % Sixth-order CD method, O(dx^6)
                Xcurr(iX) = X0(iX) + dX;
                Fplus = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - dX;
                Fminus = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                Xcurr(iX) = X0(iX) + 2.0*dX;
                Fplus2 = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - 2.0*dX;
                Fminus2 = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                Xcurr(iX) = X0(iX) + 3.0*dX;
                Fplus3 = F(t0, Xcurr, varargin{:}); % Forward perturbed
                
                Xcurr(iX) = X0(iX) - 3.0*dX;
                Fminus3 = F(t0, Xcurr, varargin{:}); % Backward perturbed
                
                dFdX = (Fplus3 - Fminus3 ...
                    + 9.0*(Fminus2 - Fplus2) ...
                    + 45.0*(Fplus - Fminus)) / (60.0*dX);
                
                if(nargout == 2)
                    Fpert(:,1,1) = Fplus;
                    Fpert(:,1,2) = Fminus;
                    Fpert(:,1,3) = Fplus2;
                    Fpert(:,1,4) = Fminus2;
                    Fpert(:,1,5) = Fplus3;
                    Fpert(:,1,6) = Fminus3;
                end
            end
            
        end % function getGradient_iX
        
        %%
        function [err_cond, err_sub] = getErrorCoeff(obj, Fcurr, Fpert)
            % Computes the coefficients of condition and subtractive cancellation
            % errors in dFdX, given perturbations in the current value of F.
            
            if(obj.order == 1) % O(dx) FD
                Fplus = Fpert(1,1,1);
                err_cond = abs(Fplus) + abs(Fcurr);
                err_sub = max([abs(Fplus), abs(Fcurr)]);
                
            elseif(obj.order == 2) % O(dx^2) CD
                Fplus = Fpert(1,1,1);
                Fminus = Fpert(:,1,2);
                err_cond = (abs(Fplus) + abs(Fminus))/2.0;
                err_sub = max([abs(Fplus), abs(Fminus)])/2.0;
                
            elseif(obj.order == 4) % O(dx^4) CD
                Fplus = Fpert(1,1,1);
                Fminus = Fpert(1,1,2);
                Fplus2 = Fpert(1,1,3);
                Fminus2 = Fpert(1,1,4);
                err_cond = (abs(Fplus2) + abs(Fminus2) ...
                    + 8.0*(abs(Fplus) + abs(Fminus)))/12.0;
                err_sub = (max([abs(Fplus2), abs(Fminus2)]) ...
                    + 8.0*max([abs(Fplus), abs(Fminus)]))/12.0;
                
            elseif(obj.order == 6) % O(dx^6) CD
                Fplus = Fpert(1,1,1);
                Fminus = Fpert(1,1,2);
                Fplus2 = Fpert(1,1,3);
                Fminus2 = Fpert(1,1,4);
                Fplus3 = Fpert(1,1,5);
                Fminus3 = Fpert(1,1,6);
                err_cond = (abs(Fplus3) + abs(Fminus3) ...
                    + 9.0*(abs(Fplus2) + abs(Fminus2)) ...
                    + 45.0*(abs(Fplus) + abs(Fminus)))/60.0;
                err_sub = (max([abs(Fplus3), abs(Fminus3)]) ...
                    + 9.0*max([abs(Fplus2), abs(Fminus2)]) ...
                    + 45.0*max([abs(Fplus), abs(Fminus)]))/60.0;
                
            end
            
        end % function getErrorCoeff
        
        %% Pre-allocate arrays as needed
        function InitArrays(obj, nF)
            
            % Initialize arrays based on number of subfunctions
            if(nF == obj.numF)
                % Arrays sizes have not changed, so just reset their values
                
                obj.err_t_prev(:) = 0.0;
                obj.dX_vec(:) = inf;
                obj.dFdX(:) = 0.0;
                obj.Cn_true(:) = 0.0;
                obj.slope_errt_best(:) = inf;
                obj.num_valid_errt(:) = 0;
                obj.num_zero_errt(:) = 0;
                obj.errt_valid(:) = false;
                obj.max_valid_dx(:) = 0.0;
                obj.dFdX_computed(:) = false;
                obj.Fpert(:) = 0.0;

            elseif(nF > 0)
                % Array sizes have changed, so re-allocate them
                
                obj.numF = nF;
                obj.err_t_prev = zeros(obj.numF, 1);
                obj.dX_vec = inf(obj.np, 1);
                obj.dFdX = zeros(obj.numF, obj.np);
                obj.Cn_true = zeros(obj.numF, 1);
                obj.slope_errt_best = inf(obj.numF, 1);
                obj.num_valid_errt = zeros(obj.numF, 1);
                obj.num_zero_errt = zeros(obj.numF, 1);
                obj.errt_valid = false(obj.numF, 1);
                obj.max_valid_dx = zeros(obj.numF, 1);
                obj.dFdX_computed = false(obj.numF, 1);
                obj.Fpert = zeros(obj.numF, obj.np, obj.order);

            else
                error('AutoDX:InvalidSize', 'Invalid function length %d', nF);
            end
            
            % Dimensions of Fpert also depend on the order
            if(size(obj.Fpert,3) ~= obj.order)
                obj.Fpert = zeros(obj.numF, obj.np, obj.order);
            end
        end
    end
    
end % class AutoDX