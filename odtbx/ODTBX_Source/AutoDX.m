function [dFdX, dX, dX_maxsafe, err_dFdX, err_F] = ...
    AutoDX(F, t, X, dX_max, iX, order, varargin)
%AutoDX Numerically compute the optimal step size dX to minimize the error
%in the finite-difference derivative dF/dX of the function F(t,X).
%   [DFDX,DX] = AUTODX(F, T, X, DX_MAX, ORDER, IX)
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


end % function AutoDX()

function [dFdX, varargout] ...
    = getGradient_Order(F, t, X, iX, dx, order, FX, varargin)
%Compute the gradient dF/dX with a given order finite-difference method

Xcurr = X; % Initialize perturbed parameter vector

args = varargin;

if(order == 1) % First-order forward difference method, O(dx)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    dFdX = (Fplus - FX)/dX;
    varargout{1} = Fplus;
    varargout{2} = FX;

elseif(order == 2) % Second-order central difference method, O(dx^2)
    Xcurr(iX) = X(iX) + dx;
    Fplus = feval(F, t, Xcurr, args{:}); % Forward perturbed
    
    Xcurr(iX) = X(iX) - dx;
    Fminus = feval(F, t, Xcurr, args{:}); % Backward perturbed
    
    dFdX = (Fplus - Fminus)/(2.0*dx);
    varargout{1} = Fplus;
    varargout{2} = Fminus;
    
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
    varargout{1} = Fplus;
    varargout{2} = Fminus;
    varargout{3} = Fplus2;
    varargout{4} = Fminus2;
    
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
    varargout{1} = Fplus;
    varargout{2} = Fminus;
    varargout{3} = Fplus2;
    varargout{4} = Fminus2;
    varargout{5} = Fplus3;
    varargout{6} = Fminus3;
else
    error('order = %d is not supported', order);
end

end % function getGradient_Order

function err_cond = getCondErrCoeff(order, varargin)
% Computes the coefficient of condition error in F

if(order == 1) % O(dx) FD
    Fplus = varargin{1};
    Fcurr = varargin{2};
    err_accum = abs(Fplus) + abs(Fcurr);
    
elseif(order == 2) % O(dx^2) CD
    Fplus = varargin{1};
    Fminus = varargin{2};
    err_accum = (abs(Fplus) + abs(Fminus))/2.0;
    
elseif(order == 4) % O(dx^4) CD
    Fplus = varargin{1};
    Fminus = varargin{2};    
    Fplus2 = varargin{3};
    Fminus2 = varargin{4};
    err_accum = (abs(Fplus2) + abs(Fminus2) ...
                 + 8.0*(abs(Fplus) + abs(Fminus)))/12.0;
    
elseif(order == 6) % O(dx^6) CD
    Fplus = varargin{1};
    Fminus = varargin{2};    
    Fplus2 = varargin{3};
    Fminus2 = varargin{4};
    Fplus3 = varargin{5};
    Fminus3 = varargin{6};
    err_accum = (abs(Fplus3) + abs(Fminus3) ...
                 + 9.0*(abs(Fplus2) + abs(Fminus2)) ...
                 + 45.0*(abs(Fplus) + abs(Fminus)))/60.0;

else
    error('order = %d is not supported', order);
end

end % function getCondErrCoeff

function err_sub = getSubErrCoeff(order, varargin)
% Computes the coefficient of subtractive cancellation error in dFdX

if(order == 1) % O(dx) FD
    Fplus = varargin{1};
    Fcurr = varargin{2};
    err_accum = max([abs(Fplus), abs(Fcurr)]);
    
elseif(order == 2) % O(dx^2) CD
    Fplus = varargin{1};
    Fminus = varargin{2};
    err_accum = max([abs(Fplus), abs(Fminus)])/2.0;
    
elseif(order == 4) % O(dx^4) CD
    Fplus = varargin{1};
    Fminus = varargin{2};    
    Fplus2 = varargin{3};
    Fminus2 = varargin{4};
    err_accum = (max([abs(Fplus2), abs(Fminus2)]) ...
                 + 8.0*max([abs(Fplus), abs(Fminus)]))/12.0;
    
elseif(order == 6) % O(dx^6) CD
    Fplus = varargin{1};
    Fminus = varargin{2};    
    Fplus2 = varargin{3};
    Fminus2 = varargin{4};
    Fplus3 = varargin{5};
    Fminus3 = varargin{6};
    err_accum = (max([abs(Fplus3), abs(Fminus3)]) ...
                 + 9.0*max([abs(Fplus2), abs(Fminus2)]) ...
                 + 45.0*max([abs(Fplus), abs(Fminus)]))/60.0;

else
    error('order = %d is not supported', order);
end

end % function getSubErrCoeff