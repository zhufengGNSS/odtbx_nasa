function [xDot,A,Q] = jatHCW(t,x,dynarg)

% JATHCW  Returns time derivatives of Hills-Clohessy-Wiltshire equations.
%
%      [xDot,A,Q] = jatHCW(t,x,dynarg) 
%
%   Returns the derivatives of the Hills-Clohessy-Wiltshire equations modeled in JAT.
%
%   INPUTS
%   VARIABLE        SIZE    	DESCRIPTION (Optional/Default)
%       t           (1x1)       Time since start of simulation (secs)
%       x        (1x6 or 6x1)   Input state in CW frame (km and km/s)
%       dynarg   (1x3 or 3x1)   A real array or a cell array of dimension
%                               at least 1 and up to 3
%                               (1) Mean orbit rate (rad/s)
%                               (2) step size (secs) (Required if Q is
%                                   specified in the Output.)
%                               (3) Process noise on velocity state (km/s) 
%                                   (default is 1e-18)
%
%   OUTPUTS
%       xDot        (6x1)       Time derivatives of state (km and seconds)
%       A           (6x6)       State Transition Matrix
%       Q           (6x6)       Process Noise Matrix (km and seconds)
%
%   keyword: JAT Forces
%   See also JATFORCES
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Derek Surka         08/14/2007   	Original
%   Sun Hur-Diaz        04/10/2009      Fixed the dimensions in comments,
%                                       added process noise input option,
%                                       fixed handling of dynarg as cell or
%                                       array

if nargin<3 
    error('Must pass orbit rate to jatHCW');
end

if iscell(dynarg) 
    orbitRate = dynarg{1};
else
    orbitRate = dynarg(1);
end

if length(t) ~= 1
    error('t must be scalar in the current version of jatHCW.');
end

x = x*1000;  % in m

xDot = jat.cm.ClohessyWiltshire.derivs(t,x,orbitRate);

if nargout > 1 
    if length(dynarg) < 2 
        error('Must pass stepsize to jatHCW to compute transition matrix');
    end
    if iscell(dynarg)
        stepSize = dynarg{2};
    else
        stepSize = dynarg(2);
    end
    A = jat.cm.ClohessyWiltshire.phiMatrix(stepSize,orbitRate);
end

if nargout == 3,
    % Note that t is a scalar for now
    if length(dynarg) < 3
        Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2),[1 1 length(t)]); % km
    else
        if iscell(dynarg)
            Qv = dynarg{3};
        else
            Qv = dynarg(3);
        end
        Q = repmat(diag([0 0 0 Qv Qv Qv]),[1 1 length(t)]); % km
    end
end

xDot = xDot/1000;  % back to km

end
