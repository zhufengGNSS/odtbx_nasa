function [interpolatedState] = HermiteInterpolator(timeArray,currentTime,neqns,newStates)

% HERMITEINTERPOLATOR Interpolates based on modified divided difference.
%
% [interpolatedState] = HermiteInterpolator(timeArray,currentTime,neqns,newStates)
%
%   HermiteInterpolator calculates an interpolated state at the given 
%   currentTime based on the timeArray and newStates.  The neqns describe 
%   the newStates.  This function expects paired state and derivative data 
%   elements, such as [x; y; z; xdot; ydot; zdot].  The number of 
%   dimensions must match half the number of states, i.e. [x; y; zdot; ydot] 
%   requires an neqns of 4.  This function uses ddHermite to interpolate 
%   based on modified divided differences.
%
%   INPUTS 
%   VARIABLE        SIZE        DESCRIPTION (Optional/Default)
%       timeArray   (N*1)       array of times for the interpolator.
%       currentTime (scalar)    user specified output time.
%       neqns       (scalar)    number of states (2, 4 or 6) in newStates
%                               at each time point, (nimimum is 2 for a
%                               single variable and its first derivative)
%       newStates   (N*neqns)   double array of states for the
%                               interpolator for each point in the
%                               timeArray
%
%   OUTPUTS 
%       interpolatedState	(N*1)       array of interpolated states.
%
%   NOTE: Be sure that the units of the timeArray, currentTime, and the
%   newStates values and derivatives are consistent.
%
%   keyword: integrator
%   See also: ddhermite, jatWorldPropagatorRK4, jatWorldPropagatorRK8, 
%   jatForces, createJATWorld
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

%   Adapted from Jat ODToolboxIntegrators.java
%   Author Stephen D Metcalfe, Emergent Space Technologies
%   Date   06/09/2010
%
%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Stephen D Metcalfe  06/09/2010   	Original
%   Stephen D Metcalfe  09/29/2010      Fixed error in selecting start per
%                                       Kevin Berry's suggestion. Added
%                                       check for NaN/Inf return from
%                                       ddhermite.
%   Stephen D Metcalfe  10/04/2010      Fixed MLINs and removed test for
%                                       current time < start time

if nargin < 4 
    error('%s\n%s\n',...
        'Not enough parameters passed',...
        'Expectiing timeArray, currentTime, neqns, newStates');
elseif neqns ~= 2 && neqns ~= 4 && neqns ~= 6
    error('%s %d\n%s\n%s\n',...
        'Not enough equations passed ', neqns,...
        'HermiteInterpolator processes state and derivative data',...
        'Expecting x,y,z,xdot,ydot,zdot data');
elseif length(timeArray) ~= size(newStates)
    error('Length of timeArray does not match length of newStates')
end

interpolatedState = zeros(neqns,1); % allocate returned states array
entries = length(timeArray);        % entries to be sent to ddhermite
s = 1;                              % if 7 entries or less use them all

if entries > 7                      % limit interpolation to 7 entries
	for s=1:entries                 % find the first time AFTER currentTime
		if(currentTime<timeArray(s))% if time in array is after currentTime
			break;                  %   leave start as this index
		end
	end
	if s < 5                        % if currentTime is near the beginning
		s = 1;                      %   use the first 7 entries
    elseif entries-s < 3            % if currentTime is near the end
		s = entries-6;              %   use the last 7 entries
    else                            % if currentTime is in the middle
		s = s-4;                    %   bracket currentTime
	end
	entries = 7;                    % pass 7 entries to ddHermite
end

x  = timeArray(s:s+entries-1,:);      % extract time array
f  = zeros(entries,1);              %#ok<NASGU> % allocate position array
fd = zeros(entries,1);              %#ok<NASGU> % allocate velocity array

deriv = floor(neqns/2);
for g=1:deriv
    f  = newStates(s:s+entries-1,g);        % extract position array
    fd = newStates(s:s+entries-1,g+deriv);  % extract velocity array

	% interpolate using divided difference
	[fo, fdo] = ddhermite(currentTime, x, f, fd, entries);
    if isnan(fo) || isinf(fo) || isnan(fdo) || isinf(fdo)
        error('ddhermite returned NaN or Infinity');
    end
    interpolatedState(g) = fo;
    interpolatedState(g+deriv) = fdo;
end
