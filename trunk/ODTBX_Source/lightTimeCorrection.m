function [t2_lt, x2_lt] = lightTimeCorrection(t1, x1, t2, x2, options, direction, tol, maxn)

% LIGHTTIMECORRECTION Calculates light time corrections including optional environmental delays
%
% [x2_lt, t2_lt] = lightTimeCorrection(x1, t1, x2, t2, options,... 
%                   direction, tol, maxn)
%
% Returns the transmission time and vehicle positions at that time given a
% history of positions and the measurement times. This function calculates
% the light travel time and iterates to find the tracking bodies position
% and velocity at that time. 
%
% x1 is considered the "user" body and its time tags will be held fixed. x2
% is the "tracking" body and its time tags will be varied. Because of this
% method, t1 can be a single value but t2 must be at least 2 for
% interpolation/extrapolation.
%
% The optional input direction specifies whether to move x2 forward in time
% (direction = 1) or backwards in time (direction = -1). For example, if x1
% is a spacecraft and x2 is a ground station, the direction=-1 would move
% x2 backwards to the transmit time (i.e. 1way-forward measurement), and 
% the direction=1 would move x2 forwards to the receive time (i.e.
% 1way-return measurement).
%
% This functionality is originally from Bo Naasz' xlmeas file and was based
% on the GEONS DATSIM Math Spec.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. This function doesn't
% use any options parameters directly, but the ones that are passed into
% the LOSRange function are:
%
%   PARAMETER           VALID VALUES             NOTES
%   epoch                datenum                 Time associated with start
%                                                of simulation
%   useGPSIono          {true, false(default)}   only for GPS sats as x2
%   useIono             {true, false(default)}   only for groundstats as x2
%   useTropo            {true, false(default)}   only for groundstats as x2
%   useChargedParticle  {true, false(default)}   only for groundstats as x2
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t1           (1xN)   Times corresponding to x1 (secs)
%      x1           (6xN)	Positions (km) of user at times t1
%      t2           (1xM)   Times corresponding to x2 (secs)
%      x2           (6xM)   Positions (km) of tracker at times t2
%      options      (1x1)   data structure
%      direction    (1x1)   %1 moves x2 forward, -1 moves x2 backwards
%      tol          (1x1)   light-time correction tolerance (default=1e-10)
%      maxn         (1x1)   maximum iterations (default = 10)
%
%   OUTPUTS
%       t2_lt       (1xN)   Corrected times corresponding to x2_lt
%       x2_lt       (6xN)   Positions (km) of tracker times t2_lt
%
% VALIDATION/REGRESSION TEST
%
%  These have been moved to lightTimeCorrection_test.m in the regression
%  testing framework to conform with the new testing format.
%
%   keyword: measurement
%   See also RRDOTLT, ODTBXOPTIONS, LOSRange
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
%   Bo Naasz            12/03/2005      Original xlmeas.m
%   Derek Surka         07/20/2007      Modified xlmeas to work with OD 
%                                       Toolbox
%   Derek Surka         08/22/2007      Modified input parameters, added 
%                                       functionality
%   Derek Surka         09/06/2007      Added kFixed and propFwd
%   Kevin Berry         05/19/2008      Modified for simplicity
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Ravi Mathur         08/28/2012      Extracted regression test

%% Determine whether this is an actual call to the program or a test

if( nargin < 6 || isempty(direction) )
    direction = -1; %1 moves x2 forward, -1 moves x2 backwards
end
if( nargin < 7 || isempty(tol) )
    tol = 1e-10;    % light time correction tolerance
end
if( nargin < 8 || isempty(maxn) )
    maxn = 10;  % maximum light-time iterations
end

% Its possible to only enter 1 position for first vehicle, but 2nd must
% have more than one for interpolation
if( length(t2)==1 )
    error('LIGHTTIMECORRECTION: Can not vary x2 because it only has 1 input time.')
end

c  = JATConstant('c')/1000;
r1 = x1(1:3,:);
r2 = interp1(t2,x2(1:3,:)',t1,'spline')'; %size(t1) may not equal size(t2)
r  = LOSRange(t1,r1,r2,options);
corr = 1;
n = 0;
while corr>tol % Iterate until r changes by less than light time correction tolerance
    dtlight = r/c;% Light time correction

    % Calculate t2_lt, the corrected time of x2
    t2_lt = t1 + direction*dtlight;

    % Interpolate to get the position of x2 at time t2_lt
    r2 = interp1(t2,x2(1:3,:)',t2_lt,'spline')';

    % Get the new range estimate
    newr = LOSRange(t1,r1,r2,options);
    corr = abs(newr - r);
    r    = newr;
    n    = n+1;
    if n > maxn
        error('XLMEAS: Maximum allowed number of iterations reached')
    end
end

x2_lt = interp1(t2,x2',t2_lt,'spline')';

end