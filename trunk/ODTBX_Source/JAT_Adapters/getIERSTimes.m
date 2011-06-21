function [a, b] = getIERSTimes()

% GETIERSTIMES Returns the applicable times for the JAT IERS Parameters.
%
% IERS is the International Earth Rotation and Reference Service.
%
% [a, b] = getIERSTimes() returns the earliest and latest valid times
% used by JAT.  The times are in Modified Julian Date UTC.  This data
% comes from jat.spacetime.FitIERS.  If the JAT data is uninitialized then
% zeros are returned.  A FitIERS instance must be created in the JVM for
% this data to be populated.
%
% Note: Times used by JAT outside of the time interval formed by a and b
% will not have accurate IERS support for polar motion or UT1-UTC 
% differences.
%
%   INPUTS
%   (none)
%
%   OUTPUTS
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%      a        (1X1)       Earliest valid IERS time point (MJD, UTC), or
%                           zero if not initialized.
%      b        (1X1)       Latest valid IERS time point (MJD, UTC), or
%                           zero if not initialized.
%
%   keyword: JAT Time IERS polar motion
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

a = jat.spacetime.FitIERS.getEarliestTime();
b = jat.spacetime.FitIERS.getLatestTime();