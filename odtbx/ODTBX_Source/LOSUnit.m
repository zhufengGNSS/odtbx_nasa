function unitVect = LOSUnit(t, r1, r2, options)
%
% Line of sight unit vector between two objects (provides angle information to the estimator)
%
%   unitVector = LOSUnit(t, r1, r2, options) returns the unit vector
% pointing from one object to another. This model assumes that the
% spacecraft sensor measures offset angles from its boresight, which are
% then converted into a unit line-of-sight vector. The unit vector is then
% rotated from the sensor frame into some external frame of interest, such
% as the inertial frame. This model returns the post-processed unit vector
% representation of the angle data, in the same reference frame as the
% inputs (r1 and r2).
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. THIS IS NOT CURRENTLY
% USED
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t         (1xN)	Times corresponding to r1 (secs)
%      r1        (3xN)	User spacecraft position (km)
%      r2        (3xN)  Tracking spacecraft/ground station position (km)
%      options   (1x1)  Data structure
%
% OUTPUTS
%      unitVect  (3xN)  unit vector representation of the angle data
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, RRDOT, RRDOTLT
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
%   Kevin Berry         09/03/2009      Original

r        = r1 - r2 ;
unitVect = unit(r);