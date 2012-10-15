function y = mJD2JD(mjd)

% MJD2JD  Converts Modified Julian Date to Julian Date.
%
%   y = mJD2JD(mjd) converts the input Modified Julian Date to Julian Date.
%   
%   INPUTS 
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%      mjd      (1X1)       Modified Julian Date
%
%   OUTPUTS 
%      y        (1X1)       Julian Date 
%
%   keyword: JAT Adapter, time 
%   See also MATLABTIME2MJD, MJD2MATLABTIME, JD2MATLABTIME, MATLABTIME2JD, 
%   JD2MJD
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
%   Derek Surka         05/18/2007   	Original
%   Rob Antonucci       04/12/2010      Simplified computation

% 2400000.5 is the difference between MJD and JD (as MJD is defined)
y = mjd + 2400000.5;

end
