function y = matlabTime2JD(time)

% MATLABTIME2JD  Returns the Julian Date of the input datenum
%
%   y = matlabTime2JD(time) returns the Julian Date of the input datenum.
%   This function is used to enable a common time interface with the JAT
%   adapters.
%
%   INPUTS 
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%     time      (1Xn)       Vector of times in datenum format
%
%   OUTPUTS 
%      y        (1Xn)       Vector of Julian dates
%
%   keyword: JAT Time 
%   See also DATENUM, MATLABTIME2MJD, MJD2MATLABTIME, JD2MATLABTIME, JD2MJD,
%   MJD2JD
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
%   Derek Surka         06/19/2007   	Original
%   Brent Wm. Barbee    04/16/2009      Added a for loop to allow for
%                                       processing of an array of input
%                                       values.
%   Rob Antonucci       04/12/2010      Simplified computation

% 2400000.5 is the difference between MJD and JD (as MJD is defined)
% 678942 is the difference between MJD (days since Nov 17. 1858) and
% Matlab's datenum() format (days since Jan 1, 0000).
y = time - 678942 + 2400000.5;


