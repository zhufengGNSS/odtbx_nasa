function y = matlabTime2MJD(time)

% MATLABTIME2MJD  Returns the Modified Julian Date of the input time
% specified in Matlab's datenum format.
%
%   y = matlabTime2MJD(time) returns the Modified Julian Date of the input
%   time specified in Matlab's datenum format. This function is used to
%   enable a common time interface with the JAT adapters.
%
%   INPUTS 
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%     time      (1Xn)       Vector of times in datenum format.
%
%   OUTPUTS 
%      y        (1Xn)       Vector of Modified Julian Dates
%
%   keyword: JAT Adapter, time 
%   See also DATENUM, MJD2MATLABTIME, JD2MATLABTIME, MATLABTIME2JD, MJD2JD,
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
%   Derek Surka         06/19/2007   	Original
%   Rob Antonucci       04/12/2010      Simplified computation

% 678942 is the difference between MJD (days since Nov 17. 1858) and
% Matlab's datenum() format (days since Jan 1, 0000).
y = time - 678942;
