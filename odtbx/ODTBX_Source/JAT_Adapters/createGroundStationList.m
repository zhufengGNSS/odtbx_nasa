function y = createGroundStationList(fileName)

% CREATEGROUNDSTATIONLIST  Returns a JAT NASA_GroundStationList object
% containing all the groundstations in the NDOSL file.
%
%  y = createGroundStationList() returns a JAT NASA_GroundStationList object
% containing all the groundstations in the NDOSL file. The information
% about each groundstation can be accessed using the getGroundStationInfo
% function.
%
% The list is stored as a Java object
%
%   INPUTS
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%      fileName  string     Name of NDOSL file to read in (optional)
%
%   OUTPUTS
%      y        (1x1)       jat.groundstations.NASA_GroundStationList object
%
%
%   keyword: JAT Adaptor, groundstation
%   See also GETGROUNDSTATIONINFO
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
%               		(MM/DD/YYYY)
%   Derek Surka          08/27/2007   	Original
%   Derek Surka          09/25/2007     Added filename input
%   Kevin Berry          07/29/2008     Switched call to NASA_GroundStation LIst

if( (nargin<1) || isempty(fileName) )
    f       = filesep;
    jatPath = whereIsJat;
    gsFile  = strcat(jatPath,'groundstations',f,'DBS_NDOSL_WGS84.txt');
else
    gsFile  = which(fileName);
    if isempty(gsFile)
        gsFile = fileName;
    end
end

y = jat.groundstations.NASA_GroundStationList();
y.readFromNDOSLFile(gsFile);

end
