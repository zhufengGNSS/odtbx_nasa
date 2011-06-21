% Validation test data for LLA2ecef.m
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
%   Author                     Date            Comment
%                          (MM/DD/YYYY)
%   Keith Speckman          06/03/2008          Original

Re = JATConstant('rEarth')/1e3;

LatLonAltValues = [	0	0	0
			0	90	0
			0	180	0
			0	-90	0
			0	45	0
			0	135	0
			0	-135	0
			0	-45	0
			90	0	0
			-90	0	0
			45	90	0
			-45	-90	0
			0	0	0.1*Re
			0	-90	0.1*Re
			90	0	0.1*Re	
			20	35	25];

for k = 1:size(LatLonAltValues,1)
	x(k,:) = [LLA2ecef(LatLonAltValues(k,1),LatLonAltValues(k,2),LatLonAltValues(k,3))]';
end
