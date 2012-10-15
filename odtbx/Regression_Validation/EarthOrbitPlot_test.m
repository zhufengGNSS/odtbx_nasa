function failed = EarthOrbitPlot_test
%
% EarthOrbitPlot_test Regression test for EarthOrbitPlot
% See also: EarthOrbitPlot.m
%
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
%
%   Ravi Mathur         08/28/2012      Extracted from EarthOrbitPlot.m

% Create an ellipsoid trajectory
a = 6378*4;
b = 6378*3.5;
e = acos(b/a);
h = a*e;
k = 0;
t = 0:2*pi/100:2*pi;
Trajectory(1,:) = h+a*cos(t);
Trajectory(2,:) = k+b*sin(t)*cos(pi/6); 
Trajectory(3,:) = k+b*sin(t)*sin(pi/6); 

EarthOrbitPlot(Trajectory);
failed = 0;

end