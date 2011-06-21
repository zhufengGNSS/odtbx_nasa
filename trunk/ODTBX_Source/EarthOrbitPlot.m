function Output = EarthOrbitPlot(Trajectory)
% EARTHORBITPLOT Creates a 3-D plot of an Earth orbit trajectory
%
%  EarthOrbitPlot(Trajectory);
%
%   EarthOrbitPlot plots a three dimentional inertial Earth orbit around
%   a spherical representation of the earth.  The plot can be zoomed and
%   rotated using Matlab's built-in plot utilities.
% 
% INPUT
%      VARIABLE       SIZE       DESCRIPTION (Optional/Default)
%      Trajectory     (3xN)      [x, y, z] inertial trajectory matrix (km)
%
% OUTPUT:
%      Trajectory Plot
%
% VALIDATION TEST
%
%   To perform a validation test, replace the Ephem input with
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
% REGRESSION TEST
%
%   To perform a regression test, replace the Ephem input with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
%
%   keywords: plot, trajectory
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
%               	   (MM/DD/YYYY)
%   Keith Speckman          09/08/2007          Original

% Determine whether this is an actual call to the program or a test

if strcmpi(Trajectory,'ValidationTest') | strcmpi(Trajectory,'RegressionTest')

	Output = EarthOrbitPlot_regression_validation_test();

else

	getPlot(Trajectory);
	Output = 0;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getPlot(xprop);


[X,Y,Z] = sphere(20);

figure
	% Plot Trajectory
	plot3(xprop(1,:),xprop(2,:),xprop(3,:));
	hold on

	% Plot Earth Sphere
	surf(6378*X,6378*Y,6378*Z)
	colormap( [(0) (.7) (1)])

	title('Earth Orbit Trajectory Plot')
	xlabel('X (Km)')
	ylabel('Y (Km)')
	zlabel('Z (Km)')
	axis equal
	hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation and Regression Test
%
% Validation is done by inspection of the plot.
% Regression is confirmed if the script runs without a failure.

function Output = EarthOrbitPlot_regression_validation_test();

% Create an ellipsoid trajectory
a = 6378*4;
b = 6378*3.5;
e = acos(b/a);
h = a*e;
k = 0;
t = [0:2*pi/100:2*pi];
Trajectory(1,:) = h+a*cos(t);
Trajectory(2,:) = k+b*sin(t)*cos(pi/6); 
Trajectory(3,:) = k+b*sin(t)*sin(pi/6); 

Output = EarthOrbitPlot(Trajectory);

end
