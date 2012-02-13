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
%   Keith Speckman     09/08/2007       Original
%   Russell Carpenter  02/13/2012       Various Improvements

% Determine whether this is an actual call to the program or a test

if strcmpi(Trajectory,'ValidationTest') || strcmpi(Trajectory,'RegressionTest')
	Output = EarthOrbitPlot_regression_validation_test();
else
	getPlot(Trajectory);
	Output = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getPlot(xprop)

figure
% Plot Trajectory
plot3(xprop(1,:),xprop(2,:),xprop(3,:));
hold on

% Plot Earth Sphere
earth(6378)

title('Earth Orbit Trajectory Plot')
xlabel('X (Km)')
ylabel('Y (Km)')
zlabel('Z (Km)')
whitebg('k')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function earth(varargin)
% Sphere with Earth map superimposed
% Optional arguments
% - Earth radius
% - Globe transparency level, 1 = opaque, through 0 = invisible
% - Number of globe panels around the equator deg/panel = 360/npanels

% Original by Khashayar Parsay
% Modified by Russell Carpenter

if nargin < 1 || isempty(varargin{1})
    Re = 1;
else
    Re = varargin{1};
end
if nargin < 2 || isempty(varargin{2})
    alpha = 1;
else
    alpha = varargin{2};
end
if nargin < 3 || isempty(varargin{3})
    npanels = 72;
else
    npanels = varargin{3};
end   
[x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata = imread('earth.jpg');
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, ...
    'EdgeColor', 'none');
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation and Regression Test
%
% Validation is done by inspection of the plot.
% Regression is confirmed if the script runs without a failure.

function Output = EarthOrbitPlot_regression_validation_test()

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

Output = EarthOrbitPlot(Trajectory);