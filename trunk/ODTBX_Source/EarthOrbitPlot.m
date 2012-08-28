function EarthOrbitPlot(trajectory)
% EARTHORBITPLOT Creates a 3-D plot of an Earth orbit trajectory
%
%  EarthOrbitPlot(Trajectory);
%
%   EarthOrbitPlot plots a three dimentional inertial Earth orbit around
%   a spherical representation of the earth.  The plot can be zoomed and
%   rotated using Matlab's built-in plot utilities.
%
% Note: running this file will change the background color of the figure in
% which it plots to black.  Although other figures are unaffected, this
% figure will continue to have a black background, even if you close it,
% for the remainder of the Matlab session.  To change it back to white, use
% the command "whitebg."
% 
% INPUT
%      VARIABLE       SIZE       DESCRIPTION (Optional/Default)
%      Trajectory     (3xN)      [x, y, z] inertial trajectory matrix (km)
%
% OUTPUT:
%      Trajectory Plot
%
% VALIDATION/REGRESSION TEST
%
%   These tests have been moved to EarthOrbitPlot_test.m to conform to
%   the new regression testing format.
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
%   Ravi Mathur         08/28/2012      Extracted regression test

figure

% Plot trajectory
plot3(trajectory(1,:),trajectory(2,:),trajectory(3,:));
hold on

% Plot Earth Sphere
earth(6378)

title('Earth Orbit Trajectory Plot')
xlabel('X (Km)')
ylabel('Y (Km)')
zlabel('Z (Km)')
whitebg(gcf,'k')
hold off

end
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

end