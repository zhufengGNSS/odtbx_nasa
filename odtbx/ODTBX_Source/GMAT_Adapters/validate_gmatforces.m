function validate_gmatforces(filename)

% TEST SCRIPT!
%
% Script to analytically validate gmatForces_km.
%
% This script loads GMAT and exercises the derivative model using calls to 
% ODTBX's integ integrator.

% GMAT: General Mission Analysis Tool
%
% Copyright (c) 2002-2011 United States Government as represented by the
% Administrator of The National Aeronautics and Space Administration.
% All Other Rights Reserved.
%
% Developed jointly by NASA/GSFC and Thinking Systems, Inc. under NASA
% Prime Contract NNG10CP02C, Task Order 28.
%
% Author: Darrel J. Conway, Thinking Systems, Inc.
%         Adapted from script written by Allen Brown, 
%            Emergent Space Technologies, Inc.
% Created: 2011/05/17


if (nargin >= 1)
    if strcmp(filename, '') == 1
       filename = 'GmatConfig.script';
    end;
else
    filename = 'GmatConfig.script';
end

if CheckForFile(filename) == 0
    errmsg = sprintf('%s%s%s\n','The GMAT configuration script "',...
        filename,'" does not exist, so GMAT cannot be loaded');
    error(errmsg);
end

startgmat(filename);

mu = 398600.4415
r = [10000; 0; 0];
v = [0; sqrt(mu/r(1)); 0];
state = [r; v]
T = sqrt((4*pi*pi*r(1)^3)/mu)
t = 1;

w = 17;

% note, integ makes a copy of w and uses that copy
[tPlot,xPlot] = integ(@gmatforces_km,[0 T],state,[],w);

figure;
subplot(3,1,1);plot(tPlot,xPlot(1,:)); ylabel('X axis (km)');
subplot(3,1,2);plot(tPlot,xPlot(2,:)); ylabel('Y axis (km)');
subplot(3,1,3);plot(tPlot,xPlot(3,:)); ylabel('Z axis (km)');
xlabel('Time (sec)');
disp('state numerical integration difference at analytical period for circular orbit:');
diff = xPlot(:,end)-xPlot(:,1)

disp('Starting timing tests...');
tlen = [1e2 5e2 1e3 5e3 1e4 5e4 1e5 5e5];
for t = 1:length(tlen)
    tic;
    [tPlot,xPlot] = integ(@gmatforces_km,[0:tlen(t)]*(T/tlen(t)),state,[],w);
    tElapsed(t) = toc;
end
figure;
subplot(2,1,1);loglog(tlen,tElapsed,'bx-');
xlabel('points');
ylabel('Elapsed time (sec)');
subplot(2,1,2);plot(tlen,tElapsed./tlen,'bx-');
xlabel('points');
ylabel('Elapsed time/point');

closegmat;
