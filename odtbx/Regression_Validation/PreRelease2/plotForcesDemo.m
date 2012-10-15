function plotForcesDemo(filename1, filename2)

% Compare the results of two Pre-Release R2 JAT Forces Demos
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

% Assumes each file contains
%   tPlot
%   xPlot = [x; y; z; xDot; yDot; zDot]

% Load first file
load(filename1);
t1 = tPlot;
r1 = xPlot(1:3,:);
v1 = xPlot(4:6,:);
dc = dcm('ric',r1,v1);
x1RIC = zeros(3,length(t1));
for( i=1:size(dc,3) )
    x1RIC(:,i) = dc(:,:,i)*r1(:,i);
end
clear tPlot xPlot

% Load second file
load(filename2)
t2 = tPlot;
r2 = xPlot(1:3,:);
v2 = xPlot(4:6,:);
% dc = dcm('ric',r2,v2);
x2RIC = zeros(3,length(t2));
for( i=1:size(dc,3) )
    x2RIC(:,i) = dc(:,:,i)*r2(:,i);
end
clear tPlot xPlot

rDiff = r1 - r2;
vDiff = v1 - v2;
xDiff = x1RIC - x2RIC;

% Compute magnitudes of position vectors
m     = length(t1);  % Assume length(t1) = length(t2)
r1Mag = zeros(1,m);
r2Mag = zeros(1,m);
rDiffMag = zeros(1,m);

for i=1:m
    r1Mag(i) = sqrt( r1(1,i)^2 + r1(2,i)^2 + r1(3,i)^2 );
    r2Mag(i) = sqrt( r2(1,i)^2 + r2(2,i)^2 + r2(3,i)^2 );
    rDiffMag(i) = sqrt( rDiff(1,i)^2 + rDiff(2,i)^2 + rDiff(3,i)^2 );
end

labels = ['R';'I';'C'];
for i=1:3
    titlestr = [labels(i), ' Position'];
    figure('Name',titlestr);
    subplot(2,1,1)
    plot(t1,x1RIC(i,:),t2,x2RIC(i,:));
    ylabel([labels(i), ' (km)']);
    title(titlestr);
    
    subplot(2,1,2)
    plot(t1,xDiff(i,:)*1000);
    ylabel([labels(i), ' (m)']);
    xlabel('Time (secs)');
    title([labels(i), ' Difference'])
end


labels = ['x';'y';'z'];
for i=1:3
    titlestr = [labels(i), ' Position'];
    figure('Name',titlestr);
    subplot(2,1,1)
    plot(t1,r1(i,:),t2,r2(i,:));
    ylabel([labels(i), ' (km)']);
    title(titlestr);
    
    subplot(2,1,2)
    plot(t1,rDiff(i,:)*1000);
    ylabel([labels(i), ' (m)']);
    xlabel('Time (secs)');
    title([labels(i), ' Difference'])

    titlestr = [labels(i), ' Velocity'];
    figure('Name',titlestr);
    subplot(2,1,1)
    plot(t1,v1(i,:),t2,v2(i,:));
    ylabel([labels(i), ' (km/s)']);
    title(titlestr);
    
    subplot(2,1,2)
    plot(t1,vDiff(i,:)*1000);
    ylabel([labels(i), ' (m/s)']);
    xlabel('Time (secs)');
    title([labels(i), ' Difference'])
end

figure('Name','Position Magnitude');
subplot(2,1,1)
plot(t1,r1Mag,t2,r2Mag);
ylabel(['Magnitude (km)']);
title('Magnitudes of Position Vectors')

subplot(2,1,2)
plot(t1,(r1Mag-r2Mag)*1000);
xlabel('Time (secs)');
ylabel('Difference (m)');
title('Difference (|r1|-|r2|)')

figure('Name','r Difference');
plot(t1,rDiffMag*1000);
xlabel('Time (secs)');
ylabel('Difference (m)');
title('Difference (|r1-r2|)')

