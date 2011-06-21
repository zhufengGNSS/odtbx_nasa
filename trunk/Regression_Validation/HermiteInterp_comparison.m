%Comparing the ODTBX Hermite Interpolator to the Matlab interp1 function
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
%   Kevin Berry         09/28/2010   	Original

%% Set the initial state for the satellite
kep1.sma  = 500+JATConstant('meanRadius','Earth')/1000;
kep1.ecc  = 0.37;
kep1.incl = 56;
kep1.raan = 196;
kep1.argp = 344;
kep1.tran = 347;
gm = JATConstant('muEarth')/1e9;
x0 = kep2cart(kep1,gm);

%% Propagate the test orbit for 1 period
T = 2*pi*sqrt(kep1.sma^3/gm);
[tspan,Xref] = integ(@r2bp,[0 T],x0,[],gm);
    
%% truncate the data set
cut_step = 10;
%cut_ind = 1:cut_step:length(tspan); %doesn't includes the last index
cut_ind = union(1:cut_step:length(tspan),length(tspan)); %includes the last index
tcut = tspan(cut_ind);
Xcut = Xref(:,cut_ind);

%% Test interp1 using the Cubic Spline method
Xspline = interp1(tcut,Xcut',tspan,'spline')';

%% Test interp1 using the Piecewise cubic Hermite method
Xpchip = interp1(tcut,Xcut',tspan,'pchip')';

%% Test ddhermite which is a divided difference Hermite method
Xddherm = zeros(size(Xref));
for n=1:length(tspan)
    Xddherm(:,n) = HermiteInterpolator(tcut,tspan(n),6,Xcut')';
end

%% Test ddhermite in java which is a divided difference Hermite method
Xddjava = zeros(size(Xref));
for n=1:length(tspan)
    Xddjava(:,n) = jat.matlabInterface.JATIntegrators.HermiteInterpolator(tcut,tspan(n),6,Xcut')';
end

%% Plot the resulting trajectories
figure
plot3(Xref(1,:),Xref(2,:),Xref(3,:),'k');
hold on
plot3(Xcut(1,:),Xcut(2,:),Xcut(3,:),'r');
plot3(Xspline(1,:),Xspline(2,:),Xspline(3,:),'m');
plot3(Xpchip(1,:),Xpchip(2,:),Xpchip(3,:),'c');
plot3(Xddherm(1,:),Xddherm(2,:),Xddherm(3,:),'b');
plot3(Xddjava(1,:),Xddjava(2,:),Xddjava(3,:),'g');
legend('truth','cut','interp1/spline','interp1/pchip','ddhermite m','ddhermite java');
title('Interpolated Trajectory Results');
xlabel('X (Km)');ylabel('Y (Km)');zlabel('Z (Km)');
hold off

%% Plot the magnitudes of the differences
figure
plot(tspan,sqrt(sum((Xref(1:3,:)-Xspline(1:3,:)).^2)),'m');
hold on
plot(tspan,sqrt(sum((Xref(1:3,:)-Xpchip(1:3,:)).^2)),'c');
plot(tspan,sqrt(sum((Xref(1:3,:)-Xddherm(1:3,:)).^2)),'b');
plot(tspan,sqrt(sum((Xref(1:3,:)-Xddjava(1:3,:)).^2)),'g');
plot(tspan(cut_ind),0,'*r')
legend('interp1/spline','interp1/pchip','ddhermite m','ddhermite java','data points');
title('Magnitude of Errors in Interpolated Trajectories');
xlabel('Time (s)');ylabel('RSS errors (Km)');
hold off

