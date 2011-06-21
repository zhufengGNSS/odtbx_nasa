function [fail] = jatIonoDelayModel_test(ploton)
% Validation test for jatIonoDelayModel.  This test compares with an 
% independent formulation of jatIonoDelayModel.
%
% Inputs:
%   ploton  (optional) any non-empty value causes the comparison plot to be
%                      displayed (not usually provided when running as
%                      part of a test suite)
%
% Outputs:
%   fail    1=test failure, 0=test success
%
% Equations and data generated from Campbell's thesis
% https://fftbwall.gsfc.nasa.gov/odtbx/trac.cgi/attachment/wiki/IncrementSix/WillCampbell.zip
% equations and data: S. Hur-Diaz
% regression test adaptation: A. Brown
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

% check ploton input, set to 1 for plotting, 0 for no plotting
if exist('ploton','var') && ~isempty(ploton)
    ploton = 1;
else
    ploton = 0;
end

fail = 0; % default to not failed

%
% Comparison data from Campbell thesis:
%
% Ground station data
lam = 0; % longitue of ground station in deg
phi = 0; % latitude of ground station in deg
f = 7650e6; % frequency (Hz)

% Satellite data
el = 90; % elevation of satellite in deg
az = 0;  % azimuth of satellite in deg

d2r=pi/180;

rE=6378; % Earth radius in km
h = 350; % mean ionospheric height in km
T=32; % hours
thd=14; % hours
N = 3.192e18; % Hz^2-m
D = 3.192e19; % Hz^2-m

%
% equations coded from Campbell thesis:
%
zeta=asin(rE*cos(el*d2r)/(rE+h)); % rad

phip=asin(sin(phi*d2r)*sin(el*d2r+zeta)+cos(phi*d2r)*cos(az*d2r)*cos(el*d2r+zeta)); % rad
lamp = lam+asin((sin(az*d2r)*cos(el*d2r+zeta))/cos(phip))/d2r; % deg

st=0:1:24;
drho=zeros(length(st),1);
for i=1:length(st)
    UT=st(i);
    tp = mod(lamp/15 + UT,24); % hours
    chi = min(pi/2, 2*pi*abs(tp-thd)/T);

    drho(i)=(N+D*cos(chi))/f^2/cos(zeta);
end

if ploton
    figure
    plot(st,drho,'xb-'),grid
    title('Iono Delay vs. Station local time')
    xlabel('hrs'),ylabel('meters')
end

%
% Setup case to run the above using jatIonoDelayModel.m
%
ephem.Lat = 0;
ephem.Lon = 0;
ephem.Height = 0;
ephem.StationInfo = 'LatLonHeight';
ephem.SatCoords = 'ECEF';
ephem.SatPos = [1 0 0]'*7000e3;
ephem.SignalFreq = 7650e6;

drhoiono=zeros(length(st),1);
for i=1:length(st)
    
    ephem.Epoch = datenum(2000,9,19,st(i),0,0);
    drhoiono(i) = jatIonoDelayModel(ephem);
    
end

if ploton
    hold on
    plot(st,drhoiono,'or--'),grid
    title('Iono Delay from jatIonoDelayModel vs. UTC')
    xlabel('hrs'),ylabel('meters')
    hold off
end

delta=drho-drhoiono;
maxerr=max(abs(delta));

% Test criteria: less than 1 mm difference in calculated delay:
if maxerr >= 0.001
    fail = 1;
    disp(sprintf('jatIonoDelayModel_test failed with a maximum error of %f, which is greater than 1mm.',maxerr));
end

end
