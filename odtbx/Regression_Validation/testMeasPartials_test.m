function failed = testMeasPartials_test
% This function tests the H matrix derivations of gsmeas, tdrssmeas, ddormeas, lnrmeas, and gpsmeas.
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur             08/27/2012      Rename to conform to new
%                                           regression test format

%% Set the random number seed
RandStream.setGlobalStream(RandStream('shr3cong','Seed',546.4737))

%% Test gsmeas H matrix with RangeRate
epoch  = datenum('Jan 1 2010');
gsList = createGroundStationList();
gsID   = {'DS16','DS46','DS66'};
nGS    = length(gsID);
gsECEF = zeros(3,nGS);
for n=1:nGS
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange',true);
measOptions = setOdtbxOptions(measOptions,'useRangerate',true);
measOptions = setOdtbxOptions(measOptions,'useDoppler',false);
measOptions = setOdtbxOptions(measOptions,'useUnit',true);
measOptions = setOdtbxOptions(measOptions,'useAngles',true);
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',-90); %Ignore the horizon
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);

Hpd1 = testMeasPartials(@gsmeas,measOptions);

%% Test gsmeas H matrix with Doppler
epoch  = datenum('Jan 1 2010');
gsList = createGroundStationList();
gsID   = {'DS16','DS46','DS66'};
nGS    = length(gsID);
gsECEF = zeros(3,nGS);
for n=1:nGS
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange',true);
measOptions = setOdtbxOptions(measOptions,'useRangerate',false);
measOptions = setOdtbxOptions(measOptions,'useDoppler',true);
measOptions = setOdtbxOptions(measOptions,'useUnit',true);
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',-90); %Ignore the horizon
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);

Hpd2 = testMeasPartials(@gsmeas,measOptions);

%% Test tdrssmeas H matrix with RangeRate
epoch  = datenum([2008  9 26  0  0   .000]);
tdrss.type  = 'Keplerian';
tdrss.sat(1).epoch = epoch;
tdrss.sat(1).sma  = 42165.3431; %semi-major axis
tdrss.sat(1).ecc  = 0.00026248; %eccentricity
tdrss.sat(1).incl = 1.7350*pi/180; %inclination
tdrss.sat(1).raan = 278.2107*pi/180; %right ascension of the ascending node
tdrss.sat(1).argp = 270.1285*pi/180; %argument of periapse
tdrss.sat(1).mean = 135.9127*pi/180; %mean anomaly
tdrss.sat(2).epoch = epoch;
tdrss.sat(2).sma  = 42166.4487; %semi-major axis
tdrss.sat(2).ecc  = 0.00030059; %eccentricity
tdrss.sat(2).incl = 8.5680*pi/180; %inclination
tdrss.sat(2).raan = 61.9693*pi/180; %right ascension of the ascending node
tdrss.sat(2).argp = 92.0025*pi/180; %argument of periapse
tdrss.sat(2).mean = 37.3454*pi/180; %mean anomaly
tdrss.sat(3).epoch = epoch;
tdrss.sat(3).sma  = 42167.6194; %semi-major axis
tdrss.sat(3).ecc  = 0.00025831; %eccentricity
tdrss.sat(3).incl = 9.8973*pi/180; %inclination
tdrss.sat(3).raan = 54.8622*pi/180; %right ascension of the ascending node
tdrss.sat(3).argp = 138.8066*pi/180; %argument of periapse
tdrss.sat(3).mean = 125.6643*pi/180; %mean anomaly
for n=1:length(tdrss.sat)
    S.M = tdrss.sat(n).mean;
    S.ecc = tdrss.sat(n).ecc;
    tdrss.sat(n).tran = kepanom(S, 'nu'); %true anomaly
end
epoch  = datenum([2010  9 26  0  0   .000]);
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', true);
measOptions = setOdtbxOptions(measOptions,'useDoppler', false);
measOptions = setOdtbxOptions(measOptions, 'EarthAtmMaskRadius', 0); %Ignore the Earth
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);

Hpd3 = testMeasPartials(@tdrssmeas,measOptions);

%% Test tdrssmeas H matrix with Doppler
epoch  = datenum([2008  9 26  0  0   .000]);%datenum('1 Oct 2008 09:27:52.832');
tdrss.type  = 'Keplerian';
tdrss.sat(1).epoch = epoch;
tdrss.sat(1).sma  = 42165.3431; %semi-major axis
tdrss.sat(1).ecc  = 0.00026248; %eccentricity
tdrss.sat(1).incl = 1.7350*pi/180; %inclination
tdrss.sat(1).raan = 278.2107*pi/180; %right ascension of the ascending node
tdrss.sat(1).argp = 270.1285*pi/180; %argument of periapse
tdrss.sat(1).mean = 135.9127*pi/180; %mean anomaly
tdrss.sat(2).epoch = datenum([2008  9 26  0  0   .000]);
tdrss.sat(2).sma  = 42166.4487; %semi-major axis
tdrss.sat(2).ecc  = 0.00030059; %eccentricity
tdrss.sat(2).incl = 8.5680*pi/180; %inclination
tdrss.sat(2).raan = 61.9693*pi/180; %right ascension of the ascending node
tdrss.sat(2).argp = 92.0025*pi/180; %argument of periapse
tdrss.sat(2).mean = 37.3454*pi/180; %mean anomaly
tdrss.sat(3).epoch = datenum([2008  9 26  0  0   .000]);
tdrss.sat(3).sma  = 42167.6194; %semi-major axis
tdrss.sat(3).ecc  = 0.00025831; %eccentricity
tdrss.sat(3).incl = 9.8973*pi/180; %inclination
tdrss.sat(3).raan = 54.8622*pi/180; %right ascension of the ascending node
tdrss.sat(3).argp = 138.8066*pi/180; %argument of periapse
tdrss.sat(3).mean = 125.6643*pi/180; %mean anomaly
for n=1:length(tdrss.sat)
    S.M = tdrss.sat(n).mean;
    S.ecc = tdrss.sat(n).ecc;
    tdrss.sat(n).tran = kepanom(S, 'nu'); %true anomaly
end
epoch  = datenum([2008  9 26  0  0   .000]);%datenum('1 Oct 2008 09:27:52.832');
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', true);
measOptions = setOdtbxOptions(measOptions,'useDoppler', false);
measOptions = setOdtbxOptions(measOptions, 'EarthAtmMaskRadius', 0); %Ignore the Earth
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);

Hpd4 = testMeasPartials(@tdrssmeas,measOptions);

%% Test ddormeas H matrix
epoch    = datenum('01 Jul 2007 14:00:00.000');
gsID = {'DS15','DS45','DS65'};
gsList  = createGroundStationList; %generates a list of all of the ground stations
gsECEF    = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end
ddor.rate = true;
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %datenum format
measOptions = setOdtbxOptions(measOptions,'ddor',ddor);

Hpd5 = testMeasPartials(@ddormeas,measOptions);

%% Test lnrmeas H matrix
epoch  = datenum('1 Dec 2016 00:00:00.000');
relay.type            = 'STKEphem';
relay.sat(1).filename = 'Relay1.e'; %Lunar Relay Satellite #1
relay.sat(2).filename = 'Relay2.e'; %Lunar Relay Satellite #2
relay.MoonMaskRadius  = 0; %Ignore the Moon
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'relay', relay);

Hpd6 = testMeasPartials(@lnrmeas,measOptions);

%% Test gpsmeas H matrix
measOptions = setOdtbxOptions(measOptions, ...
                              'epoch', datenum('Jan 1 2006'), ...
                              'useRange', true, ...
                              'useDoppler', true, ...
                              'PrecnNutnExpire', 0); % Always update precession and nutation

% Parameters specific to link budget
link_budget.AntennaPattern    = {'omni.txt'};
link_budget.AtmosphereMask    = - JATConstant('rEarth','WGS84') / 1000; %removes the Earth
link_budget.RecAcqThresh      = -Inf;            % no threshold
link_budget.RecTrackThresh    = -Inf;            % no threshold
link_budget.TX_AntennaPointing= -1; % 1 for zenith pointing, -1 for nadir pointing
measOptions = setOdtbxOptions(measOptions, 'linkbudget', link_budget);

% Link budget parameters not explicitly set (will be defaulted)
% link_budget.GPSBand           = 'L1';
% link_budget.RXAntennaMask     = truth.rcv_ant_mask;
% link_budget.NoiseTemp         = truth.Ts;                         % K
% link_budget.AtmAttenuation    = truth.Ae;                         % dB
% link_budget.TransPowerLevel   = truth.sv_power;                   % enum
% link_budget.TransPowerOffset  = truth.xmit_power_offset;          % truth units? gpsmeas: dB
% link_budget.GPSBlock          = 6;                                % use a GPS 2D antenna with L1, see gpsmeas
% link_budget.TXAntennaMask     = truth.xmit_ant_mask;              % rad
% link_budget.ReceiverNoise     = truth.Nf;                         % dB
% link_budget.RecConversionLoss = truth.L;                          % dB
% link_budget.SystemLoss        = truth.As;                         % dB
% link_budget.LNAGain           = truth.Ga;                         % dB
% link_budget.CableLoss         = truth.Ac;                         % dB
% link_budget.DynamicTrackRange = truth.dyn_range;                  % dB

% Don't display tons of repeated warning messages
warning('off', 'ODTBX:GPSMEAS:noBodyQuat')
warning('off', 'ODTBX:linkbudget_default:UsingDefault');
Hpd7 = testMeasPartials(@gpsmeas,measOptions);

%% Uncomment this section to reset regression test data
% Hpd1_sav = Hpd1;
% Hpd2_sav = Hpd2;
% Hpd3_sav = Hpd3;
% Hpd4_sav = Hpd4;
% Hpd5_sav = Hpd5;
% Hpd6_sav = Hpd6;
% Hpd7_sav = Hpd7;
% save('testMeasPartials_TestData.mat','Hpd1_sav','Hpd2_sav','Hpd3_sav','Hpd4_sav','Hpd5_sav','Hpd6_sav','Hpd7_sav');

%% Check the results
tol = 0.05;
load('testMeasPartials_TestData.mat')
err1 = max(max(max(abs(Hpd1-Hpd1_sav))));
err2 = max(max(max(abs(Hpd2-Hpd2_sav))));
err3 = max(max(max(abs(Hpd3-Hpd3_sav))));
err4 = max(max(max(abs(Hpd4-Hpd4_sav))));
err5 = max(max(max(abs(Hpd5-Hpd5_sav))));
err6 = max(max(max(abs(Hpd6-Hpd6_sav))));
err7 = max(max(max(abs(Hpd7-Hpd7_sav))));
failed = max([err1,err2,err3,err4,err5,err6,err7])>tol;
