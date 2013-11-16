%GSMEAS_EXAMPLE This demonstrates the use of the gsmeas measurement model.
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

close all

d2r             = pi/180;

t = [0 60 120 180];
%x = [18000 10000 10000 0 0 0]';

%Over Goldstone
x = [5158.712338     -1514.921648      3777.724544       2.754779       9.375862      -0.001963;
5310.829717      -949.016866      3768.064835       2.313986       9.479618      -0.319662;
5436.235143      -378.334492      3739.451984       1.865419       9.535073      -0.633026;
5534.648919       194.234583      3692.270615       1.415291       9.542710      -0.937961]';

%% Set up Ground Station information
epoch  = datenum('Jan 1 2006');
gsList = createGroundStationList();
gsID   = {  'DS12'
            'DS16' 
            'DS46'
            'DS66' };
nGS    = length(gsID);
gsECEF = zeros(3,nGS);
for n=1:nGS
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end

%% Set up normal ODTBX options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'frequencyTransmit', 2106406.250);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',0);
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
measOptions = setOdtbxOptions(measOptions,'gsID', gsID);

%% Link budget information
% Parameters specific to link budget
link_budget.AntennaPattern    = {'omni.txt'};           % Specify receive antenna pattern for each antenna
    %  Specify antenna pattern for each antenna, existing antennas are:
    %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
    %     omni.txt                  - zero dB gain,  180 degree half beamwidth
    %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
    %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
    %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
link_budget.RXAntennaMask     = 180*d2r;     % Cut off angle for the receive antenna
link_budget.AtmosphereMask    = 0;                      % Mask altitude, km
link_budget.NoiseTemp         = 300;                    % System noise temp of receive antenna (K)
    % System noise temp [K], space pointing antenna = 290
    % System noise temp [K], earth pointing antenna = 300
link_budget.AtmAttenuation    = 0.0;                    % Attenuation due to atmosphere, should be negative (dB)
link_budget.TransPowerLevel   = 2;                      % Transmitter power level (1=min, 2=typical, 3=max)
link_budget.TransPowerOffset  = 0;                      % Global transmitter power offset (dB)
link_budget.TXAntennaMask     = 70 * d2r;               % Cut off angle for the transmit antenna (rad)
    %  The actual mask used is the lesser of this mask and the limit of the defined pattern
    %  Note:  mask = 70 deg includes entire defined pattern
    %         mask = 42 deg includes only main and first side lobes
    %         mask = 26 deg includes only main lobe
link_budget.ReceiverNoise     = -3;                     % Noise figure of receiver/LNA (dB)
link_budget.RecConversionLoss = 0;                      % Receiver implementation, A/D conversion losses (dB)
link_budget.SystemLoss        = 0;                      % System losses, in front of LNA (dB)
link_budget.LNAGain           = 40;                     % Gain provided by the LNA (dB)
link_budget.CableLoss         = -2;                     % Cable losses after LNA (dB)
link_budget.RecAcqThresh      = 32;                     % Receiver acquisition threshold (dB-Hz)
link_budget.RecTrackThresh    = 32;                     % Receiver tracking threshold (dB-Hz)
link_budget.DynamicTrackRange = 15;                     % Receiver dynamic range (dB)
link_budget.Frequency = 146.520e6;                      % Hz 
%link_budget.RXpattern = 'ao40_hga_measured_10db.txt';
link_budget.TXpattern = 'ao40_hga_measured_10db.txt';
link_budget.TX_AntennaPointing= 1; % 1 for zenith pointing, -1 for nadir pointing
measOptions = setOdtbxOptions(measOptions, 'linkbudget', link_budget);

%% Set variations and run gsmeas
measOptions = setOdtbxOptions(measOptions,'useLightTime',true);
measOptions = setOdtbxOptions(measOptions,'useIonosphere',false);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);
[y1,H1,R1] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useIonosphere',true);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',true);
[y2,H2,R2] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useLightTime',false);
measOptions = setOdtbxOptions(measOptions,'useIonosphere',false);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);
[y3,H3,R3] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useIonosphere',true);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',true);
[y4,H4,R4] = gsmeas(t, x, measOptions);

%% Display outputs
disp(' ')
disp(' ')
disp([y1, y2, y3, y4])
