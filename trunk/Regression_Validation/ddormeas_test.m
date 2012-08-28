function failed = ddormeas_test()
%
% ddormeas_test Regression test for ddormeas
% See also: ddormeas.m
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
%   Ravi Mathur         08/28/2012      Extracted from ddormeas.m

disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
epoch    = datenum('01 Jul 2007 14:00:00.000');
x0 = [1.5e8 0 0 0 5e-2 0]';% km & km/sec
P0 = diag([10 10 10 .1 .1 .1].^2); %#ok<NASGU> %initial covariance matrix (km^2 & km^2/sec^2)

% Set the Tracking Stations
gsID = {'DS15',... Goldstone 34 Meter
        'DS45',... Canberra 34 Meter
        'DS65',... Madrid 34 Meter
        }; %The NASA ID for each ground station
gsSigma = [1e-6 1e-9... (km range & km/s rangerate sigma) Goldstone 34 Meter
           1e-6 1e-9... (km range & km/s rangerate sigma) Canberra 34 Meter
           1e-6 1e-9... (km range & km/s rangerate sigma) Madrid 34 Meter
           ]; %The measurement sigmas for each ground station

% Set the measurement options
gsList  = createGroundStationList; %generates a list of all of the ground stations
gsECEF    = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end
ddor.rate = true;

% Initialize options structure. In this case it is a measurement structure
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %datenum format
measOptions = setOdtbxOptions(measOptions,'rSigma',gsSigma);
measOptions = setOdtbxOptions(measOptions,'ddor',ddor);

% Set propagator information
tspan  = 0:6000:86400;
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
% default information for propagating Earth-centric orbits using JAT
% force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

% Set Expected Results
y_expected = [
                         0    -1.34971696639199e-005     9.09494701772928e-013     6.16669884865372e-006                         0     2.40550487673035e-007
       -0.0881998013637713    -1.37598456273745e-005      -0.00842719790489355    -8.31648270929288e-006        0.0229141913841886      5.5283473699666e-006
        -0.143038607860944    -3.28935742790178e-006       -0.0812128905763529    -1.37676275811608e-005        0.0389013225230883    -2.47178109346298e-006
        -0.123861529403257       9.023266792148e-006        -0.146281489874127    -5.92356581078859e-006       -0.0262979671406356    -1.99124986330633e-005
       -0.0527340503595042     1.25798347026248e-005        -0.137660174973007     9.01487714432817e-006        -0.192994662347701    -3.38448750284676e-005
      0.000338087622367311     3.11895617408728e-006       -0.0486858088784174     1.87682090613836e-005        -0.399891655461943    -3.18757302564991e-005
       -0.0305411990620996    -1.37288855254595e-005        0.0577744859347149     1.39978103044686e-005        -0.538331014643518    -1.17829225805727e-005
        -0.155215144309295    -2.60527629277344e-005        0.0910142582843037    -4.53809954024803e-006        -0.525460187659519     1.60023425420065e-005
        -0.313688684326962    -2.39666815885454e-005      -0.00258448906970443    -2.59160925819613e-005        -0.364185146097952     3.50332746428538e-005
        -0.411790444298276    -6.73245689288389e-006        -0.197421004382704    -3.62513974665533e-005        -0.146098953617184     3.40167136666419e-005
         -0.38242396581245     1.63858949784293e-005        -0.399902802635552    -2.81583801182888e-005       0.00637499022741395     1.46100838588964e-005
        -0.230134243601242     3.21928975547009e-005        -0.506190530161803    -5.81548299080938e-006        0.0191982840624405    -9.71786149971974e-006
       -0.0282301305951478     3.23080950674613e-005        -0.466922345807461     1.79141247614341e-005       -0.0875235588455325    -2.29954392734217e-005
         0.129034686206069     1.84616987743615e-005        -0.314673366203351     3.02275675136154e-005        -0.218558604502505    -1.76625191324119e-005
         0.185424678748859     7.02347231531437e-007        -0.136478502226964     2.67860760906881e-005        -0.273003076374152     5.89128784256221e-007]';

% Display Expected Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Expected DDOR Measurements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-5s %14s %14s %14s %14s %14s %14s\n','Time','gs1 - gs2','gs1 - gs2','gs1 - gs3','gs1 - gs3','gs2 - gs3','gs2 - gs3');
fprintf('%3s %13s %15s %13s %15s %13s %15s\n','(s)','(km)','(km/s)','(km)','(km/s)','(km)','(km/s)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:length(y_expected)
    fprintf('%-6i %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n',t(n),y_expected(:,n))
end
disp(' ')
disp(' ')

% Run tdrssmeas
y = ddormeas(t,x,measOptions);

% Display Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Calculated DDOR Measurements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-5s %14s %14s %14s %14s %14s %14s\n','Time','gs1 - gs2','gs1 - gs2','gs1 - gs3','gs1 - gs3','gs2 - gs3','gs2 - gs3');
fprintf('%3s %13s %15s %13s %15s %13s %15s\n','(s)','(km)','(km/s)','(km)','(km/s)','(km)','(km/s)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:length(y)
    fprintf('%-6i %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n',t(n),y(:,n))
end

passed = tol > max( max( abs( y - y_expected ) ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end

end