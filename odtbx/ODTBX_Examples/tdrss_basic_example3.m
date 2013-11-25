
tic

%% Input TDRS States as SPEphem type ephemeris files
tdrss.type                     = 'SPEphem';
tdrss.sat(1).filename          = 'TD0_SPephem_2008270.txt'; %TDRS-East
tdrss.sat(2).filename          = 'TD6_SPephem_2008270.txt'; %TDRS-West
tdrss.sat(3).filename          = 'TD4_SPephem_2008270.txt'; %TDRS-Spare
 
%% Input TDRS Antenna Data
tdrss.sat(1).antenna.maxrange  = 70000; %km
tdrss.sat(1).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(2).antenna.maxrange  = 70000; %km
tdrss.sat(2).antenna.halfangle = 13*pi/180; %radians
tdrss.sat(3).antenna.maxrange  = 70000; %km
tdrss.sat(3).antenna.halfangle = 13*pi/180; %radians

%% Set propagator information
epoch  = datenum([2008  9 26  0  0   .000]);
tspan  = 0:10:360;
x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec
options = [];

jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
    % default information for propagating Earth-centric orbits using JAT 
    % force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

%% Measurement Options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', true); 
measOptions = setOdtbxOptions(measOptions,'EarthAtmMaskRadius',6478.12); %km
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);
measOptions = setOdtbxOptions(measOptions,'rSigma',[1e-3 1e-3 1e-3 1e-4 1e-4 1e-4]);

%% Set Ground Station (Optional, internal default is {'WSGT','GTSS'})
% Alternatively you can define the ground station locations in ECEF directly
gsECEF = [-1539.37859986625         -5070.00944053254
         -5160.93110611279          3569.07519673414
          3408.22918525621          1491.68888437059];
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);

%% Run tdrssmeas
[y,H,R] = tdrssmeas_basic(t,x,measOptions);
R;
H;
y;

toc