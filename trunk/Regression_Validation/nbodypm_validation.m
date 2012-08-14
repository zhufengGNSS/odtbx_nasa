cspice_kclear; %Unloads the kernels

%% Load the default kernels
cspice_furnsh(which('de421.bsp')); %SPK Planetary ephemeris kernel
cspice_furnsh(which('naif0010.tls')); %LSK Leapseconds kernel
cspice_furnsh(which('pck00010.tpc')); %PCK Planetary constants kernel
cspice_furnsh(which('de-403-masses.tpc'))

%% Set Earth Centered user satellite
epoch = datenum('30 DEC 2010 00:00:0.000');
mu_cb = 398600.436;
x0 = [-14816.738456485748
    315.96572152709962
    -6505.1400327072142
    -1.5930781963821792
    -3.5018098894198597
    4.3395209673543702];
tspan = 0:10:86400*30;

%% Test Run #1, Earth Centered with r2bp
disp('r2bp')
tic
[~,x1a] = integ(@r2bp,tspan,x0,[],mu_cb);
toc

disp('nbodypm')
CentralBody = 'Earth';
PointMasses = [];
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB', mu_cb);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
tic
[~,x1b] = integ(@nbodypm,tspan,x0,[],nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x1b(1:3,:)-x1a(1:3,:)).^2)),'b');
title('Comparison of nbodypm with r2pb Earth Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');

%% Test Run #2, Earth Centered with jatForces_km
disp('jatForces_km')

eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'OdeSolver',@ode113,'OdeSolvOpts',...
    odeset('reltol',5e-14,'abstol',1e-18,'initialstep',1,'Vectorized','on'));

jOptions = odtbxOptions('force');
jOptions    = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jOptions    = setOdtbxOptions(jOptions, 'earthGravityModel', '2Body');
jOptions    = setOdtbxOptions(jOptions, 'useSolarGravity', true);
jOptions    = setOdtbxOptions(jOptions, 'useLunarGravity', true);
jOptions    = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', false);
jOptions    = setOdtbxOptions(jOptions, 'useAtmosphericDrag', false);
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
tic
[~,x2a] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);
toc

disp('nbodypm')
CentralBody = 'Earth';
PointMasses = {'SUN','MOON'};
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB', JATConstant('muEarth')*1e-9);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_PM', {JATConstant('muSun')*1e-9,JATConstant('muMoon')*1e-9});
tic
[~,x2b] = integ(@nbodypm,tspan,x0,eOpts,nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x2b(1:3,:)-x2a(1:3,:)).^2)),'b');
title('Comparison of nbodypm with jatForces Earth Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');

%% Test Run #3, Earth Centered with STK
disp('STK')

[~,x3a] = read_stkephem('NbodyCheck1.e');

eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'OdeSolver',@ode113,'OdeSolvOpts',...
    odeset('reltol',5e-14,'abstol',1e-18,'initialstep',1,'Vectorized','on'));

disp('nbodypm')
CentralBody = 'Earth';
PointMasses = {'SUN','MOON','JUPITER BARYCENTER','VENUS BARYCENTER',...
    'SATURN BARYCENTER','MARS BARYCENTER','MERCURY BARYCENTER',...
    'URANUS BARYCENTER','NEPTUNE BARYCENTER'};
GM_PM = {1.327122e11,4.902801076e3,1.267127648383e8,3.24858592079e5...
    3.794058536168e7,4.282837190128e4,2.203209e4,...
    5.794557628118e6,6.836534878892e6};
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB', 398600.4418);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_PM', GM_PM);
tic
[~,x3b] = integ(@nbodypm,tspan,x0,eOpts,nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x3b(1:3,:)-x3a(1:3,:)./1000).^2)),'b');
title('Comparison of nbodypm with STK Earth Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');

%% Set Sun Centered user satellite
mu_cb = 1.327122e11;
x0 = [-43440452.721443251
    926364.66817300871
    -19072100.677102622
    -16.976701153366058
    -37.317164462710828
    46.244258072723632];
tspan = 0:10:86400*30;

%% Test Run #4, Sun Centered with r2bp
disp('r2bp')
tic
[~,x4a] = integ(@r2bp,tspan,x0,[],mu_cb);
toc

disp('nbodypm')
CentralBody = 'Sun';
PointMasses = [];
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB', mu_cb);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
tic
[~,x4b] = integ(@nbodypm,tspan,x0,[],nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x4b(1:3,:)-x4a(1:3,:)).^2)),'b');
title('Comparison of nbodypm with r2pb Sun Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');

%% Test Run #5, Sun Centered with STK
disp('STK')

[~,x5a] = read_stkephem('NbodyCheck2.e');

eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'OdeSolver',@ode113,'OdeSolvOpts',...
    odeset('reltol',5e-14,'abstol',1e-18,'initialstep',1,'Vectorized','on'));

disp('nbodypm')
CentralBody = 'SUN';
PointMasses = {'Earth','MOON','JUPITER BARYCENTER','VENUS BARYCENTER',...
    'SATURN BARYCENTER','MARS BARYCENTER','MERCURY BARYCENTER',...
    'URANUS BARYCENTER','NEPTUNE BARYCENTER'};
GM_PM = {398600.4418,4.902801076e3,1.267127648383e8,3.24858592079e5...
    3.794058536168e7,4.282837190128e4,2.203209e4,...
    5.794557628118e6,6.836534878892e6};
% PointMasses = [];
% GM_PM = [];
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB',1.327122e11);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_PM', GM_PM);
tic
[~,x5b] = integ(@nbodypm,tspan,x0,eOpts,nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x5b(1:3,:)-x5a(1:3,:)./1000).^2)),'b');
title('Comparison of nbodypm with STK Sun Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');

%% Set Jupiter Centered user satellite
mu_cb = 126712767.863;
RJ = 71492;
kep1.sma  = 4*RJ;
kep1.ecc  = 0.37;
kep1.incl = 56*pi/180;
kep1.raan = 196*pi/180;
kep1.argp = 344*pi/180;
kep1.tran = 347*pi/180;
x0 = kep2cart(kep1,mu_cb);
tspan = 0:10:86400*30;

%% Test Run #6, Jupiter Centered with r2bp
disp('r2bp')
tic
[~,x6a] = integ(@r2bp,tspan,x0,[],mu_cb);
toc

disp('nbodypm')
CentralBody = 'Jupiter barycenter';
PointMasses = [];
nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'GM_CB', mu_cb);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
tic
[~,x6b] = integ(@nbodypm,tspan,x0,[],nbodyopt);
toc

% Plot the magnitudes of the differences
figure
plot(tspan./86400,sqrt(sum((x6b(1:3,:)-x6a(1:3,:)).^2)),'b');
title('Comparison of nbodypm with r2pb Jupiter Centered');
xlabel('Time (days)');ylabel('RSS errors (Km)');