classdef test_tdrssmeas_basic < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        epoch
        t
        x
        measOptions_keplerian
        measOptions_SPEphem
        measOptions_LinkBudget
    end
    
    %% Setup methods run at start of test suite
    methods (TestMethodSetup)
        function create_epoch (testCase)
            testCase.epoch  = datenum([2008  9 26  0  0   .000]); %UTC epoch
        end
        function create_tandx(testCase)
            % Set propagator information
            tspan  = 0:600:3600;
            x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec
            jOptions = odtbxOptions('force');
            jOptions = setOdtbxOptions(jOptions, 'epoch', testCase.epoch); %datenum format
            jatWorld = createJATWorld(jOptions); %creates a java object that stores the
            % default information for propagating Earth-centric orbits using JAT
            % force models.
            eOpts = odtbxOptions('estimator');
            eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
            [testCase.t,testCase.x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);
        end 
        function create_keplerian_measOptions(testCase)
            tdrss.type  = 'Keplerian'; % Input TDRS States as Keplarian Elements
            % TDRS East
            tdrss.sat(1).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
            tdrss.sat(1).sma  = 42165.3431; %semi-major axis
            tdrss.sat(1).ecc  = 0.00026248; %eccentricity
            tdrss.sat(1).incl = 1.7350*pi/180; %inclination
            tdrss.sat(1).raan = 278.2107*pi/180; %right ascension of the ascending node
            tdrss.sat(1).argp = 270.1285*pi/180; %argument of periapse
            tdrss.sat(1).mean = 135.9127*pi/180; %mean anomaly
            % TDRS West
            tdrss.sat(2).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
            tdrss.sat(2).sma  = 42166.4487; %semi-major axis
            tdrss.sat(2).ecc  = 0.00030059; %eccentricity
            tdrss.sat(2).incl = 8.5680*pi/180; %inclination
            tdrss.sat(2).raan = 61.9693*pi/180; %right ascension of the ascending node
            tdrss.sat(2).argp = 92.0025*pi/180; %argument of periapse
            tdrss.sat(2).mean = 37.3454*pi/180; %mean anomaly
            % TDRS Spare
            tdrss.sat(3).epoch = datenum([2008  9 26  0  0   .000]); %UTC epoch
            tdrss.sat(3).sma  = 42167.6194; %semi-major axis
            tdrss.sat(3).ecc  = 0.00025831; %eccentricity
            tdrss.sat(3).incl = 9.8973*pi/180; %inclination
            tdrss.sat(3).raan = 54.8622*pi/180; %right ascension of the ascending node
            tdrss.sat(3).argp = 138.8066*pi/180; %argument of periapse
            tdrss.sat(3).mean = 125.6643*pi/180; %mean anomaly
            % Convert mean anomaly to true anomaly
            for n=1:length(tdrss.sat)
                S.M = tdrss.sat(n).mean;
                S.ecc = tdrss.sat(n).ecc;
                tdrss.sat(n).tran = kepanom(S, 'nu'); %true anomaly
            end

            % Input TDRS Antenna Data
            tdrss.sat(1).antenna.maxrange  = 70000; %km
            tdrss.sat(1).antenna.halfangle = 13*pi/180; %radians
            tdrss.sat(2).antenna.maxrange  = 70000; %km
            tdrss.sat(2).antenna.halfangle = 13*pi/180; %radians
            tdrss.sat(3).antenna.maxrange  = 70000; %km
            tdrss.sat(3).antenna.halfangle = 13*pi/180; %radians

            % Measurement Options
            testCase.measOptions_keplerian = odtbxOptions('measurement');
            testCase.measOptions_keplerian = setOdtbxOptions(testCase.measOptions_keplerian,'epoch',testCase.epoch);
            testCase.measOptions_keplerian = setOdtbxOptions(testCase.measOptions_keplerian,'useRange', true);
            testCase.measOptions_keplerian = setOdtbxOptions(testCase.measOptions_keplerian,'useRangeRate', true); 
            testCase.measOptions_keplerian = setOdtbxOptions(testCase.measOptions_keplerian,'EarthAtmMaskRadius',6478.12); %km
            testCase.measOptions_keplerian = setOdtbxOptions(testCase.measOptions_keplerian,'tdrss', tdrss);

        end
        function create_SPEphem_measOptions(testCase)
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
            
            %% Measurement Options
            testCase.measOptions_SPEphem = odtbxOptions('measurement');
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'epoch',testCase.epoch);
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'useRange', true);
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'useRangeRate', true); 
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'EarthAtmMaskRadius',6478.12); %km
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'tdrss', tdrss);
            testCase.measOptions_SPEphem = setOdtbxOptions(testCase.measOptions_SPEphem,'rSigma',[1e-3 1e-3 1e-3 1e-4 1e-4 1e-4]);
        end
        function create_linkbudget_measOptions(testCase)
            d2r = pi/180;
            %% Link budget information
            % Parameters specific to link budget
            link_budget.AntennaPattern    = {'omni.txt','omni.txt'};           % Specify receive antenna pattern for each antenna
                %  Specify antenna pattern for each antenna, existing antennas are:
                %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
                %     omni.txt                  - zero dB gain,  180 degree half beamwidth
                %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
                %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
                %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
            link_budget.RXAntennaMask     = 90*d2r;     % Cut off angle for the receive antenna
            link_budget.AtmosphereMask    = 50;                      % Mask altitude, km
            link_budget.NoiseTemp         = 300;                    % System noise temp of receive antenna (K)
                % System noise temp [K], space pointing antenna = 290
                % System noise temp [K], earth pointing antenna = 300
            link_budget.AtmAttenuation    = 0.0;                    % Attenuation due to atmosphere, should be negative (dB)
            link_budget.TransPowerLevel   = 3;                      % Transmitter power level (1=min, 2=typical, 3=max)
            link_budget.TransPowerOffset  = 20;                      % Global transmitter power offset (dB)
            link_budget.TXAntennaMask     = 90 * d2r;               % Cut off angle for the transmit antenna (rad)
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
            %link_budget.RXpattern = 'ao40_hga_measured_10db.txt';
            link_budget.TXpattern = 'omni.txt';
            link_budget.TX_AntennaPointing= -1; % 1 for zenith pointing, -1 for nadir pointing
            testCase.measOptions_LinkBudget = testCase.measOptions_keplerian;
            testCase.measOptions_LinkBudget = setOdtbxOptions(testCase.measOptions_LinkBudget, 'linkbudget', link_budget);
        end
    end
    
    %% Teardown methods run at end of test suite
%     methods (TestMethodTeardown)
%         function (testCase)
%         end
%     end
    
    %% Test methods 
    methods (Test)
        function testKeplerianTdrss(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles
            
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,testCase.measOptions_keplerian);
            y_expected = [
          77143.8676829857        -0.539571844028578                       NaN                       NaN          77941.4862100122          1.18932965459402
          78228.5123284356          3.95435946413444                       NaN                       NaN          79843.3003057668          4.83742685697774
          81412.4520066554          6.18418467353267           81376.442501345         -6.55684204133613          83240.3532719398           6.0550782392167
          85144.6557209073          5.92600530129099          77380.2041620424         -6.41305679389542          86685.8437713931          5.17006503654171
                       NaN                       NaN          74065.7701163095         -4.35680062059234                       NaN                       NaN
                       NaN                       NaN          72406.4859646846         -1.04614988349788                       NaN                       NaN
                       NaN                       NaN          72840.7066213739          2.41549480603489                       NaN                       NaN]';

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
        function testSPEphemtdrss(testCase)
            % Run tdrssmeas_basic with SPEEphem defining location of TDRS
            % vehicles
            
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,testCase.measOptions_SPEphem);
            y_expected = [
          77141.7126496068        -0.533986462818289                       NaN                       NaN          77944.5800109052          1.19399493329375
          78227.3862367079          3.95235229822291                       NaN                       NaN          79846.4798758801          4.83348175573453
           81408.605593707          6.17789080804265          81357.0885331791         -6.56222765578761           83239.667852481          6.04700140965134
          85136.7408450822           5.9191799007441          77357.9569469188         -6.41708387938077          86680.1433593953          5.16188406100045
                       NaN                       NaN          74041.8977948749          -4.3579205198466                       NaN                       NaN
                       NaN                       NaN          72383.1456297923         -1.04316879997779                       NaN                       NaN
                       NaN                       NaN           72820.301180623          2.42203147750877                       NaN                       NaN]';

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                 -0.741305988765005         0.670731336064616       -0.0241848266918201                         0                         0                         0
                 -5.95287374988194e-05     -5.80321235194117e-05      0.000215219481825493        -0.741305988765005         0.670731336064616       -0.0241848266918201
                      0.98010431545312         0.171408190433771        -0.100073788188475                         0                         0                         0
                  5.83796581752408e-06      6.46083446774104e-05        0.0001678383045485          0.98010431545312         0.171408190433771        -0.100073788188475
                    -0.657667595110735         0.727901715857909         0.193990789452421                         0                         0                         0
                 -3.28273099178312e-05     -8.55090011694494e-05      0.000209559901373018        -0.657667595110735         0.727901715857909         0.193990789452421];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(diag(getOdtbxOptions(testCase.measOptions_SPEphem,'rSigma').^2),[1 1 7]);
            %R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
        function testDefaultTdrss(testCase)
            % Run tdrssmeas_basic with default location of TDRS vehicles
            
            measOptions = setOdtbxOptions(testCase.measOptions_keplerian,'tdrss', []);
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);
            y_expected = [
          77199.4845235451        -0.337294744524049                       NaN                       NaN           82173.209153457         0.495031732268476
          78394.8321710625          4.10267169509572                       NaN                       NaN          82430.0740878081         0.264940546811722
          81634.5482517139          6.22105420136666          82236.2264321011         -6.31004450040832          82412.1470303643        -0.353819237029461
          85360.1906205278          5.87354284922713          78264.1668908176          -6.5827391249593          82012.5297684694        -0.942729709498263
                       NaN                       NaN          74729.0908131103         -4.90501322110226          81348.3200662233         -1.19971288410771
                       NaN                       NaN          72670.7581625975         -1.78279237792486          80669.2662669526        -0.978704807665696
                       NaN                       NaN          72681.5663061616          1.78957445116117          80274.2230359565        -0.260768889986982]';

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
     -0.740348700979397            0.672223029178649                         0                         0                         0                         0
     -5.57990123396747e-05      -6.1453899239493e-05      0.000217286491366576        -0.740348700979397         0.672223029178649                         0
         0.986251215662958         0.165252956413304                         0                         0                         0                         0
     -1.03799075772032e-05       6.1948643392994e-05       0.00016349135895124         0.986251215662958         0.165252956413304                         0
         0.161946416552939         -0.98679955318477                         0                         0                         0                         0
      7.00827830948725e-05      1.15014802627767e-05      0.000187230214430795         0.161946416552939         -0.98679955318477                         0];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
        function testDefinedGndStation(testCase)
            % Run tdrssmeas_basic with default location of TDRS vehicles
            % and defined ground station locations
            measOptions = setOdtbxOptions(testCase.measOptions_keplerian,'tdrss', []);
            measOptions = setOdtbxOptions(measOptions,'gsID',{'WSGT','GTSS'});
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);
            y_expected = [
          77199.4845235451        -0.337294744524049                       NaN                       NaN           82173.209153457         0.495031732268476
          78394.8321710625          4.10267169509572                       NaN                       NaN          82430.0740878081         0.264940546811722
          81634.5482517139          6.22105420136666          82236.2264321011         -6.31004450040832          82412.1470303643        -0.353819237029461
          85360.1906205278          5.87354284922713          78264.1668908176          -6.5827391249593          82012.5297684694        -0.942729709498263
                       NaN                       NaN          74729.0908131103         -4.90501322110226          81348.3200662233         -1.19971288410771
                       NaN                       NaN          72670.7581625975         -1.78279237792486          80669.2662669526        -0.978704807665696
                       NaN                       NaN          72681.5663061616          1.78957445116117          80274.2230359565        -0.260768889986982]';

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
     -0.740348700979397            0.672223029178649                         0                         0                         0                         0
     -5.57990123396747e-05      -6.1453899239493e-05      0.000217286491366576        -0.740348700979397         0.672223029178649                         0
         0.986251215662958         0.165252956413304                         0                         0                         0                         0
     -1.03799075772032e-05       6.1948643392994e-05       0.00016349135895124         0.986251215662958         0.165252956413304                         0
         0.161946416552939         -0.98679955318477                         0                         0                         0                         0
      7.00827830948725e-05      1.15014802627767e-05      0.000187230214430795         0.161946416552939         -0.98679955318477                         0];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
        function testScheduledTdrss(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles and a tracking schedule for the 
            % TDRS vehicles
            tdrss = getOdtbxOptions(testCase.measOptions_keplerian,'tdrss');
            tdrss.sat(1).schedule = {
                     '26 Sep 2008 00:00:00.000'  '26 Sep 2008 00:20:00.000'
                     '26 Sep 2008 00:24:00.000'  '26 Sep 2008 01:30:00.000'
                     };
            tdrss.sat(2).schedule = {
                     '26 Sep 2008 00:00:00.000'  '26 Sep 2008 00:20:00.000'
                     '26 Sep 2008 00:30:00.000'  '26 Sep 2008 01:30:00.000'
                     };
            tdrss.sat(3).schedule = {
                     '26 Sep 2008 00:00:00.000'  '26 Sep 2008 00:25:00.000'
                     '26 Sep 2008 00:40:00.000'  '26 Sep 2008 01:30:00.000'
                     };
            measOptions = setOdtbxOptions(testCase.measOptions_keplerian,'tdrss', tdrss);
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);
            y_expected = 1.0e4*[
           7.714386768298576   7.822851232843555                 NaN   8.514465572090732                 NaN                 NaN                 NaN
          -0.000053957184403   0.000395435946413                 NaN   0.000592600530129                 NaN                 NaN                 NaN
                         NaN                 NaN                 NaN                 NaN   7.406577011630946   7.240648596468461   7.284070662137386
                         NaN                 NaN                 NaN                 NaN  -0.000435680062059  -0.000104614988350   0.000241549480603
           7.794148621001217   7.984330030576678   8.324035327193982                 NaN                 NaN                 NaN                 NaN
           0.000118932965459   0.000483742685698   0.000605507823922                 NaN                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
        end
        function testBadrSigma(testCase)
            % Make sure tdrssmeas_basic throws a warning when given bad
            % value for rSigma and gives expected outputs
            measOptions = setOdtbxOptions(testCase.measOptions_keplerian,'rSigma',[1e-3 1e-3 1e-3 1e-4 1e-4]);
            % Run tdrssmeas_basic
            [y,H,R]=testCase.verifyWarning(@()tdrssmeas_basic(testCase.t,testCase.x,measOptions),...
                'tdrssmeas_basic:invalid_rSigma');
            y_expected = [
          77143.8676829857        -0.539571844028578                       NaN                       NaN          77941.4862100122          1.18932965459402
          78228.5123284356          3.95435946413444                       NaN                       NaN          79843.3003057668          4.83742685697774
          81412.4520066554          6.18418467353267           81376.442501345         -6.55684204133613          83240.3532719398           6.0550782392167
          85144.6557209073          5.92600530129099          77380.2041620424         -6.41305679389542          86685.8437713931          5.17006503654171
                       NaN                       NaN          74065.7701163095         -4.35680062059234                       NaN                       NaN
                       NaN                       NaN          72406.4859646846         -1.04614988349788                       NaN                       NaN
                       NaN                       NaN          72840.7066213739          2.41549480603489                       NaN                       NaN]';

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end

        function testBadTdrsType(testCase)
            % Make sure tdrssmeas_hifi throws an error when given invalid
            % tdrss type
            measOptions = testCase.measOptions_keplerian;
            measOptions.tdrss.type = [];
            % Run tdrssmeas_basic
            testCase.verifyError(@()tdrssmeas_basic(testCase.t,testCase.x,measOptions),...
                'tdrssmeas_basic:invalid_tdrss_type');
        end
        
        function testKeplerianTdrss_IonoModel_basic(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles
            measOptions = setOdtbxOptions(testCase.measOptions_keplerian,'useIonosphere', true);%@IonoModel_hifi); %
            
            % Run tdrssmeas_basic
            warning('OFF','IonoModel_basic:r2_not_ground_station')
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);
            warning('ON','IonoModel_basic:r2_not_ground_station')
                        
            y_expected = 1.0e4*[
               7.714387156088839   7.822851620349047   8.141245587900293   8.514465959069305                 NaN                 NaN                 NaN
              -0.000053957184889   0.000395435945950   0.000618418466914   0.000592600529714                 NaN                 NaN                 NaN
                             NaN                 NaN   8.137645984733497   7.738022198036385   7.406578839159222   7.240650468109646   7.284072576263919
                             NaN                 NaN  -0.000655684124160  -0.000641305601938  -0.000435679987208  -0.000104614916175   0.000241549550030
               7.794149033040943   7.984330443024411   8.324035739946940   8.668584790093176                 NaN                 NaN                 NaN
               0.000118932966224   0.000483742686292   0.000605507824343   0.000517006503901                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
      
        function testKeplerianTdrss_TropoModel_basic(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles
            measOptions = testCase.measOptions_keplerian;
            measOptions = setOdtbxOptions(measOptions,'useTroposphere', true);
            
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);

            y_expected =  1.0e+04 *[
               7.714387785554491   7.822852247657278   8.141246213168793   8.514466582417715                 NaN                 NaN                 NaN
              -0.000053957192763   0.000395435938488   0.000618418459871   0.000592600523098                 NaN                 NaN                 NaN
                             NaN                 NaN   8.137644586296044   7.738020752831365   7.406577348737398   7.240648934068559   7.284071000246082
                             NaN                 NaN  -0.000655684202604  -0.000641305677815  -0.000435680060438  -0.000104614986680   0.000241549482324
               7.794149907364663   7.984331322894367   8.324036624008851   8.668585676934677                 NaN                 NaN                 NaN
               0.000118932987646   0.000483742703140   0.000605507836401   0.000517006510992                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            
        end
         
        function testKeplerianTdrss_ChargedParticleModel_basic(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles
            measOptions = testCase.measOptions_keplerian;
            measOptions = setOdtbxOptions(measOptions,'useChargedParticle', true);
            
            % Run tdrssmeas_basic
            [y,H,R] = tdrssmeas_basic(testCase.t,testCase.x,measOptions);
            y_expected = 1.0e4*[
               7.714386768633569   7.822851233178548   8.141245201000530   8.514465572425731                 NaN                 NaN                 NaN
              -0.000053957184403   0.000395435946413   0.000618418467353   0.000592600530129                 NaN                 NaN                 NaN
                             NaN                 NaN   8.137644250469531   7.738020416539276   7.406577011965985   7.240648596803502   7.284070662472429
                             NaN                 NaN  -0.000655684204134  -0.000641305679390  -0.000435680062059  -0.000104614988350   0.000241549480603
               7.794148621336222   7.984330030911686   8.324035327528991   8.668584377474321                 NaN                 NaN                 NaN
               0.000118932965459   0.000483742685698   0.000605507823922   0.000517006503654                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
        end
        
        function testKeplerianTdrss_withLinkBudget(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles and a LinkBudget
            
            % Run tdrssmeas_basic
            [y,H,R,AntLB] = testCase.verifyWarning(@()tdrssmeas_basic(testCase.t,testCase.x,testCase.measOptions_LinkBudget),...
                'ODTBX:GSMEAS:noBodyQuat');
            y_expected = 1.0e+04 *[
               7.714386768298576   7.822851232843555   8.141245200665534                 NaN                 NaN                 NaN                 NaN
              -0.000053957184403   0.000395435946413   0.000618418467353                 NaN                 NaN                 NaN                 NaN
                             NaN                 NaN   8.137644250134501   7.738020416204241   7.406577011630946   7.240648596468461   7.284070662137386
                             NaN                 NaN  -0.000655684204134  -0.000641305679390  -0.000435680062059  -0.000104614988350   0.000241549480603
               7.794148621001217   7.984330030576678   8.324035327193982                 NaN                 NaN                 NaN                 NaN
               0.000118932965459   0.000483742685698   0.000605507823922                 NaN                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.Halpha_r
            Halpha_r_expected = [
           0.733520633104795   0.988903348096852   1.512136150067561   1.998411347070086   2.406639359534008   2.726048381027634   2.811370606754448
           2.940683259181947   2.420913571388128   1.869292174355140   1.349965574218727   0.876371305200500   0.556215091235198   0.665561016499735
           0.850658147651562   1.216304116596038   1.716558697225260   2.162853587863256   2.520328628613372   2.734667162894562   2.662318050650024];
            ABSTOL = 10*eps(max(max(abs(Halpha_r_expected))));
            testCase.verifyEqual(AntLB{2}.Halpha_r,Halpha_r_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HCN0
            HCN0_expected = 1.0e2*[
           0.336853700703305   0.334324574181007   0.327311927781234  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000   0.330217472555737   0.337796659077158   0.341862432137444   0.340798923008644
           0.335663578231637   0.331358622920084  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HCN0_expected))));
            testCase.verifyEqual(AntLB{2}.HCN0,HCN0_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HAd
            HAd_expected = 1.0e2*[
          -1.871434173824728  -1.873963300347027  -1.880975946746799  -1.888533197239103  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -1.886405366875405  -1.878070401972296  -1.870491215450875  -1.866425442390590  -1.867488951519390
          -1.872624296296397  -1.876929251607950  -1.884136453899395  -1.890885091797588  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAd_expected))));
            testCase.verifyEqual(AntLB{2}.HAd,HAd_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HAr
            HAr_expected = [
             0     0     0  -100  -100  -100  -100
          -100  -100  -100     0     0     0     0
             0     0  -100  -100  -100  -100  -100];
            ABSTOL = 10*eps(max(max(abs(HAr_expected))));
            testCase.verifyEqual(AntLB{2}.HAr,HAr_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HAP
            HAP_expected = 1.0e2*[
          -1.671434173824728  -1.673963300347027  -1.680975946746800  -1.688533197239103  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -1.686405366875405  -1.678070401972296  -1.670491215450875  -1.666425442390590  -1.667488951519390
          -1.672624296296397  -1.676929251607950  -1.684136453899395  -1.690885091797588  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAP_expected))));
            testCase.verifyEqual(AntLB{2}.HAP,HAP_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HRP
            HRP_expected = 1.0e2*[
          -1.671434173824728  -1.673963300347027  -1.680975946746800  -2.688533197239103  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -2.686405366875405  -1.678070401972296  -1.670491215450875  -1.666425442390590  -1.667488951519390
          -1.672624296296397  -1.676929251607950  -2.684136453899395  -2.690885091797588  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HRP_expected))));
            testCase.verifyEqual(AntLB{2}.HRP,HRP_expected,'AbsTol',ABSTOL);
            % Check AntLB{2}.HAt
            HAt_expected = [
             0     0     0  -300  -300  -300  -300
          -300  -300  -300     0     0     0     0
             0     0  -300  -300  -300  -300  -300];
            ABSTOL = 10*eps(max(max(abs(HAt_expected))));
            testCase.verifyEqual(AntLB{2}.HAt,HAt_expected,'AbsTol',ABSTOL);
        end
        
        function testKeplerianTdrss_withLinkBudget_with_2D_antennas_and_Quat(testCase)
            % Run tdrssmeas_basic with keplarian orbital elements defining
            % location of TDRS vehicles and a LinkBudget
            qatt = repmat([0;0;0;1],1,length(testCase.t));
            measOptions = testCase.measOptions_LinkBudget;
            link_budget = getOdtbxOptions(measOptions,'linkbudget');
            link_budget.AntennaPattern    = {'sensysmeas_ant_2D.txt'};
            link_budget.TXpattern = 'sensysmeas_ant_2D.txt';
            measOptions = setOdtbxOptions(measOptions, 'linkbudget', link_budget);
            % Run tdrssmeas_basic
            [y,H,R,AntLB] = tdrssmeas_basic(testCase.t,testCase.x,measOptions,qatt);
            
            y_expected = 1.0e+04 *[
               7.714386768298576                 NaN                 NaN                 NaN                 NaN                 NaN                 NaN
              -0.000053957184403                 NaN                 NaN                 NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN                 NaN   7.240648596468461   7.284070662137386
                             NaN                 NaN                 NaN                 NaN                 NaN  -0.000104614988350   0.000241549480603
                             NaN                 NaN                 NaN                 NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN                 NaN                 NaN                 NaN];

            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            % Check H
            % only check first part of H to save space and time
            H_expected = [
          -0.742821985720926   0.669022891051028  -0.024977365340341                   0                   0                   0
          -0.000059532435171  -0.000058064975708   0.000215202993524  -0.742821985720926   0.669022891051028  -0.024977365340341
           0.979885503831094   0.173074629499780  -0.099345719611854                   0                   0                   0
           0.000005608658853   0.000064592868335   0.000167850515620   0.979885503831094   0.173074629499780  -0.099345719611854
          -0.659488549547777   0.726425021515655   0.193343065899321                   0                   0                   0
          -0.000032733403000  -0.000085507916948   0.000209615926767  -0.659488549547777   0.726425021515655   0.193343065899321];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(6),[1 1 7]);
            ABSTOL = 10*eps(max(max(max(R_expected))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.Halpha_r
            Halpha_r_expected = [
           1.545816363625434   1.663389550265595   1.722010269386996   1.714665325433029   1.664709976656622   1.593747260971031   1.516111857396414
           1.471286460484761   1.566241778171114   1.632430754077231   1.650447450600562   1.611057912565839   1.523495631509910   1.420277267361348
           1.765364695346714   1.873774207589260   1.916324552583821   1.895652107140453   1.836916886643970   1.761240567281440   1.682724817661689];
            ABSTOL = 10*eps(max(max(abs(Halpha_r_expected))));
            testCase.verifyEqual(AntLB{1}.Halpha_r,Halpha_r_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HCN0
            HCN0_expected = 1.0e2*[
           0.326418153819899  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000   0.332116885770487   0.349707407671871
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HCN0_expected))));
            testCase.verifyEqual(AntLB{1}.HCN0,HCN0_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAd
            HAd_expected = 1.0e2*[
          -1.871434173824728  -1.873963300347027  -1.880975946746799  -1.888533197239103  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -1.886405366875405  -1.878070401972296  -1.870491215450875  -1.866425442390590  -1.867488951519390
          -1.872624296296397  -1.876929251607950  -1.884136453899395  -1.890885091797588  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAd_expected))));
            testCase.verifyEqual(AntLB{1}.HAd,HAd_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAr
            HAr_expected = 1.0e2*[
          -0.050447379551392  -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000  -0.049104341214295
          -0.041683747494141  -0.050098801792906  -1.000000000000000  -1.000000000000000  -1.000000000000000  -0.049758145081386  -0.031096599060640
          -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAr_expected))));
            testCase.verifyEqual(AntLB{1}.HAr,HAr_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAP
            HAP_expected = 1.0e2*[
          -1.631422341156743  -1.633972065220034  -1.641021355369793  -1.648574806073459  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -1.646445742852872  -1.638117933077341  -1.630512530903683  -1.626412843676161  -1.627483867795522
          -1.632618336618197  -1.636957927195456  -1.644181006378969  -1.650911733141066  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAP_expected))));
            testCase.verifyEqual(AntLB{1}.HAP,HAP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HRP
            HRP_expected = 1.0e2*[
          -1.681869720708135  -2.686254267995532  -2.697500898678496  -2.704526975302987  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -2.697040259136452  -2.689611678074793  -2.680441370534671  -1.676170988757547  -1.658580466856162
          -2.691748642175225  -2.696964775282579  -2.704315601526389  -2.710939810233895  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HRP_expected))));
            testCase.verifyEqual(AntLB{1}.HRP,HRP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAt
            HAt_expected = 1.0e2*[
           0.040011832667986  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000   0.040012598714429   0.040005083723868
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAt_expected))));
            testCase.verifyEqual(AntLB{1}.HAt,HAt_expected,'AbsTol',ABSTOL);
        end
    end
    
end

