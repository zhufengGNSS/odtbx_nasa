function [y,H,R] = tdrssmeasCon(t,X,options)
% TDRSSMEAS  Makes TDRSS based measurements.
%
%   [y,H,R] = TDRSSMEAS(t,x,options) creates TDRSS measurements
% based on the information in OPTIONS. The measurement types that can be
% returned are range, range-rate, and doppler. See the OD Toolbox functions
% named LOSRANGE, LOSRANGERATE, and LOSDOPPLER for details of each
% measurement type.
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   epoch                 datenum              UTC time associated with the
%                                              start of the simulation
%   gsElevationConstraint degs                 Spacecraft elevation must be
%                                              above this value to return
%                                              measurement
%   useRange            {true(default), false} Return range measurement
%   rangeType           {'1way','1wayFWD','1wayRTN','2way'(default)}
%                                              1way = 1wayFWD
%   useRangeRate        {true(default), false} Return rangerate measurement
%   useDoppler          {true, false(default)} Return doppler measurement
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only used for Tropo and
%                                              for Doppler
%   rSigma              {(1xM),ones(1,M)*1e-3) Measurement covariance
%   useLightTime        {true, false(default)} Include light time delays
%   useIonosphere       {true, false(default)} Includes Ionospheric delays
%   useTroposphere      {true, false(default)} Includes Tropospheric delays
%   EarthAtmMaskRadius  {scalar>0, 6478.12}    km, Effective Earth
%                                              Occultation Radius
%   trdss               structure(km, radians) Structure of all TDRS data
%       Keplerian Version:
%           tdrss.type         = 'Keplerian';
%           tdrss.sat(i).epoch = datenum([2008  9 26  0  0   .000]); %UTC
%           tdrss.sat(i).sma   = 42165.3431; %semi-major axis in km
%           tdrss.sat(i).ecc   = 0.00026248; %eccentricity
%           tdrss.sat(i).incl  = 1.7350*pi/180; %inclination in radians
%           tdrss.sat(i).raan  = 278.2107*pi/180; %right ascension of the
%                                                 ascending node in radians
%           tdrss.sat(i).argp  = 270.1285*pi/180; %argument of periapse in
%                                                 radians
%           tdrss.sat(i).tran  = 135.9336*pi/180; %true anomaly in radians
%       Ephem Version:
%           tdrss.type            = 'SPEphem';
%           tdrss.sat(1).filename = 'TD0_SPephem_2008270.txt'; %TDRS-East
%           tdrss.sat(2).filename = 'TD6_SPephem_2008270.txt'; %TDRS-West
%           tdrss.sat(3).filename = 'TD4_SPephem_2008270.txt'; %TDRS-Spare
%       Optional Antenna Data:
%           tdrss.sat(i).antenna.maxrange = 70000; %max dist to user in km
%           tdrss.sat(i).antenna.halfangle = 13*pi/180; %Nadir facing field
%                                                       of view in radians
%       Optional Schedule Data:
%           tdrss.sat(i).schedule = {
%                     '30 Jun 2015 23:27:36.155'  '1 Jul 2015 00:03:23.077'
%                     '01 Jul 2015 00:46:46.000'  '1 Jul 2015 00:48:04.709'
%                     };
%
% The options parameters associated with each of the desired measurement
% types are passed to the appropriate function.
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each tdrs are grouped together in rows.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs)
%      x            (6xN)   spacecraft state [position;velocity] (km)
%      options      (1x1)   data structure
%
%   OUTPUTS
%      y            (MxN)   measurements
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%
% VALIDATION TEST
%
%  To perform a validation test, pass in 'ValidationTest' as the
%  only input argument and specify only one output argument.
%
% REGRESSION TEST
%
%  To perform a regression test, pass in 'RegressionTest' as the
%  only input argument and specify only one output argument.
%
%   keyword: measurement
%   See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, ODTBXOPTIONS, GSMEAS,
%   READ_SPEPHEM

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Kevin Berry         12/19/2008      Original
%   Kevin Berry         12/22/2008      Replaced elevation.m with
%                                       jatStaAzEl.m
%   Kevin Berry         12/29/2008      Corrected a mistake in the
%                                       light time section
%   Kevin Berry         05/28/2009      Added Earth occultation
%                                       Added antenna constraints
%   Kevin Berry         06/25/2009      Added time scale comments
%   Brandon Farzad      08/06/2009      Added Scheduling capabilities
%   Kevin Berry         02/03/2010      Changed to faster gsmeas options

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    y = tdrssmeas_test();
else
    [y,H,R] = Get_tdrssmeas(t,X,options);
end
end


%% Main function
function [y,H,R] = Get_tdrssmeas(t,X,options)


%% Get values from options
epoch = getOdtbxOptions(options, 'epoch', NaN); %UTC
uselt = getOdtbxOptions(options, 'useLightTime', false);
tdrss = getOdtbxOptions(options, 'tdrss', []);
RMin = getOdtbxOptions(options, 'EarthAtmMaskRadius', 100+JATConstant('meanRadius','Earth')/1000);
%Default is the Mean Radius of the Earth plus 100 km of Atmosphere
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
numtypes     = useRange + useRangeRate + useDoppler;
ephemer   = getOdtbxOptions(options, 'EphemerisError',[]); % ephem has const,growth,Ei

if isnan(epoch); error('An epoch must be set in the options structure.'); end
if isempty(tdrss); error('tdrss state information must be set in the options structure'); end
if size(t,1)==length(t), t=t'; end

x=X(1:6,:);
if isfield(ephemer,'const') && size(X,1)==7 % assume ephem only appears in truth options...
    xf=X(7,1); % PULL FIRST ELEMENT BECAUSE IT SHOULD BE CONSTANT!!!!!!
%     if ~isfield(ephem,'phi')
%         ephem.phi=rand*2*pi;
%         ephem.theta=rand*2*pi;
%     end
    omega=0.7272205217e-4;
    eflag=true;
    hsize=7;
else
    eflag=false;
    hsize=6;
end
theta=cell(size(tdrss.sat));
dx=cell(size(tdrss.sat));

%% Get TDRS state
if strcmpi(tdrss.type, 'keplerian')
    if uselt
        c       = JATConstant('c')/1000;
        ltDT    = sqrt(sum(x(1:3,:).^2))/c; %
        tdrss.t = union(t-ltDT, union(t, t+ltDT)); %seconds from scenario epoch
    else
        tdrss.t = t; %seconds from scenario epoch
    end
    for n=1:length(tdrss.sat);
        %Set the time span for each satellite
        tdrss.sat(n).tspan = tdrss.t;
        epdiff = (epoch-tdrss.sat(n).epoch)*86400; %secs from tdrss epoch to scenario epoch

        GM         = JATConstant('muEarth')/10^9;
        KOEf = kepprop2b(tdrss.sat(n),tdrss.sat(n).tspan+epdiff,GM);
%         theta{n}=KOEf.tran;
        tdrss.x{n} = kep2cart(KOEf,GM);
    end

elseif strcmpi(tdrss.type, 'spephem')
    for n=1:length(tdrss.sat) % More accurate
        [tdrss.t tdrss.x{n}] = read_spephem(tdrss.sat(n).filename);
        tdrss.sat(n).epoch   = tdrss.t(1); %UTC
        epdiff               = (epoch-tdrss.sat(n).epoch)*86400; %secs from tdrss epoch to scenario epoch
        tdrss.sat(n).tspan   = (tdrss.t - tdrss.t(1))*86400 - epdiff;
        KOEf = kepprop2b(tdrss.sat(n),tdrss.sat(n).tspan+epdiff,GM);
%         theta{n} = KOEf.tran;
    end
else
    error([tdrss.type ' is not a supported state type'])
end

if eflag
%     thetaC=[0 2*pi/3 4*pi/3];
    for n=1:length(tdrss.sat)
        if ephemer.const
%             theta{n}(:)=thetaC(n);
%             E=epherr(ephemer.Ei,omega,theta{n},ephemer.growth);
            zt=zeros(size(t));
            E=epherr(ephemer.Ei,omega,zt,ephemer.growth,ephemer.theta(n),ephemer.phi(n));
        else
%             E=epherr(ephemer.Ei,omega,theta{n},ephemer.growth);
            E=epherr(ephemer.Ei,omega,t,ephemer.growth,ephemer.theta(n),ephemer.phi(n));
        end
        [~,dx{n}] = rot2cart(tdrss.sat(n).tspan,tdrss.x{n},[],E); % RIC -> XYZ
        dx{n}=dx{n}{:};
        tdrss.x{n} = tdrss.x{n} + xf*dx{n};
    end
end

%% Loop over all TDRS Satellites
%initialize Y and H and R
y    = nan(length(tdrss.sat)*numtypes,length(t));
H    = zeros(length(tdrss.sat)*numtypes,hsize,length(t));
R = [];

% gsIDs    = {'WSGT','GTSS'};
% gsList   = getOdtbxOptions(options, 'gsList', []);
% if isempty(gsList); gsList = createGroundStationList(); end;
% gsECEF = zeros(3,length(gsIDs));
% for n=1:length(gsIDs)
%     gsECEF(:,n) = getGroundStationInfo(gsList,gsIDs{n},'ecefPosition',epoch);
% end
gsECEF = [-1539.37859986625         -5070.00944053254
         -5160.93110611279          3569.07519673414
          3408.22918525621          1491.68888437059];
options  = setOdtbxOptions(options,'gsElevationConstraint', 0);
for n=1:length(tdrss.x)
    %% Set the ground station
    Ephem.satPos      = tdrss.x{n}(1:3,1)*1000; %ECI satellite coordinates (m)
    Ephem.SatCoords   = 'ECI';
    Ephem.Epoch       = epoch+tdrss.sat(n).tspan(1)/86400;
    Ephem.StationInfo = 'ECEF';
    Ephem.staPos      = gsECEF*1000;
    [~,Elv]  = jatStaAzEl(Ephem);
    gsECEFn  = gsECEF(:, max(Elv) == Elv );
    options  = setOdtbxOptions(options,'gsECEF',gsECEFn);
    %% Put Tdrss states into an array
    for j=size(tdrss.x,2):-1:1
        xt(:,:,j)=tdrss.x{1,j}(:,:);
    end

    %% set the tracking schedule
    % Check schedule for times when tracking is done
    if isfield(tdrss.sat(n),'schedule')
        tSch = NaN(size(tdrss.sat(n).schedule,1),2);%Peallocated for speed
        %this is where the schedule is put into a format that is easier to
        %work with
        for j=1:size(tdrss.sat(n).schedule,1)
            tSch (j,:) = ((datenum(tdrss.sat(n).schedule(j,:))-epoch)*86400);
        end
        tind = [];
        for m=1:size(tSch,1)
            tind = union(tind, find( tSch(m,1)<=t & t<=tSch(m,2) ));
        end
    else
        tind = 1:length(t);
    end
    % t, x, and xt are restricted only to the times when tracking is done
    t1   = t(tind);
    x1   = x(:,tind);
    x2  = xt(:,tind,n);
    if ~isempty (t1)
        %% Get the measurements for the Ground and Space leg
        if uselt
            [ys, Hs, ~, t2_lt, x2_lt] = rrdotlt(t1, x1, tdrss.sat(n).tspan, tdrss.x{n},...
                options); %calculate the Space leg
            if length(t2_lt)==2 %then it was a 2way measurement
                options       = setOdtbxOptions(options,'rangeType','1wayFWD');
                [yg1,Hg1]      = gsmeas(t2_lt{1}, x2_lt{1}, options); %calculate the FWD Ground leg
                options       = setOdtbxOptions(options,'rangeType','1wayRTN');
                [yg2,Hg2]     = gsmeas(t2_lt{2}, x2_lt{2}, options); %calculate the RTN Ground leg
                yg            = (yg1+yg2)/2; %average the ground legs
                Hg            = (Hg1+Hg2)/2; % average the ground legs
            else %else it was a 1way measurement and gsmeas will know which way
                [yg, Hg] = gsmeas(t2_lt{1}, x2_lt{1}, options); %calculate the Ground legl
            end
        else
            if isequal(tdrss.sat(n).tspan,t1)
                x2 = tdrss.x{n};
            else
                x2 = interp1(tdrss.sat(n).tspan,tdrss.x{n}',t1,'spline')'; %interpolate the TDRS position to time t1
            end
            [ys, Hs] = rrdot(t1,x1,x2,options); %calculate the Space leg
            [yg, Hg] = gsmeas(t1,x2,options); %calculate the Ground leg
        end

        %% Apply the antenna constraints
        [u1 l1] = unit(x1(1:3,:)-x2(1:3,:)); %unit vector and distance from tdrss to user
        [u2 l2]= unit(-x2(1:3,:)); %unit vector and distance from tdrss to center of Earth
        if isfield(tdrss.sat(n),'antenna') %if there is antenna data for tdrss spacecraft n
            if isfield(tdrss.sat(n).antenna,'maxrange') %if the max range is defined
                index = l1>tdrss.sat(n).antenna.maxrange;
                ys(:,index) = NaN;
            end
            if isfield(tdrss.sat(n).antenna,'halfangle') %if the half angle is defined
                index = acos(dot(u1,u2))>tdrss.sat(n).antenna.halfangle;
                ys(:,index) = NaN;
            end
        end

        %% Apply the Earth occultation constraint
        LimbAngle = asin(RMin./l2); %rad. angles to the edge of the atmosphere
        LimbDist  = sqrt(l2.^2-RMin^2); %km distances from tdrss to the edge of the atmosphere
        index1 = l1>LimbDist; %meaning the user is further from tdrss than the edge of the atmosphere
        index2 = acos(dot(u1,u2))<LimbAngle; %meaning the user is within the cone defined by tdrss and the atmosphere
        index = index1 & index2; %meaning the user is behind the Earth
        ys(:,index) = NaN;


        %% Combine with results from previous satellites
        %     y = [y; ys + yg]; %combine the two measurements
        %     H = [H; Hs]; %the ground leg drops out of the partial
        %     derivatives
        indstart = 1 + numtypes*(n-1);
        indstop  = numtypes*n;
        y(indstart:indstop, tind)    = ys+ yg;
        Hf=zeros(numtypes,1,length(tind));
        if eflag
            for nt=1:numtypes
                for dN=1:length(tind)
                    Hf(nt,1,dN) = (Hg(nt,:,dN)-Hs(nt,:,dN))*dx{n}(:,tind(dN));
                end
            end
            H(indstart:indstop,7,tind) = Hf;
        end
        H(indstart:indstop, 1:6, tind) = Hs;%the ground leg drops out of the partial derivatives
    end
end

%% Set the measurement covariance output
if nargout > 2,
    %get sigma out of the options
    sigmaDefault    = ones(1,size(y,1))*1e-3;
    sigma           = getOdtbxOptions(options, 'rSigma', sigmaDefault );
    if( length(sigma)~=length(sigmaDefault) )
        disp('WARNING: In tdrssmeas, length(rSigma) does not match total number of measurements');
        disp('   Setting rSigma = default (1e-3) for all measurements');
        sigma = sigmaDefault;
    end
    sigma           = diag(sigma);
    R = repmat(sigma.^2,[1,1,size(y,2)]);
end

end

%% Validation Test

function failed = tdrssmeas_test()
%gsList = createGroundStationList();
disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
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

% Set propagator information
epoch  = datenum([2008  9 26  0  0   .000]); %UTC epoch
tspan  = 0:600:3600;
x0     = [6878;0.00;0.00;0.00;0.00;8.0]; %in km & km/sec
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
% default information for propagating Earth-centric orbits using JAT
% force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

% Measurement Options
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %UTC
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useLightTime', false);
measOptions = setOdtbxOptions(measOptions,'EarthAtmMaskRadius',6478.12); %km
measOptions = setOdtbxOptions(measOptions,'tdrss', tdrss);

% Set Expected Results
y_expected = [
    77143.8676829857  2835.46917821229               NaN               NaN  77941.4862100122 -6249.96951871453;
    78228.5123284356 -20780.2992395048               NaN               NaN  79843.3003057668 -25420.8497100346;
    81412.4520066554 -32498.1098036057   81376.442501345   34456.437489037  83240.3532719398 -31819.6509120546;
    85144.6557209073 -31141.3680452223  77380.2041620424  33700.8409138789  86685.8437713931 -27168.8751418441;
    NaN               NaN  74065.7701163095  22895.1417907037               NaN               NaN;
    NaN               NaN  72406.4859646846  5497.55474322249               NaN               NaN;
    NaN               NaN  72840.7066213739 -12693.5108798617               NaN               NaN]';

% Display Expected Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Expected TDRSS Measurements (East, West, and Spare) by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','East Range','East Doppler','West Range','West Doppler','Spare Range','Spare Doppler');
fprintf('%4s %13s %14s %14s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:7
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),y_expected(:,n))
end
disp(' ')
disp(' ')

% Run tdrssmeas
[y,H,R] = tdrssmeas(t,x,measOptions); %#ok<NASGU>
% Even though I don't use H and R here, I still need to generate them to
% test the portion of tdrssmeas that calculates them. If they aren't
% requested as outputs, then the "nargout > 2" check won't be true.

% Display Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Calculated TDRSS Measurements (East, West, and Spare) by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','East Range','East Doppler','West Range','West Doppler','Spare Range','Spare Doppler');
fprintf('%4s %13s %14s %14s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:7
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),y(:,n))
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
