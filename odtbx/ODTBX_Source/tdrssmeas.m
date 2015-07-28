function [y,H,R] = tdrssmeas(t,x,options)
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
%           tdrss.sat(i).epoch = datenum('30 Apr 2012 18:00:00.000'); %UTC
%           tdrss.sat(i).sma   = 42164.1363; %semi-major axis in km
%           tdrss.sat(i).ecc   = 0; %eccentricity
%           tdrss.sat(i).incl  = 0; %inclination in radians
%           tdrss.sat(i).raan  = 88*pi/180; %right ascension of the
%                                            ascending node in radians
%           tdrss.sat(i).argp  = 0; %argument of periapse in radians
%           tdrss.sat(i).tran  = 0; %true anomaly in radians
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
% VALIDATION/REGRESSION TEST
%
%   These tests have been moved to EarthOrbitPlot_test.m to conform to
%   the new regression testing format.
%
%   keyword: measurement
%   See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, ODTBXOPTIONS, GSMEAS,
%   READ_SPEPHEM
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
%   Kevin Berry         04/30/2012      Added default TDRSS locations
%   Ravi Mathur         08/28/2012      Extracted regression test
%   Russell Carpenter   07/15/2015      Changed default ground stations
%                                       and fixed some bugs

%% Get values from options
gsID         = getOdtbxOptions(options, 'gsID', [] );
gsList       = getOdtbxOptions(options, 'gsList', []);
gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
epoch = getOdtbxOptions(options, 'epoch', NaN); %UTC
uselt = getOdtbxOptions(options, 'useLightTime', false);
tdrss = getOdtbxOptions(options, 'tdrss', []);
RMin = getOdtbxOptions(options, 'EarthAtmMaskRadius', 100+JATConstant('meanRadius','Earth')/1000);
%Default is the Mean Radius of the Earth plus 100 km of Atmosphere
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
numtypes     = useRange + useRangeRate + useDoppler;

if isnan(epoch); error('An epoch must be set in the options structure.'); end
if size(t,1)==length(t), t=t'; end
if isempty(tdrss) % Use default tdrss locations
    % From http://nssdc.gsfc.nasa.gov/multi/tdrs.html, The operational 
    % spacecraft are located at 41°, 174° and 275° West longitude
    w = [0;0;JATConstant('wEarth')];
    D = jatDCM('ecef2eci', epoch, 0);
    M = rotransf(-D*w,D);
    
    tdrss.type = 'Keplerian';
    tdrss.sat(1).epoch = epoch; %UTC
    tdrss.sat(1).sma   = 42164.1363; %semi-major axis in km
    tdrss.sat(1).ecc   = 0; %eccentricity
    tdrss.sat(1).incl  = 0; %inclination in radians
    X_ecef=[42164.1363*cosd(-41);42164.1363*sind(-41);0;0;0;0];
    X_eci = M(1:6,1:6)*X_ecef;
    tdrss.sat(1).raan  = atan2(X_eci(2),X_eci(1)); %right ascension of the ascending node in radians
    tdrss.sat(1).argp  = 0; %argument of periapse in radians
    tdrss.sat(1).tran  = 0; %true anomaly in radians

    tdrss.sat(2) = tdrss.sat(1);    
    X_ecef=[42164.1363*cosd(-174);42164.1363*sind(-174);0;0;0;0];
    X_eci = M(1:6,1:6)*X_ecef;
    tdrss.sat(2).raan  = atan2(X_eci(2),X_eci(1)); %right ascension of the ascending node in radians
    
    tdrss.sat(3) = tdrss.sat(1);  
    X_ecef=[42164.1363*cosd(-275);42164.1363*sind(-275);0;0;0;0];
    X_eci = M(1:6,1:6)*X_ecef;
    tdrss.sat(3).raan  = atan2(X_eci(2),X_eci(1)); %right ascension of the ascending node in radians
end

%% Get TDRS state
if strcmpi(tdrss.type, 'keplerian')
    % A TDRS ephem is created at user times
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
        tdrss.x{n} = kep2cart(kepprop2b(tdrss.sat(n),tdrss.sat(n).tspan+epdiff,GM),GM);
    end

elseif strcmpi(tdrss.type, 'spephem')
    % TDRS ephem times not necessarily same as user times!
    for n=1:length(tdrss.sat) % More accurate
        [tdrss.t,tdrss.x{n}] = read_spephem(tdrss.sat(n).filename);
        tdrss.sat(n).epoch   = tdrss.t(1); %UTC
        epdiff               = (epoch-tdrss.sat(n).epoch)*86400; %secs from tdrss epoch to scenario epoch
        tdrss.sat(n).tspan   = (tdrss.t - tdrss.t(1))*86400 - epdiff;
    end
else
    error([tdrss.type ' is not a supported state type'])
end

%% Get Ground Station States
if isempty(gsECEF)
    if isempty(gsList)
        gsList = createGroundStationList(); 
    end
    if isempty(gsID)
        gsID = {'WHSK','GWMK'};% Changed from {'WSGT','GTSS'} per Cheryl Gramling;
    end
    gsECEF = zeros(3,length(gsID));
    for n=1:length(gsID)
        gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
    end

end

%% Loop over all TDRS Satellites
%initialize Y and H and R
y    = nan(length(tdrss.sat)*numtypes,length(t));
H    = zeros(length(tdrss.sat)*numtypes,6,length(t));
R = [];

options  = setOdtbxOptions(options,'gsElevationConstraint', 0);
for n=1:length(tdrss.x)
    %% Set the ground station to the one with max elevation
    Ephem.satPos      = tdrss.x{n}(1:3,1)*1000; %ECI satellite coordinates (m)
    Ephem.SatCoords   = 'ECI';
    Ephem.Epoch       = epoch+tdrss.sat(n).tspan(1)/86400;
    Ephem.StationInfo = 'ECEF';
    Ephem.staPos      = gsECEF*1000;
    [~,Elv]  = jatStaAzEl(Ephem);
    gsECEFn  = gsECEF(:, max(Elv) == Elv );
    options  = setOdtbxOptions(options,'gsECEF',gsECEFn);
    %% Put Tdrss states into an array
    %for j=size(tdrss.x,2):-1:1
    %    xt(:,:,j)=tdrss.x{1,j}(:,:);
    %end
    xt=tdrss.x{1,n}(:,:);

    %% set the tracking schedule
    % Check schedule for times when tracking is done
    if isfield(tdrss.sat(n),'schedule')
        tSch = NaN(size(tdrss.sat(n).schedule,1),2);%Peallocated for speed
        %this is where the schedule is put into a format that is easier to
        %work with
        for j=1:size(tdrss.sat(n).schedule,1)
            tSch (j,:) = ((datenum(tdrss.sat(n).schedule(j,:))-epoch)*86400);
        end
        tind1 = [];
        for m=1:size(tSch,1)
            tind1 = union(tind, find( tSch(m,1)<=t & t<=tSch(m,2) ));
        end
    else
        tind1 = 1:length(t);
    end
    % assume TDRS ephem starts before and ends after user ephem
    iafter1 = find(tdrss.sat(n).tspan<=t(1),1,'last');
    ibefend = find(tdrss.sat(n).tspan>=t(end),1,'first');
    tind2 = iafter1:ibefend;
    % t, x, and xt are restricted only to the times when tracking is done
    t1   = t(tind1);
    x1   = x(:,tind1);
    t2   = tdrss.sat(n).tspan(tind2);
    x2  = xt(:,tind2);%,n);
    if ~isempty (t1)
        %% Get the measurements for the Ground and Space leg
        if uselt
            %[ys, Hs, ~, t2_lt, x2_lt] = rrdotlt(t1, x1, tdrss.sat(n).tspan, tdrss.x{n},...
            %    options); %calculate the Space leg
            [ys, Hs, ~, t2_lt, x2_lt] = rrdotlt(t1, x1, t2, x2, options); %calculate the Space leg
            if length(t2_lt)==2 %then it was a 2way measurement
                options       = setOdtbxOptions(options,'rangeType','1wayFWD');
                yg1     = gsmeas(t2_lt{1}, x2_lt{1}, options); %calculate the FWD Ground leg
                options       = setOdtbxOptions(options,'rangeType','1wayRTN');
                yg2           = gsmeas(t2_lt{2}, x2_lt{2}, options); %calculate the RTN Ground leg
                yg            = (yg1+yg2)/2; %average the ground legs
            else %else it was a 1way measurement and gsmeas will know which way
                yg = gsmeas(t2_lt{1}, x2_lt{1}, options); %calculate the Ground legl
            end
        else
            %if isequal(tdrss.sat(n).tspan,t1)
            %    x2 = tdrss.x{n};
            %else
            %    x2 = interp1(tdrss.sat(n).tspan,tdrss.x{n}',t1,'spline')'; %interpolate the TDRS position to time t1
            %end
            if ~isequal(t1,t2)
                x2 = interp1(t2,x2',t1,'spline')';
            end  
            [ys, Hs] = rrdotang(t1,x1,x2,options); %calculate the Space leg
            yg = gsmeas(t1,x2,options); %calculate the Ground leg
        end

        %% Apply the antenna constraints
        if uselt
            % Just apply for uplink (FWD ground leg) since if not satsified
            % on uplink there can be no return, and cases for which it
            % works on uplink and not downlink are less likely.
            x2 = x2_lt{1}; 
        end
        [u1,l1] = unit(x1(1:3,:)-x2(1:3,:)); %unit vector and distance from tdrss to user
        [u2,l2]= unit(-x2(1:3,:)); %unit vector and distance from tdrss to center of Earth
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
        y(indstart:indstop, tind1)    = ys+ yg;
        H(indstart:indstop, :, tind1) = Hs;%the ground leg drops out of the partial derivatives
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