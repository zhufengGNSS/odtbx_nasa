function [y,H,R,AntLB] = tdrssmeas_basic(t,x,options,qatt)
% TDRSSMEAS_BASIC  Makes TDRSS based measurements.
%
%   [y,H,R] = TDRSSMEAS_BASIC(t,x,options) creates TDRSS measurements
% based on the information in OPTIONS. The measurement types that can be
% returned are range, and range-rate. See the OD Toolbox functions
% named LOSRANGE, and LOSRANGERATE for details of each measurement type.
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   gsID                string                 optional, Used with createGroundStationList()
%   gsList                                     optional, Defaults to list from createGroundStationList()
%   gsECEF              3xNumberOfGS           optional, causes gsID and
%                                              gsList to be ignored
%   epoch                 datenum              UTC time associated with the
%                                              start of the simulation
%   gsElevationConstraint degs                 Spacecraft elevation must be
%                                              above this value to return
%                                              measurement
%   useRange            {true(default), false} Return range measurement
%   useRangeRate        {true(default), false} Return rangerate measurement
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only used for Tropo and
%                                              for Ionosphere and Charged
%                                              Particles
%   rSigma              {(1xM),ones(1,M)*1e-3) Measurement covariance
%   useIonosphere       {true, false(default), @function_handle} Includes Ionospheric delays
%   useTroposphere      {true, false(default), @function_handle} Includes Tropospheric delays
%   useChargedParticle  {true, false(default), @function_handle} Includes Charged Particle Delays
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
%   Regression tests exist in test_tdrssmeas_basic.m.
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
%   Jason Schmidt       11/08/2013      Converted into tdrssmeas_basic
%                                       version that does pure goemetric
%                                       measurements and only does range
%                                       and rangerate measurements

%% Get configuration values from options and define defaults if required
% gsID and gsList are optional, gsID and gsList go hand in hand, if one is defined, so must the other
gsID         = getOdtbxOptions(options, 'gsID', [] ); 
gsList       = getOdtbxOptions(options, 'gsList', []);
% gsECEF is optional, If gsECEF is defined, gsID and gsList will be ignored
gsECEF       = getOdtbxOptions(options, 'gsECEF', []); 
% epoch in datenum format, required
epoch = getOdtbxOptions(options, 'epoch', NaN); 
if isnan(epoch); error('An epoch must be set in the options structure.'); end
% tdrss struct is optional, default KOEs will be used if not defined
tdrss = getOdtbxOptions(options, 'tdrss', []); 
% RMin Default is the Mean Radius of the Earth plus 100 km of Atmosphere
RMin = getOdtbxOptions(options, 'EarthAtmMaskRadius', 100+JATConstant('meanRadius','Earth')/1000); 
% Output both range and range rate by default
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
numtypes     = useRange + useRangeRate;
% Frequency of the signal
f = getOdtbxOptions(options,'frequencyTransmit',[]);
if isempty(f)
    f = 1475000000;%default to Single Access 1 (SA1) K-band frequency
    options = setOdtbxOptions(options,'frequencyTransmit',f);
end

% make sure t is a row vector
if size(t,1)==length(t)
    t=t'; 
end
N=length(t);

% Populate tdrss with default Keplerian Orbital Elements for TDRS if empty
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

% sigma (standard deviation of measurements) will default to 1e-3
sigmaDefault    = repmat(1e-3,1,length(tdrss.sat)*numtypes);
sigma           = getOdtbxOptions(options, 'rSigma', sigmaDefault );

%Initialize Y and H for speed
y    = nan(length(tdrss.sat)*numtypes,length(t));
H    = zeros(length(tdrss.sat)*numtypes,6,length(t));
%Hardcode ground station elevation constraint to 0 degrees
options  = setOdtbxOptions(options,'gsElevationConstraint', 0);

% populate link budget options
if isfield(options,'linkbudget') && ~isempty(options.linkbudget)
    dolinkbudget = true;
    link_budget = getOdtbxOptions(options, 'linkbudget',[]);
    
    d2r = pi/180;
    % Set some default values
    link_budget = linkbudget_default(link_budget, 'AntennaPattern', {'omni.txt'});
        %  Specify antenna pattern for each antenna, existing antennas are:
        %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
        %     omni.txt                  - zero dB gain,  180 degree half beamwidth
        %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
        %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
        %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
    num_ant = length(link_budget.AntennaPattern); %hasn't been tested for >4 antennas
    link_budget = linkbudget_default(link_budget, 'RXAntennaMask', 180*d2r);
    link_budget = linkbudget_default(link_budget, 'AtmosphereMask', 0); % km
        %  Troposphere mask radius ~50 km
        %  Ionosphere mask radius ~(500-1000 km)
    link_budget = linkbudget_default(link_budget, 'NoiseTemp', 300); % K
        % System noise temp [K], space pointing antenna = 290
        % System noise temp [K], earth pointing antenna = 300
    link_budget = linkbudget_default(link_budget, 'AtmAttenuation', 0.0); % dB
        % attenuation due to atmosphere (should be negative) [dB]
    link_budget = linkbudget_default(link_budget, 'TransPowerLevel', 2); % 1-minimum, 2-typical, 3-max
    link_budget = linkbudget_default(link_budget, 'TransPowerOffset', 0.0); % dB, global offset
    link_budget = linkbudget_default(link_budget, 'TXAntennaMask', 70*d2r );  % in rad
        %  The actual mask used is the lesser of this mask and the limit of the defined pattern
        %  Note:  mask = 70 deg includes entire defined pattern
        %         mask = 42 deg includes only main and first side lobes
        %         mask = 26 deg includes only main lobe
    link_budget = linkbudget_default(link_budget, 'ReceiverNoise', -3 );  % dB, Noise figure of receiver/LNA
    link_budget = linkbudget_default(link_budget, 'RecConversionLoss', -1.5 );  % dB
        % Receiver implementation, A/D conversion losses [dB]
        %   Novatel: L = -4.0 	
        %   Plessey: L = -1.5		
    link_budget = linkbudget_default(link_budget, 'SystemLoss', 0 ); % dB, System losses, in front of LNA
    link_budget = linkbudget_default(link_budget, 'LNAGain', 40 ); % dB, LNA gain (trimble pre-amp spec = 42-48)
    link_budget = linkbudget_default(link_budget, 'CableLoss', -2 ); % dB, Cable losses (after LNA)
    link_budget = linkbudget_default(link_budget, 'RecAcqThresh', 32 ); % dB-Hz, Receiver acquisition threshold
    link_budget = linkbudget_default(link_budget, 'RecTrackThresh', link_budget.RecAcqThresh ); % dB-Hz, Receiver tracking threshold
    link_budget = linkbudget_default(link_budget, 'DynamicTrackRange', 15 ); % dB
        % Dynamic tracking range of receiver, or maximum difference  
        % in power levels tracked simultaneously. If the difference
        % in snrs between two satellites is more that link_budget.DynamicTrackRange,
        % the weaker of the two will not be considered visible. 
    % Reassign the options structure with any changed/default link budget values
    options = setOdtbxOptions(options, 'linkbudget', link_budget);
    
    
    if nargin < 4
        wstr = ['S/C body quaternion is not specified in the input.  ',...
            'Body will be assumed to be aligned with the same coordinate ',...
            'frame as the position and velocity states.'];
        warning('ODTBX:GSMEAS:noBodyQuat',wstr);
        qatt = repmat([0;0;0;1],1,N);
    end
else
    dolinkbudget = false;
end

%% Get TDRS states
if strcmpi(tdrss.type, 'keplerian') % Simple, but less accurate ephemeris
    tdrss.t = t; %seconds from scenario epoch

    for n=1:length(tdrss.sat);
        %Set the time span for each tdrs vehicle
        tdrss.sat(n).tspan = tdrss.t;
        epdiff = (epoch-tdrss.sat(n).epoch)*86400; %secs from tdrss epoch to scenario epoch
        %Define position and velocity of each tdrs vehicle
        GM         = JATConstant('muEarth')/10^9;
        tdrss.x{n} = kep2cart(kepprop2b(tdrss.sat(n),tdrss.sat(n).tspan+epdiff,GM),GM);
    end

elseif strcmpi(tdrss.type, 'spephem') % More accurate ephemeris than keplarian
    for n=1:length(tdrss.sat) 
        % populate time and position for each tdrs vehicle
        [tdrss.t, tdrss.x{n}] = read_spephem(tdrss.sat(n).filename);
        tdrss.sat(n).epoch   = tdrss.t(1); %UTC
        epdiff               = (epoch-tdrss.sat(n).epoch)*86400; %secs from tdrss epoch to scenario epoch
        tdrss.sat(n).tspan   = (tdrss.t - tdrss.t(1))*86400 - epdiff;
    end
else
    error('tdrssmeas_basic:invalid_tdrss_type',[tdrss.type ' is not a supported state type'])
end

%% Get Ground Station States
if isempty(gsECEF)
    if isempty(gsID)
        gsID = {'WSGT','GTSS'}; %default ground stations
        gsECEF = [-1539.37859986625         -5070.00944053254
                  -5160.93110611279          3569.07519673414
                   3408.22918525621          1491.68888437059]; % hardcoded for speed
    end
    if isempty(gsECEF)
        if isempty(gsList)
            gsList = createGroundStationList(); 
        end
        gsECEF = zeros(3,length(gsID));
        for n=1:length(gsID)
                gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
        end
    end
end

%% Loop over all TDRS Satellites 
% to obtain the measurements (y) and measurement sensitivity matrix (H)

for n=1:length(tdrss.x)
    %% Set the ground station used by this TDRS vehicle
    % Because TDRS vehicles are in GEO, the ground station will not change
    Ephem.satPos      = tdrss.x{n}(1:3,1)*1000; %ECI satellite coordinates (m)
    Ephem.SatCoords   = 'ECI';
    Ephem.Epoch       = epoch+tdrss.sat(n).tspan(1)/86400;
    Ephem.StationInfo = 'ECEF';
    Ephem.staPos      = gsECEF*1000;
    [~,Elv]  = jatStaAzEl(Ephem);
    % Use the ground station that has the highest elevation angle to the
    % TDRS vehicle
    gsECEFn  = gsECEF(:, max(Elv) == Elv );
    options  = setOdtbxOptions(options,'gsECEF',gsECEFn);

    %% set the tracking schedule for this TDRS vehicle
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

    if ~isempty (t1)
        %% Get the measurements for the Ground and Space leg
        if isequal(tdrss.sat(n).tspan,t1)
            x2 = tdrss.x{n};
        else
            x2 = interp1(tdrss.sat(n).tspan,tdrss.x{n}',t1,'spline')'; %interpolate the TDRS position to time t
        end
        [ys, Hs] = rrdotang(t1,x1,x2,options); %calculate the Space leg

        
        gsmeas_options = options;
        gsmeas_options.linkbudget = []; % turn off link budget for ground leg
        yg = gsmeas(t1,x2,gsmeas_options); %calculate the Ground leg


        %% Apply the antenna constraints
        [u1, l1] = unit(x1(1:3,:)-x2(1:3,:)); %unit vector and distance from tdrss to user
        [u2, l2]= unit(-x2(1:3,:)); %unit vector and distance from tdrss to center of Earth
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

        %% Perform link budget analysis
        if dolinkbudget
            % get Link Budget
            [AntLB_new, HVIS_new] = get_linkbudget(t, x1, x2, options, qatt);

            if n~=1
                AntLB = Append_AntLB(AntLB,AntLB_new);
                HVIS(n,:) = HVIS_new;
            else
                AntLB = AntLB_new;
                HVIS = [HVIS_new;
                        NaN(length(tdrss.x)-1,size(HVIS_new,2))];
            end
            % Make measurements that are not visible NaNs
            kk=1;
            for ii=1:size(HVIS_new,1)
                for jj=1:size(HVIS_new,2)
                    if HVIS_new(ii,jj)==0
                        ys(kk:kk+numtypes-1,jj) = NaN;
                    end
                end
                kk=kk+numtypes;
            end
        end

        %% Combine with results from previous satellites
        %     y = [y; ys + yg]; %combine the two measurements
        %     H = [H; Hs]; %the ground leg drops out of the partial
        %     derivatives
        indstart = 1 + numtypes*(n-1);
        indstop  = numtypes*n;
        y(indstart:indstop, tind)    = ys+ yg;
        H(indstart:indstop, :, tind) = Hs;%the ground leg drops out of the partial derivatives
    end
end

%% Set the measurement covariance (R)
if nargout > 2,
    if( length(sigma)~=length(sigmaDefault) )
        warning('tdrssmeas_basic:invalid_rSigma',' In tdrssmeas, length(rSigma) does not match total number of measurements. Setting rSigma = default (1e-3) for all measurements');
        sigma = sigmaDefault;
    end
    sigma           = diag(sigma);
    R = repmat(sigma.^2,[1,1,size(y,2)]);
end

end

function  AntLB = Append_AntLB(AntLB,AntLB_new)
    for ii=1:numel(AntLB)
        AntLB{ii}.Halpha_r  = [AntLB{ii}.Halpha_r;
                           AntLB_new{ii}.Halpha_r];
        AntLB{ii}.Hvis_beta = [AntLB{ii}.Hvis_beta;
                           AntLB_new{ii}.Hvis_beta]; 
        AntLB{ii}.Hvis_CN0  = [AntLB{ii}.Hvis_CN0;
                           AntLB_new{ii}.Hvis_CN0]; 
        AntLB{ii}.HCN0      = [AntLB{ii}.HCN0;
                           AntLB_new{ii}.HCN0]; 
        AntLB{ii}.HAd       = [AntLB{ii}.HAd;
                           AntLB_new{ii}.HAd]; 
        AntLB{ii}.HAr       = [AntLB{ii}.HAr;
                           AntLB_new{ii}.HAr]; 
        AntLB{ii}.HAP       = [AntLB{ii}.HAP;
                           AntLB_new{ii}.HAP];
        AntLB{ii}.HRP       = [AntLB{ii}.HRP;
                           AntLB_new{ii}.HRP];
        AntLB{ii}.HAt       = [AntLB{ii}.HAt;
                           AntLB_new{ii}.HAt];

    end
end