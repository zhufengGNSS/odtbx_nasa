function [y,H,R] = lnrmeas(t,x,options)
% LNRMEAS  Makes Lunar Relay based measurements.
%
%   [y,H,R] = LNRMEAS(t,x,options) creates Lunar Relay measurements
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
%   epoch                 datenum              UTC time associated with start
%                                              of simulation
%   useRange            {true(default), false} Use range measurement
%   rangeType           {'1way','1wayFWD','1wayRTN','2way'(default)}
%                                              - 1way=1wayFWD
%   useRangeRate        {true(default), false} Use range rate measurement
%   useDoppler          {true, false(default)} Use doppler measurement
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only used for Tropo and
%                                              for Doppler
%   rSigma              {(1xM),ones(1,M)*1e-3) Measurement covariance
%   useLightTime        {true, false(default)} Include light time delays
%   useIonosphere       {true, false(default)} Includes Ionospheric delays
%   useTroposphere      {true, false(default)} Includes Tropospheric delays
%   relay               structure(km, radians) Structure of all Lunar Relay
%                                              data
%       STK Ephem Version:
%           relay.type            = 'STKEphem';
%           relay.sat(1).filename = 'LR1.e'; %Lunar Relay Satellite #1
%           relay.sat(2).filename = 'LR2.e'; %Lunar Relay Satellite #2
%           ...
%       Optional Antenna Data:
%           relay.sat(i).antenna.maxrange = 70000; %max dist to user in km
%           relay.sat(i).antenna.halfangle = 15*pi/180; %Nadir facing field
%               of view in radians
%           relay.MoonMaskRadius = 20+JATConstant('meanRadius','Moon')/1000
%               The radius of Moon Obscuration. Default is the Mean Radius
%               of the Moon plus 20 km to account for lunar mountains
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
%   TDRSSMEAS, READ_STKEPHEM
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
%   Kevin Berry         12/29/2008      Original
%   Kevin Berry         05/28/2009      Added Moon occultation
%                                       Added antenna constraints
%   Kevin Berry         06/16/2009      Added self test
%   Kevin Berry         06/25/2009      Added time scale comments
%   Allen Brown         09/22/2009      Updated regression/validation test
%                                       data to match leap second data for
%                                       1 Jan 2009.
%   Ravi Mathur         08/26/2012      Updated regression/validation test
%                                       data to match leap second data for
%                                       1 Jul 2012.

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    y = lnrmeas_test(t);  
else
	[y,H,R] = Get_lnrmeas(t,x,options);
end
end


%% Main function
function [y,H,R] = Get_lnrmeas(t,x,options)

%% Get values from options
epoch = getOdtbxOptions(options, 'epoch', NaN);
uselt = getOdtbxOptions(options, 'useLightTime', false);
relay = getOdtbxOptions(options, 'relay', []);
if isnan(epoch); error('An epoch must be set in the options structure.'); end
if isempty(relay); error('relay state information must be set in the options structure'); end
if size(t,1)==length(t), t=t'; end

%% Get Lunar Radius used for Occultation
if isfield(relay,'MoonMaskRadius') %if the Moon Mask Radius is defined
    RMin = relay.MoonMaskRadius; %Assumed to be kilometers
else
    RMin = 20+JATConstant('meanRadius','Moon')/1000;
        %Default is the Mean Radius of the Moon plus 20 km to account for
        %lunar mountains
end

%% Get relay states
if strcmpi(relay.type, 'stkephem')
    for n=1:length(relay.sat)
        try
            [relay.t relay.x{n} hdrInfo] = read_stkephem(relay.sat(n).filename);
        catch me
            disp(sprintf('lnrmeas:read_stkephem failure.  Failed to read file %s\n.  Original error follows:\n', relay.sat(n).filename));
            rethrow(me);
        end
        if isfield(hdrInfo, 'CentralBody')
            if ~strcmpi(hdrInfo.CentralBody,'Earth')
                error([hdrInfo.CentralBody ' is not a supported central body '...
                    'for this measurement model. This model requires an Earth '...
                    'based ephemeris for each object.'])
            end
        end
        if isfield(hdrInfo, 'CoordinateSystem')
            if ~strcmpi(hdrInfo.CoordinateSystem,'J2000')
                error([hdrInfo.CoordinateSystem ' is not a supported coordinate '...
                    'system for this measurement model. This model requires a '...
                    'J2000 coordinate system for each object.'])
            end
        else
            error(['"Fixed" is not a supported coordinate system'...
                'for this measurement model. This model requires a '...
                'J2000 coordinate system for each object.'])
        end
        if ~strcmpi(hdrInfo.Format,'EphemerisTimePosVel')
            error([hdrInfo.Format ' is not a supported format for this '...
                'measurement model. This model requires EphemerisTimePosVel '...
                'for the stk ephemeris file format.'])
        end
        if isfield(hdrInfo, 'DistanceUnit')
            if strcmpi(hdrInfo.DistanceUnit,'Meters')
                relay.x{n} = relay.x{n}./1000;
            elseif strcmpi(hdrInfo.DistanceUnit,'Kilometers')
                % Do nothing... this function runs in kilometers
            else
                error([hdrInfo.DistanceUnit ' is not a supported unit '...
                    'for this measurement model. Feel free to add the '...
                    'unit to the code yourself and convert it to kilometers.'])
            end
        else
            %Default DistanceUnit in STK Ephems = Meters
            relay.x{n} = relay.x{n}./1000;
        end
        relay.sat(n).epoch = hdrInfo.Epoch;
        epdiff               = (epoch-relay.sat(n).epoch)*86400; %secs from relay epoch to scenario epoch
        relay.sat(n).tspan   = relay.t - epdiff;
    end
else
    error([relay.type ' is not a supported state type'])
end

%% Get the position of the Moon in J2000
x_Moon = zeros(3,length(t));
for m=1:length(t)
    x_Moon(1:3,m) = ephemDE405('Geocentric_Moon',epoch+t(m)/86400,'UTC');
end

%% Loop over all relay satellites
y = [];
H = [];
R = [];
for n=1:length(relay.x)
    %% Get the measurements
    if uselt
        [ys, Hs, Rs] = rrdotlt(t, x, relay.sat(n).tspan, relay.x{n},...
            options); %calculates the Space leg
    else
        x2 = interp1(relay.sat(n).tspan,relay.x{n}',t,'spline')'; %interpolate the relay's position to time t
        [ys, Hs, Rs] = rrdotang(t,x,x2,options); %calculates the Space leg
    end
    
    %% Apply the antenna constraints
    [u1 l1] = unit(x(1:3,:)-x2(1:3,:)); %unit vector and distance from relay to user
    [u2 l2] = unit(x_Moon(1:3,:)-x2(1:3,:)); %unit vector and distance from relay to center of Moon
    if isfield(relay.sat(n),'antenna') %if there is antenna data for relay spacecraft n
        if isfield(relay.sat(n).antenna,'maxrange') %if the max range is defined
            index = l1>relay.sat(n).antenna.maxrange;
            ys(:,index) = NaN;
        end
        if isfield(relay.sat(n).antenna,'halfangle') %if the half angle is defined
            index = acos(dot(u1,u2))>relay.sat(n).antenna.halfangle;
            ys(:,index) = NaN;
        end
    end
    
    %% Apply the Moon occultation constraint
    LimbAngle = asin(RMin./l2); %rad. angles to the edge of the mask
    LimbDist  = sqrt(l2.^2-RMin^2); %km distances from relay to the edge of the mask
    index1 = l1>LimbDist; %meaning the user is further from the relay than the edge of the mask
    index2 = acos(dot(u1,u2))<LimbAngle; %meaning the user is within the cone defined by the relay and the mask
    index = index1 & index2; %meaning the user is behind the Moon
    ys(:,index) = NaN;

    %% Combine with results from previous relays
    y = [y; ys];
    H = [H; Hs];
    R = [R; Rs];

end

end

%% Validation/Regression Test

function failed = lnrmeas_test(in_string)
disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
% Input Relay States as STKEphem type ephemeris files
relay.type            = 'STKEphem';
if strcmpi(in_string,'ValidationTest')
    relay.sat(1).filename = 'Relay1.e'; %Lunar Relay Satellite #1
    relay.sat(2).filename = 'Relay2.e'; %Lunar Relay Satellite #2
    if ~exist(relay.sat(1).filename)||~exist(relay.sat(2).filename)
        error('Validation test files could not be found.');
    end
elseif strcmpi(in_string,'RegressionTest')
    relay.sat(1).filename = 'Relay1_regress.e'; %Lunar Relay Satellite #1
    relay.sat(2).filename = 'Relay2_regress.e'; %Lunar Relay Satellite #2
    if ~exist(relay.sat(1).filename)||~exist(relay.sat(2).filename)
        msg = ['Regression test files could not be found. Do you have the ',...
            'regression test folder in your path? Running the validation ',...
            'test instead...'];
        warning('LNRMEAS:Regression1',msg)
        lnrmeas('ValidationTest');
        msg = ['Regression test could not be ran because the regression ',...
            'ephems could not be found. The validation case was ran ',...
            'instead. The regression test is being considered a failure ',...
            'regardless of validation case results!'];
        warning('LNRMEAS:Regression2',msg)
        failed = 1;
        return
    end
end
    
% Set the radius of Moon Obscuration
relay.MoonMaskRadius = 20+JATConstant('meanRadius','Moon')/1000;
%Using the Mean Radius of the Moon plus 20 km to account for lunar mountains

% Input relay antenna data
relay.sat(1).antenna.maxrange  = 100000; %km
relay.sat(2).antenna.maxrange  = 100000; %km

% Set propagator information
epoch  = datenum('1 Dec 2016 00:00:00.000'); %UTC
tspan  = 0:600:7200;
x0     = [-3.10570021720167e+004
    -3.80943874664242e+005
    -1.28680660459686e+005
    1.05814791201683e+000
    1.44219161241440e+000
    6.24467648736737e-001]; %in km & km/sec

jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %UTC in datenum format
jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM2');
jOptions = setOdtbxOptions(jOptions, 'gravDeg', 2, 'gravOrder', 0); %max 20x20
jOptions = setOdtbxOptions(jOptions, 'useSolarGravity', true);
jOptions = setOdtbxOptions(jOptions, 'useLunarGravity', true);
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
measOptions = setOdtbxOptions(measOptions,'relay', relay);

% Load Expected Results
% NOTE: These results are based on the most recent July 2012 leap second. If a
% more recent leap second is added, then these results must be updated.
load lnrmeas_expected.mat;

% Display Expected Results
fprintf('%s\n',char(ones(1,66)*'-'));
disp('Expected Relay Measurements by time:')
fprintf('%s\n',char(ones(1,66)*'-'));
fprintf('%-6s %14s %14s %14s %14s\n','Time','#1 Range','#1 Doppler','#2 Range','#2 Doppler');
fprintf('%4s %13s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,66)*'-'));
for n=1:13
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f\n',t(n),y_expected(:,n))
end         
disp(' ')
disp(' ')

% Run lnrmeas
[y,H,R] = lnrmeas(t,x,measOptions);

% Display Results
fprintf('%s\n',char(ones(1,66)*'-'));
disp('Calculated Relay Measurements by time:')
fprintf('%s\n',char(ones(1,66)*'-'));
fprintf('%-6s %14s %14s %14s %14s\n','Time','#1 Range','#1 Doppler','#2 Range','#2 Doppler');
fprintf('%4s %13s %14s %14s %14s\n','(s)','(km)','(Hz)','(km)','(Hz)');
fprintf('%s\n',char(ones(1,66)*'-'));
for n=1:13
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f\n',t(n),y(:,n))
    %fprintf('%-6i %20.11f %20.11f %20.11f %20.11f\n',t(n),y(:,n))
end     

passed = tol > max( max( abs( y - y_expected ) ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed! Please check that the most recent leap second is July 2012.')
else
    disp(' ')
    disp('Test Passed.')
end

end
