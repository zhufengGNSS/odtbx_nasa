function rangeRate = LOSRangeRate(t, r1, v1, r2, v2, options)

% LOSRANGERATE  Line of sight range rate measurement between two objects 
%
%   rangeRate = LOSRangeRate(r1, v1, r2, v2) returns the 
% range rate measurement between two vehicles and is based on the 
% GEONS DATSIM Math Spec.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The parameters that are
% valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   useGPSIono          {true, false(default)} only for GPS sats as x2
%   useIono             {true, false(default), @function_handle} true only for groundstats as x2
%   useTropo            {true, false(default), @function_handle} true only for groundstats as x2
%   useChargedParticle  {true, false(default), @function_handle} true only for groundstats as x2
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only for Doppler and
%                                              measurement errors
%   epoch                datenum               UTC time associated with   
%                                              start of simulation.
%   INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t         (1xN)	Times corresponding to r1 (secs)
%      r1        (3xN)	User spacecraft position (km)
%      v1        (3xN)  User spacecraft velocity (km/s)
%      r2        (3xN)  Tracking spacecraft/ground station position (km)
%      v2        (3xN)  Tracking spacecraft/ground station velocity (km/s)
%      options   (1x1)  Data structure
%
% OUTPUTS
%      rangeRage (1xN)  range rate (km/s)
%
% VALIDATION/REGRESSION TEST
%
%  These have been moved to LOSRangeRate_test.m in the regression testing
%  framework to conform with the new testing format.
%
% keyword: measurement
% See also LOSRANGE, LOSDOPPLER, RRDOT, RRDOTLT
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
%   Derek Surka         08/24/2007      Original
%   Kevin Berry         05/08/2008      generalized for 3xN inputs
%   Kevin Berry         05/22/2008      Added JAT Iono and Tropo models
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         09/08/2008      Added Charged Particle Model
%   Kevin Berry         12/23/2008      Replaced elevation.m with
%                                       jatStaAzEl.m
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy in 
%                                       GPSIono
%   Ravi Mathur         08/28/2012      Extracted regression test
%   Jason Schmidt       11/7/2013       Encapsulated basic Tropo, Iono and
%                                       ChPart models in their own
%                                       functions.  Added support for user
%                                       defined models.

useGPSIono = getOdtbxOptions(options, 'useGPSIonosphere', false );
useIono    = getOdtbxOptions(options, 'useIonosphere', false );
useTropo   = getOdtbxOptions(options, 'useTroposphere', false );
useChPart  = getOdtbxOptions(options, 'useChargedParticle', false );

%% Define Models
% If use* is set to true, use basic model, if false, leave model undefined,
% if a function handle, use that function handle
if isa(useIono,'function_handle')
    IonoModel = useIono;
    useIono = true;
elseif useIono==true
    IonoModel = @IonoModel_basic; %Ionospheric range rate errors (only valid for ground-based tracking) 
end
    
if isa(useTropo,'function_handle')
    TropoModel = useTropo;
    useTropo = true;
elseif useTropo==true
    TropoModel = @TropoModel_basic; % Tropospheric range rate errors (only valid for ground-based tracking)
end

if isa(useChPart,'function_handle')
    ChPartModel = useChPart;
    useChPart = true;
elseif useChPart==true
    ChPartModel = @ChPartModel_basic; % Charged Particle range rate errors (only used for ground-based
                                      % tracking where the satellite is outside the magnetosphere)
end

%% Calculate Geometric rangerate and other constants
c           = JATConstant('c')/1000; %Speed of light in km
r           = r1 - r2 ;
losVector   = unit( r );
rangeRate   = dot( losVector, (v1-v2) ) ./ (1 - dot(losVector,v2)/c );
drdotTropo  = 0;
drdotIono   = 0;
drdotChPart = 0;

%% Calculate Measurement Errors Due to Troposphere
if useTropo 
    %Approximated based on a finite difference of the range error.      
    drTropo = TropoModel(t,r1,r2,options); 
    drTropo2 = TropoModel(t+1,r1+v1,r2+v2,options);
    drdotTropo = drTropo2 - drTropo; % change in Tropospheric effect
    rangeRate  = rangeRate + drdotTropo;
end

%% Calculate Measurement Errors Due to Ionosphere
if useGPSIono % GPS Ionospheric range rate errors (only valid for GPS)
    %Approximated based on a finite difference of the range error.
    %r1 and r2 states must be in an ECI.
    epoch     = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    epoch_GPS = convertTime('GPS','UTC',epoch); %GPS Time
    
    drIono = zeros(1,length(range)); %initialize before loop
    for n = 1:length(t)
        drIono(n) = jatGPSIonoModel(r1(:,n)*1000, r2(:,n)*1000, epoch_GPS+t(n)/86400)/1000;
    end

    drIono2 = zeros(1,length(range)); %initialize before loop
    for n = 1:length(t)
        drIono2(n) = jatGPSIonoModel((r1(:,n)+v1(:,n))*1000,...
            (r2(:,n)+v2(:,n))*1000, epoch_GPS+(t(n)+1)/86400)/1000;
    end
    
    drdotIono = drIono2 - drIono; % change in Ionospheric effect
elseif useIono    
    %Approximated based on a finite difference of the range error.
    %r1 and r2 states must be in an ECI.
    drIono = IonoModel(t,r1,r2,options); 
    drIono2 = IonoModel(t+1,r1+v1,r2+v2,options);
    drdotIono = drIono2 - drIono; % change in Ionospheric effect
end

%% Calculate Measurement Errors Due to Charged Particles above the Ionosphere
if useChPart 
    %Approximated based on a finite difference of the range error.
    %r1 states must be in an ECI.
    drdotChPart = ChPartModel(t,r1,r2,options);
    drdotChPart2 = ChPartModel(t+1,r1+v1,r2+v2,options); 
    drdotChPart = drdotChPart2 - drdotChPart; % change in Charged Particle effect
end

rangeRate = rangeRate + drdotTropo + drdotIono + drdotChPart;
end