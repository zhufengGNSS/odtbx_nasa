function [ range ] = LOSRange( t, r1, r2, options )
% LOSRANGE  Line of sight range measurement between two objects
%
%   range = LOSRange(r1, r2, options) returns the range measurement between
% two vehicles. This function can also include the troposphere and
% ionosphere effects if requested.  Note that a minimum elevation of -15
% degrees is applied and is held at this value if elevation is less than
% -15 deg.
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
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t         (1xN)	Times corresponding to r1 (secs)
%      r1        (3xN)	User spacecraft position (km)
%      r2        (3xN)  Tracking spacecraft/ground station position (km)
%      options   (1x1)  Data structure
%
% OUTPUTS
%      range     (1xN)  range (km)
%
% VALIDATION/REGRESSION TEST
%
%  These have been moved to LOSRange_test.m in the regression testing
%  framework to conform with the new testing format.  
%
% keyword: measurement
% See also LOSRANGERATE, LOSDOPPLER, RRDOT, RRDOTLT
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
%   Bo Naasz            12/03/2005      Original xlmeas.m
%   Derek Surka         07/20/2007      Modified xlmeas to work with OD
%                                       Toolbox
%   Derek Surka         08/22/2007      Modified input parameters, added
%                                       functionality
%   Kevin Berry         05/08/2008      Replaced with a simpler r1-r2 model
%   Kevin Berry         05/22/2008      Added JAT Iono and Tropo models
%   Kevin Berry         09/05/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         09/08/2008      Added Charged Particle Model
%   Kevin Berry         12/23/2008      Replaced elevation.m with
%                                       jatStaAzEl.m
%   Kevin Berry         06/25/2009      Fixed time scale discrepancy in 
%                                       GPSIono
%   Ravi Mathur         08/28/2012      Extracted regression test
%   Jason Schmidt       11/6/2013       Encapsulated basic Tropo, Iono and
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
    IonoModel = @IonoModel_basic;
end
    
if isa(useTropo,'function_handle')
    TropoModel = useTropo;
    useTropo = true;
elseif useTropo==true
    TropoModel = @TropoModel_basic;
end

if isa(useChPart,'function_handle')
    ChPartModel = useChPart;
    useChPart = true;
elseif useChPart==true
    ChPartModel = @ChPartModel_basic;
end


%% Calculate Geometric Range and Initialize Values
r      = r1 - r2 ;
range  = sqrt(sum(r.^2));

drTropo = 0;
drIono = 0;
drChPart = 0;

%% Calculate Measurement Errors Due to Troposphere
if useTropo
    drTropo = TropoModel(t,r1,r2,options);
end

%% Calculate Measurement Errors Due to Ionosphere
if useGPSIono
    %r1 and r2 states must be in an ECI.
    epoch     = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    epoch_GPS = convertTime('GPS','UTC',epoch); %GPS Time
    drIono    = zeros(1,length(range)); %initialize before loop
    for n = 1:length(range)
        drIono(n) = jatGPSIonoModel(r1(:,n)*1000, r2(:,n)*1000, epoch_GPS+t(n)/86400)/1000;
    end
elseif useIono
    drIono = IonoModel(t,r1,r2,options);
end

%% Calculate Measurement Errors Due to Charged Particles above the Ionosphere
if useChPart
    drChPart = ChPartModel(t,r1,r2,options);
end

range = range + drTropo + drIono + drChPart;

end

