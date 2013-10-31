function varargout = LOSRangeCon(obj, t, r1, r2, options)

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
%   useIono             {true, false(default)} only for groundstats as x2
%   useTropo            {true, false(default)} only for groundstats as x2
%   useChargedParticle  {true, false(default)} only for groundstats as x2
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
%   Phillip Anderson    07/81/2013      Updated to support object-oriented
%                                       code

useGPSIono = getOdtbxOptions(options, 'useGPSIonosphere', false );
useIono    = getOdtbxOptions(options, 'useIonosphere', false );
useTropo   = getOdtbxOptions(options, 'useTroposphere', false );
useChPart  = getOdtbxOptions(options, 'useChargedParticle', false );

cons_iono_flag = any(strncmpi(obj.loc_cons.param,'ION',3));
cons_trop_flag = any(strncmpi(obj.loc_cons.param,'TRP',3));
% numiono=sum(strncmpi(loc_cons,'ION',3));

r      = r1 - r2 ;
range  = sqrt(sum(r.^2));

if cons_iono_flag || cons_trop_flag
    % assume order is 
    % H(1,:) = dr_iono
    % H(2,:) = dr_tropo
    ndr = cons_iono_flag + cons_trop_flag;
    H = nan(ndr,length(range));
end

if useIono % Ionospheric range errors (only valid for ground-based tracking)
    %r1 and r2 states must be in ECI.
    epoch             = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SignalFreq  = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    Ephem.SatCoords   = 'ECEF';
    Ephem.StationInfo = 'ECEF';
    init2fixed        = jatDCM('eci2ecef', epoch+t/86400);
    drIono            = zeros(1,length(range)); %initialize before loop    
    for n=1:length(range)
        Ephem.Epoch  = epoch+t(n)/86400; %UTC
        Ephem.SatPos = init2fixed(:,:,n)*r1(:,n)*1000;
        Ephem.StaPos = init2fixed(:,:,n)*r2(:,n)*1000;
        drIono(n)    = jatIonoDelayModel(Ephem)/1000;
    end
    range = range + drIono;
    if cons_iono_flag
        H(1,:)=drIono;
    end
end

if useTropo % Tropospheric range errors (only valid for ground-based tracking)
    %r1 and r2 states must be in ECI.
    f                 = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    c                 = JATConstant('c'); %speed of light in meters
    epoch             = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SatCoords   = 'ECEF';
    Ephem.StationInfo = 'ECEF';
    init2fixed        = jatDCM('eci2ecef', epoch+t/86400);
    drTropo           = zeros(1,length(range)); %initialize before loop    
    for n = 1:length(range)
        Ephem.Epoch          = epoch+t(n)/86400; %UTC
        Ephem.satPos         = init2fixed(:,:,n)*r1(:,n)*1000;
        Ephem.staPos         = init2fixed(:,:,n)*r2(:,n)*1000;    
        [azAngle,elAngle]    = jatStaAzEl(Ephem);
        % Check if the elevation is too negative.  Currently, this is -15
        % deg.  Consider making this a user option.
        min_elAngle = -15*pi/180;
        elAngle( elAngle < min_elAngle ) = min_elAngle; 
        drTropo(n)           = jatTropoModel(range(n)*1000, elAngle, c/f)/1000;
    end    
    range   = range + drTropo;
    
    if cons_trop_flag && cons_iono_flag
        H(2,:)=drTropo;
    elseif cons_trop_flag
        H(1,:)=drTropo;
    end
end

if useGPSIono % GPS Ionospheric range errors (only valid for GPS tracking)
    %r1 and r2 states must be in an ECI.
    epoch     = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    epoch_GPS = convertTime('GPS','UTC',epoch); %GPS Time
    drIono    = zeros(1,length(range)); %initialize before loop
    for n = 1:length(range)
        drIono(n) = jatGPSIonoModel(r1(:,n)*1000, r2(:,n)*1000, epoch_GPS+t(n)/86400)/1000;
    end
    range  = range + drIono;
end

if useChPart % Charged Particle range errors (only used for ground-based
             % tracking where the satellite is outside the magnetosphere)
    %r1 states must be in an ECI.
    epoch                     = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SignalFreq          = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    Ephem.Epoch               = epoch+t/86400; %UTC
    Ephem.satPos              = r1*1000;
    drChPart                  = jatChargedParticleModel(Ephem)/1000;
    drChPart(isnan(drChPart)) = 0;
    range                     = range + drChPart';
end

varargout{1}=range;
if cons_iono_flag || cons_trop_flag
    varargout{2}=H;
end

end

%% Validation Test
%
%function failed = losrange_validation_test()
%
%disp(' ')
%disp('Performing Test....')
%disp(' ')
%fprintf('%20s%28s%25s%20s\n','Pos 1 (km)','Pos 2 (km)','R-Expected (km)','R-Calculated (km)')
%fprintf('%s\n\n',char(ones(1,92)*'-'));
%
%tol = 1e-7;
%options = odtbxOptions('measurement');
%options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006'));
%options = setOdtbxOptions(options,'useGPSIonosphere',false);
%options = setOdtbxOptions(options,'useIonosphere',false);
%options = setOdtbxOptions(options,'useTroposphere',false);
%options = setOdtbxOptions(options,'useChargedParticle',false);
%
%t=(1:9)*60*60;
%e1 = [10000; 0; 0]; e2 = [0; 10000; 0]; e3 = [0; 0; 10000];
%r1 = [e1 e1 e1 e2 e2 e2 e3 e3 e3];
%r2 = [e1 e2 e3 e1 e2 e3 e1 e2 e3];    
%
%ExRange = [0 sqrt(2e8) sqrt(2e8) sqrt(2e8) 0 sqrt(2e8) sqrt(2e8) sqrt(2e8) 0];
%Range = GetRange(t, r1, r2, options);
%
%fprintf('%8.2f %8.2f %8.2f %10.2f %8.2f %8.2f %16.6f %16.6f\n',...
%    [r1; r2; ExRange; Range]);
%
%passed = tol > max( abs( ExRange - Range ) );
%failed = ~passed;
%if failed
%    disp(' ')
%    disp('Test Failed!')
%else
%    disp(' ')
%    disp('Test Passed.')
%end
%end
