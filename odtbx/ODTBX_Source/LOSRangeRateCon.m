function varargout = LOSRangeRateCon(obj,t, r1, v1, r2, v2, options)

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
%   useIono             {true, false(default)} only for groundstats as x2
%   useTropo            {true, false(default)} only for groundstats as x2
%   useChargedParticle  {true, false(default)} only for groundstats as x2
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
% keyword: measurement
% See also LOSRANGE, LOSDOPPLER, RRDOT, RRDOTLT

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

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    rangeRate = losrangerate_validation_test();  
else
    if nargout==2
        [rangeRate,H]=GetRangeRate(obj,t,r1,v1,r2,v2,options);
        varargout{2}=H;
    elseif nargout==1
        rangeRate=GetRangeRate(obj,t,r1,v1,r2,v2,options);
    end
    varargout{1}=rangeRate;
end
end


%% Main function
function varargout = GetRangeRate(obj, t, r1, v1, r2, v2, options)

useGPSIono = getOdtbxOptions(options, 'useGPSIonosphere', false );
useIono    = getOdtbxOptions(options, 'useIonosphere', false );
useTropo   = getOdtbxOptions(options, 'useTroposphere', false );
useChPart  = getOdtbxOptions(options, 'useChargedParticle', false );
% solve        = getOdtbxOptions(options, 'solvefor',[]);
% dyn_cons     = getOdtbxOptions(options, 'dynamicConsider',[]);
% loc_cons     = getOdtbxOptions(options, 'localConsider',[]);

cons_iono_flag = any(strncmpi(obj.loc_cons.param,'ION',3));
cons_trop_flag = any(strncmpi(obj.loc_cons.param,'TRP',3));

c           = JATConstant('c')/1000; %Speed of light in km
r           = r1 - r2 ;
losVector   = unit( r );
rangeRate   = dot( losVector, (v1-v2) ) ./ (1 - dot(losVector,v2)/c );

if cons_iono_flag || cons_trop_flag
    % assume order is 
    % H(1,:) = dr_iono
    % H(2,:) = dr_tropo
    ndr = cons_iono_flag + cons_trop_flag;
    H = nan(ndr,length(rangeRate));
end

if useIono %Ionospheric range rate errors (only valid for ground-based tracking)    
    %Approximated based on a finite difference of the range error.
    %r1 and r2 states must be in an ECI.

    epoch = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SignalFreq = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    Ephem.SatCoords = 'ECEF';
    Ephem.StationInfo = 'ECEF';
    init2fixed  = jatDCM('eci2ecef', epoch+t/86400);
    drIono = zeros(1,length(t)); %initialize before loop    
    for n=1:length(t)
        Ephem.Epoch = epoch+t(n)/86400; %UTC
        Ephem.SatPos = init2fixed(:,:,n)*r1(:,n)*1000;
        Ephem.StaPos = init2fixed(:,:,n)*r2(:,n)*1000;
        drIono(n) = jatIonoDelayModel(Ephem)/1000;
    end
    
    init2fixed  = jatDCM('eci2ecef', epoch+t/86400+1);
    drIono2 = zeros(1,length(t)); %initialize before loop    
    for n=1:length(t)
        Ephem.Epoch = epoch+(t(n)+1)/86400; %UTC
        Ephem.satPos = init2fixed(:,:,n)*(r1(:,n)+v1(:,n))*1000;
        Ephem.staPos = init2fixed(:,:,n)*(r2(:,n)+v2(:,n))*1000;
        drIono2(n) = jatIonoDelayModel(Ephem)/1000;
    end    
    
    drdotIono = drIono2 - drIono; % change in Ionospheric effect
    rangeRate = rangeRate + drdotIono;
    if cons_iono_flag
        H(1,:)=drdotIono;
    end
end

if useTropo % Tropospheric range rate errors (only valid for ground-based tracking)
    %Approximated based on a finite difference of the range error.
    %r1 and r2 states must be in ECI.
    range             = sqrt(sum(r.^2));
    f                 = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    epoch             = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SatCoords   = 'ECI';
    Ephem.StationInfo = 'ECEF';
    init2fixed        = jatDCM('eci2ecef', epoch+t/86400);
    drTropo           = zeros(1,length(range)); %initialize before loop    
    for n = 1:length(range)
        Ephem.Epoch          = epoch+t(n)/86400; %UTC
        Ephem.satPos         = r1(:,n)*1000;
        Ephem.staPos         = init2fixed(:,:,n)*r2(:,n)*1000;    
        [azAngle,elAngle]    = jatStaAzEl(Ephem);
        elAngle( elAngle<0 ) = NaN; % can't have a negative elevation
        drTropo(n)           = jatTropoModel(range(n)*1000, elAngle, c/f)/1000;
    end        
    
    init2fixed = jatDCM('eci2ecef', epoch+t/86400+1);
    drTropo2   = zeros(1,length(range)); %initialize before loop    
    for n = 1:length(range)
        Ephem.Epoch          = epoch+t(n)/86400; %UTC
        Ephem.satPos         = r1(:,n)*1000;
        Ephem.staPos         = init2fixed(:,:,n)*r2(:,n)*1000;    
        [azAngle,elAngle]    = jatStaAzEl(Ephem);
        elAngle( elAngle<0 ) = NaN; % can't have a negative elevation
        drTropo2(n)          = jatTropoModel(range(n)*1000, elAngle, c/f)/1000;
    end  

    drdotTropo = drTropo2 - drTropo; % change in Tropospheric effect
    rangeRate  = rangeRate + drdotTropo;
    if cons_trop_flag && cons_iono_flag
        H(2,:)=drdotTropo;
    elseif cons_trop_flag
        H(1,:)=drdotTropo;
    end
end

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
    rangeRate = rangeRate + drdotIono;
end

if useChPart % Charged Particle range rate errors (only used for ground-based
             % tracking where the satellite is outside the magnetosphere)
    %r1 states must be in an ECI.
    epoch = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SignalFreq = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    Ephem.Epoch = epoch+t/86400; %UTC
    Ephem.satPos = r1*1000;
    drChPart = jatChargedParticleModel(Ephem)/1000;
    drChPart(isnan(drChPart)) = 0;
    
    Ephem.Epoch = epoch+(t+1)/86400; %UTC
    Ephem.satPos = (r1+v1)*1000;
    drChPart2 = jatChargedParticleModel(Ephem)/1000;
    drChPart2(isnan(drChPart2)) = 0;
    
    drdotChPart = drChPart2 - drChPart; % change in Charged Particle effect
    rangeRate = rangeRate + drdotChPart';
end

varargout{1}=rangeRate;
if cons_iono_flag || cons_trop_flag
    varargout{2}=H;
end
end
%% Validation Test

function failed = losrangerate_validation_test()

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%12s%17s%19s%17s%20s%20s\n','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (km/s)','Calculated (km/s)')
fprintf('%s\n\n',char(ones(1,105)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006'));
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);

t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [e1 e1 e1 e1 e1 e1 e1 e1 e1];
r2 = [e2 e2 e2 e2 e2 e2 e2 e2 e2];    
v1 = [e1 e1 e1 e2 e2 e2 e3 e3 e3];
v2 = [e1 e2 e3 e1 e2 e3 e1 e2 e3];

ExRangeRate = [0 1.41421022674001 0.707106781186547 -1.41421689802191 ...
    0 -0.707106781186547 -0.707108449010957 0.707105113370005 0];
RangeRate = GetRangeRate(t, r1, v1, r2, v2, options);

fprintf('%4.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %16.6f %16.6f\n',...
    [r1; r2; v1; v2; ExRangeRate; RangeRate]);

passed = tol > max( abs( ExRangeRate - RangeRate ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
