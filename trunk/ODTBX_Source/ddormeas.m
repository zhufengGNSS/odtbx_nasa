function [y,H,R] = ddormeas(t,x,options)
% Makes delta-differenced oneway range measurements.
%
% [y,H,R] = DDORMEAS(tspan,x,options) creates delta-differenced oneway
% range measurements based on the information in OPTIONS. The measurement
% type that this function returns is an interferometric range difference
% with that of a "quasar". Rather than use an actual quasar database, the
% "quasar" used in this model is assumed to be directly behind the
% satellite at t=0 and fixed in ECI space for all times after.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch)
%      x            (6xN)   spacecraft state [position;velocity] (km)
%      options      (1x1)   data structure, see below
%
%   OUTPUTS
%      y            (MxN)   measurements
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each groundstation combination are grouped
% together in rows. Thus, using the default options settings and passing in
% 3 ground stations, the output y will look like:
%   y = [   ddor_gs1&2(t1)      ddor_gs1&2(t2)...;
%           ddor_gs1&3(t1)      ddor_gs1&3(t2)...;
%           ddor_gs2&3(t1)      ddor_gs2&3(t2) ]
% If ddor rate is included, then the output looks like this:
%   y = [   ddor_gs1&2(t1)      ddor_gs1&2(t2)...;
%           ddorRate_gs1&2(t1)  ddorRate_gs1&2(t2)...;
%           ddor_gs1&3(t1)      ddor_gs1&3(t2)...;
%           ddorRate_gs1&3(t1)  ddorRate_gs1&3(t2)...;
%           ddor_gs2&3(t1)      ddor_gs2&3(t2)...;
%           ddorRate_gs2&3(t1)  ddorRate_gs2&3(t2) ]
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
% For station locations based on information in the NASA Directory Of
% Station Locations (NDOSL), set these parameters:
%   PARAMETER           VALID VALUES         NOTES
%   gsList                JAT object         See CREATEGROUNDSTATIONLIST
%   gsID                  See NDOSL file     Cell array of ground station
%                                             IDs in the  NDOSL
%
% For faster perfomance (or for unlisted ground stations), predefine the 
% ECEF station locations as columns of a matrix using this parameter:
%   PARAMETER           VALID VALUES         NOTES
%   gsECEF                [(3xNs) matrix]    ECEF locations of ground
%                                             stations (km). Ns is the 
%                                             number of stations.
%
% The following parameters apply to both location input types:
%   PARAMETER           VALID VALUES         NOTES
%   epoch             datenum                UTC time associated with 
%                                             start of simulation
%   rSigma            {(1xM),ones(1,M)*1e-3) Measurement uncertainty for
%                                             each measurement type for
%                                             each ground station
%   ddor              structure(km, radians) Structure of all TDRS data
%       Optional ddor rate calculation:
%         ddor.rate   {true, false(default)} Return ddor rate measurement
%       Optional Schedule Data:
%         ddor.Sched  [(ncx4) matrix]        Tracking Schedule. 
%
%      Schedule restricts the ddor measurement model to only provide  
%      measurements for specific ground stations during specific time  
%      intervals. The format for each row is:
%      [ gs_index1 gs_index2 start_time stop_time ].
%      The gs_index corresponds to index of the ground stations in gsID or
%      the column number in gsECEF. The start_time and stop_time must be in
%      seconds from epoch. The matrix can be any length (nc = number of
%      contacts)
%    Example:
%     TrackSched = {
%         1 2 '01 Jul 2007 14:00:00' '01 Jul 2007 18:00:00'
%         1 3 '02 Jul 2007 07:30:00' '02 Jul 2007 11:00:00'
%         2 1 '02 Jul 2007 14:00:00' '02 Jul 2007 18:00:00'
%         1 3 '03 Jul 2007 07:30:00' '03 Jul 2007 11:00:00'
%         1 2 '03 Jul 2007 14:00:00' '03 Jul 2007 18:00:00'
%         3 1 '04 Jul 2007 07:30:00' '04 Jul 2007 11:00:00'
%         };
%     ddor.Sched     =cell2mat(TrackSched(:,1:2)); %ground station numbers
%     ddor.Sched(:,3)=(datenum(TrackSched(:,3))-epoch)*86400;%start(epsecs)
%     ddor.Sched(:,4)=(datenum(TrackSched(:,4))-epoch)*86400;%end(epsecs)
%
% The options parameters associated with each of the desired measurement
% types are passed to the appropriate function.
%
% The ground stations to use are input as part of OPTIONS. The
% groundstation list must also be input as part of OPTIONS. This list can
% be created by CREATEGROUNDSTATIONLIST.
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
% See also ODTBXOPTIONS, GSMEAS, CREATEGROUNDSTATIONLIST
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
%   Kevin Berry         11/09/2009      Original
%   Kevin Berry         03/10/2010      Added DDOR-Rate Option
%                                       Added the option to input station
%                                         locations in ECEF instead of 
%                                         using the NDOSL list at every 
%                                         time step

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    y = ddormeas_test();
else
    [y,H,R] = Get_ddormeas(t,x,options);
end
end

%% Main function
function [y,H,R] = Get_ddormeas(t,x,options)

%% Get values from options
epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
gsID         = getOdtbxOptions(options, 'gsID', [] );
gsList       = getOdtbxOptions(options, 'gsList', []);
gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
ddor         = getOdtbxOptions(options, 'ddor',[]); %ddor structure
if ~isfield(ddor,'rate'); 
    ddor.rate = false; 
end
if isempty(gsECEF)
    if( isempty(gsList) && ~isempty(gsID) ); 
        gsList = createGroundStationList(); 
    end
    gsECEF = zeros(3,length(gsID));
    for n=1:length(gsID)
            gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
    end
end

numtypes       = 1+ddor.rate;
numtrackpairs  = sum(1:size(gsECEF,2)-1);
M              = numtypes*numtrackpairs;
N              = length(t);
if size(t,1)==N, t=t'; end

if isnan(epoch); error('An epoch must be set in the options structure.'); end
if ~isfield(ddor,'Sched')
    %If no schedule is given, all ground station combinations will be valid
    %for all times.
    ddor.Sched = zeros(numtrackpairs,4);
    ind = 0;
    for n=1:size(gsECEF,2)
        for m=n+1:size(gsECEF,2)
            ind = ind+1;
            ddor.Sched(ind,1:4) = [n m -inf inf];
        end
    end
end

%get sigma out of the options
sigmaDefault = ones(1,size(gsECEF,2)*numtypes)*1e-3; %default sigma
gsSigmas     = getOdtbxOptions(options, 'rSigma', sigmaDefault );

%% Persistent variables for file data caching
persistent QuasarCache;
UnitOptions = odtbxOptions('measurement');
UnitOptions = setOdtbxOptions(UnitOptions,'useUnit',     true);
UnitOptions = setOdtbxOptions(UnitOptions,'useRange',    false);
UnitOptions = setOdtbxOptions(UnitOptions,'useRangeRate',false);
UnitOptions = setOdtbxOptions(UnitOptions,'useDoppler',  false);


%% If required, initialize any persistent variables
if isempty(QuasarCache)
    QuasarCache = dataCache('create');
end

% Get the quasar for this satellite
quasar_eci = dataCache('get',QuasarCache,'quasar');
if isempty(quasar_eci)
    % Not already in the cache.  Set it now and save it in the cache.
    gs = getGSstates(gsECEF,epoch,t(1)); %(see subfunction below)
    for n=1:size(gsECEF,2)
        quasar_eci(:,n) = unit( x(1:3,1)-gs(1:3,1,n) );
    end
    QuasarCache = dataCache('add', QuasarCache, 'quasar', quasar_eci);
end

%% Loop over all ground station combinations
y             = nan(M,N);
H             = zeros(M,6,N);
sigma_squared = zeros(M,M);
ind           = 0;
for n=1:size(gsECEF,2)
    for m=n+1:size(gsECEF,2)
        ind = ind+1;
        ind_nm = n==ddor.Sched(:,1) & m==ddor.Sched(:,2);
        ind_mn = m==ddor.Sched(:,1) & n==ddor.Sched(:,2);
        gSch = ddor.Sched(ind_nm | ind_mn,3:4); %gs1 to gs2 = gs2 to gs1
        tind = [];
        for nm=1:size(gSch,1)
            tind = union(tind, find( gSch(nm,1)<=t & t<=gSch(nm,2) ));
        end
        ltind = length(tind);
        
        if ltind~=0
            t1  = t(tind);
            gs  = getGSstates(gsECEF(:,[n m]),epoch,t1); %(see subfunction below)
            gs1 = gs(:,:,1);
            gs2 = gs(:,:,2);
            b   = gs1(1:3,:) - gs2(1:3,:);
                      
            [u,H_u] = rrdotang(t1,x(:,tind),gs1,UnitOptions);

            iRange_user = dot(b,u);
            iRange_quas = b'*quasar_eci(:,n);
            
            ind1 = 1 + numtypes*(ind-1);
            y(ind1, tind)    = iRange_user(:) - iRange_quas(:);
            
            for nm=1:ltind
                H(ind1, :, tind(nm)) = b(:,nm)'*H_u(:,:,nm);
            end
            
            if ddor.rate
                bdot   = gs1(4:6,:) - gs2(4:6,:);
                rho    = sqrt(sum((x(1:3,tind)-gs1(1:3,:)).^2));
                vdiff  = x(4:6,tind) - gs1(4:6,:);
                udot   = zeros(3,ltind);
                for nm=1:ltind
                    udot(:,nm) = H_u(1:3,1:3,nm)*vdiff(:,nm);
                end
                
                iRR_1 = dot(bdot,u);
                iRR_2 = bdot'*quasar_eci(:,n);
                iRR_3 = dot( b , udot );
                ind2  = numtypes*ind;
                y(ind2, tind) = iRR_1(:) - iRR_2(:) + iRR_3(:);
                
                for nm=1:ltind
                    Htemp = H_u(1:3,1:3,nm)*(vdiff(:,nm)*u(:,nm)')...
                        + H_u(1:3,1:3,nm)*(u(:,nm)'*vdiff(:,nm))...
                        +(u(:,nm)*vdiff(:,nm)')*H_u(1:3,1:3,nm);
                    H(ind2, 1:3, tind(nm)) = bdot(:,nm)'*H_u(1:3,1:3,nm) - b(:,nm)'*Htemp./rho(nm);
                    H(ind2, 4:6, tind(nm)) = b(:,nm)'*H_u(1:3,1:3,nm);
                end                
            end                
        end        
    end
end

% Extra Output
if nargout > 2,
    if numtypes == 1
        ind = 0;
        for n=1:size(gsECEF,2)
            for m=n+1:size(gsECEF,2)
                ind = ind+1;
                sigma_squared(ind,ind) = gsSigmas(n)^2+gsSigmas(m)^2;
            end
        end
    elseif numtypes == 2
        ind = [-1 0];
        for n=1:size(gsECEF,2)
            for m=n+1:size(gsECEF,2)
                ind = ind+2;
                sigma_squared(ind,ind) = diag(gsSigmas(2*n-1:2*n).^2+gsSigmas(2*m-1:2*m).^2);
            end
        end
    else
        error('A 3rd ddor measurement type has not been defined here')
    end
    R = repmat(sigma_squared,[1,1,length(t)]);
end

end %end main function

%% Subfunction for getting ECI states of ground stations at times t
function gx = getGSstates(gsECEF,epoch,t)

D = jatDCM('ecef2eci', epoch+t/86400);
w = [0;0;JATConstant('wEarth')];

gx = zeros(6,length(t),size(gsECEF,2));
for nt = 1:length(t);
    M          = rotransf(-D(:,:,nt)*w,D(:,:,nt));
    gx(:,nt,:) = M(1:6,1:6)*[gsECEF;zeros(size(gsECEF))];
end

end %end function

%% Validation Test

function failed = ddormeas_test()
disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;
epoch    = datenum('01 Jul 2007 14:00:00.000');
x0 = [1.5e8 0 0 0 5e-2 0]';% km & km/sec
P0 = diag([10 10 10 .1 .1 .1].^2); %#ok<NASGU> %initial covariance matrix (km^2 & km^2/sec^2)

% Set the Tracking Stations
gsID = {'DS15',... Goldstone 34 Meter
        'DS45',... Canberra 34 Meter
        'DS65',... Madrid 34 Meter
        }; %The NASA ID for each ground station
gsSigma = [1e-6 1e-9... (km range & km/s rangerate sigma) Goldstone 34 Meter
           1e-6 1e-9... (km range & km/s rangerate sigma) Canberra 34 Meter
           1e-6 1e-9... (km range & km/s rangerate sigma) Madrid 34 Meter
           ]; %The measurement sigmas for each ground station

% Set the measurement options
gsList  = createGroundStationList; %generates a list of all of the ground stations
gsECEF    = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end
ddor.rate = true;

% Initialize options structure. In this case it is a measurement structure
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
measOptions = setOdtbxOptions(measOptions,'epoch',epoch); %datenum format
measOptions = setOdtbxOptions(measOptions,'rSigma',gsSigma);
measOptions = setOdtbxOptions(measOptions,'ddor',ddor);

% Set propagator information
tspan  = 0:6000:86400;
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', epoch); %datenum format
jatWorld = createJATWorld(jOptions); %creates a java object that stores the
% default information for propagating Earth-centric orbits using JAT
% force models.
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'ValidationCase', 0);
[t,x] = integ(@jatForces_km,tspan,x0,eOpts,jatWorld);

% Set Expected Results
y_expected = [
                         0    -1.34971696639199e-005     9.09494701772928e-013     6.16669884865372e-006                         0     2.40550487673035e-007
       -0.0881998013637713    -1.37598456273745e-005      -0.00842719790489355    -8.31648270929288e-006        0.0229141913841886      5.5283473699666e-006
        -0.143038607860944    -3.28935742790178e-006       -0.0812128905763529    -1.37676275811608e-005        0.0389013225230883    -2.47178109346298e-006
        -0.123861529403257       9.023266792148e-006        -0.146281489874127    -5.92356581078859e-006       -0.0262979671406356    -1.99124986330633e-005
       -0.0527340503595042     1.25798347026248e-005        -0.137660174973007     9.01487714432817e-006        -0.192994662347701    -3.38448750284676e-005
      0.000338087622367311     3.11895617408728e-006       -0.0486858088784174     1.87682090613836e-005        -0.399891655461943    -3.18757302564991e-005
       -0.0305411990620996    -1.37288855254595e-005        0.0577744859347149     1.39978103044686e-005        -0.538331014643518    -1.17829225805727e-005
        -0.155215144309295    -2.60527629277344e-005        0.0910142582843037    -4.53809954024803e-006        -0.525460187659519     1.60023425420065e-005
        -0.313688684326962    -2.39666815885454e-005      -0.00258448906970443    -2.59160925819613e-005        -0.364185146097952     3.50332746428538e-005
        -0.411790444298276    -6.73245689288389e-006        -0.197421004382704    -3.62513974665533e-005        -0.146098953617184     3.40167136666419e-005
         -0.38242396581245     1.63858949784293e-005        -0.399902802635552    -2.81583801182888e-005       0.00637499022741395     1.46100838588964e-005
        -0.230134243601242     3.21928975547009e-005        -0.506190530161803    -5.81548299080938e-006        0.0191982840624405    -9.71786149971974e-006
       -0.0282301305951478     3.23080950674613e-005        -0.466922345807461     1.79141247614341e-005       -0.0875235588455325    -2.29954392734217e-005
         0.129034686206069     1.84616987743615e-005        -0.314673366203351     3.02275675136154e-005        -0.218558604502505    -1.76625191324119e-005
         0.185424678748859     7.02347231531437e-007        -0.136478502226964     2.67860760906881e-005        -0.273003076374152     5.89128784256221e-007]';

% Display Expected Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Expected DDOR Measurements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-5s %14s %14s %14s %14s %14s %14s\n','Time','gs1 - gs2','gs1 - gs2','gs1 - gs3','gs1 - gs3','gs2 - gs3','gs2 - gs3');
fprintf('%3s %13s %15s %13s %15s %13s %15s\n','(s)','(km)','(km/s)','(km)','(km/s)','(km)','(km/s)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:length(y_expected)
    fprintf('%-6i %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n',t(n),y_expected(:,n))
end
disp(' ')
disp(' ')

% Run tdrssmeas
y = ddormeas(t,x,measOptions);

% Display Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Calculated DDOR Measurements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-5s %14s %14s %14s %14s %14s %14s\n','Time','gs1 - gs2','gs1 - gs2','gs1 - gs3','gs1 - gs3','gs2 - gs3','gs2 - gs3');
fprintf('%3s %13s %15s %13s %15s %13s %15s\n','(s)','(km)','(km/s)','(km)','(km/s)','(km)','(km/s)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:length(y)
    fprintf('%-6i %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n',t(n),y(:,n))
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