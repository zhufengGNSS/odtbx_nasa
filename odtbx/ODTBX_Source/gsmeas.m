function [y,H,R] = gsmeas(t,x,options)
% GSMEAS  Makes ground station based measurements.
%
% [y,H,R] = GSMEAS(tspan,x,options) creates ground station measurements
% based on the information in OPTIONS. The measurement types that can be 
% returned are LOSRANGE, LOSRANGERATE, LOSDOPPLER. See the OD Toolbox 
% function of the same name for details of each measurement type.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch)
%      x            (6xN)   ECI J2000 spacecraft state [pos;vel] (km)
%      options      (1x1)   data structure, see below
%
%   OUTPUTS
%      y            (MxN)   measurements
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each groundstation are grouped together in
% rows. Thus, using the default options settings and passing in 3 ground
% stations, the output y will look like:
%   y = [   range_gs1(t1)       range_gs1(t2)...;
%           rangeRate_gs1(t1)   rangeRate_gs1(t2)...;
%           range_gs2(t1)       range_gs2(t2)...;
%           rangeRate_gs2(t1)   rangeRate_gs2(t2)...;
%           range_gs3(t1)       range_gs3(t2)...;
%           rangeRate_gs3(t1)   rangeRate_gs3(t2)  ]
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
%   epoch                 datenum            UTC time associated with 
%                                             start of simulation
%   gsElevationConstraint degs               Spacecraft elevation must be
%                                             above this value to return
%                                             measurement
%   useRange          {true(default), false} Return range measurement
%   rangeType         {'1way','1wayFWD','1wayRTN','2way'(default)} 
%                                             Note: 1way=1wayFWD
%   useRangeRate      {true(default), false} Return range rate measurement
%   useDoppler        {true, false(default)} Return doppler measurement
%   useUnit           {true, false(default)} Return unit vector measuerment
%   frequencyTransmit {scalar>0, 1.57542e9}  Hz, Only used for Tropo and 
%                                             for Doppler
%   rSigma            {(1xM),ones(1,M)*1e-3) Measurement covariance
%   useLightTime      {true, false(default)} Include light time delay
%   useGPSIonosphere  {true, false(default)} Includes GPS Ionospheric delay
%   useIonosphere     {true, false(default)} Includes Ionospheric delay
%   useTroposphere    {true, false(default)} Includes Tropospheric delay
%   Schedule          [(ncx3) matrix]        Ground Tracking Schedule. 
%      Schedule restricts the gsmeas measurement model to only provide  
%      measurements for specific ground stations during specific time  
%      intervals. Theformat for each row is:
%      [ gs_index start_time stop_time ].
%      The gs_index corresponds to index of the ground stations in gsID or 
%      gsECEF. The start_time and stop_time must be in seconds from epoch.
%      The matrix can be any length (nc = number of contacts)
%      Example: 
%       TrackSched = {1 '28-Jan-2010 06:56:35' '28-Jan-2010 12:56:35'
%                     2 '28-Jan-2010 12:56:35' '28-Jan-2010 19:56:35'
%                     1 '28-Jan-2010 19:56:35' '29-Jan-2010 02:56:35'};
%       Sched = cell2mat(TrackSched(:,1)); %ground station numbers
%       Sched(:,2) = (datenum(TrackSched(:,2))-epoch)*86400; %start(epsecs)
%       Sched(:,3) = (datenum(TrackSched(:,3))-epoch)*86400; %end(epsecs)
%
% The options parameters associated with each of the desired measurement
% types are passed to the appropriate function.
%
% The ground stations to use are input as part of OPTIONS. The
% groundstation list must also be input as part of OPTIONS. This list can
% be created by CREATEGROUNDSTATIONLIST.
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, ODTBXOPTIONS,
% CREATEGROUNDSTATIONLIST
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
%   Kevin Berry         09/07/2007      Original
%   Derek Surka         09/17/2007      Modified to use options input and 
%                                        elevation constraint
%   Kevin Berry         04/10/2008      Added Doppler functionality
%   Kevin Berry         05/21/2008      Simplified with my new rrdot and
%                                        rrdotlt functions
%   Allen Brown         02/25/2009      Updated documentation
%   Kevin Berry         06/25/2009      Added time scale comments
%   Kevin Berry         07/01/2009      Added angle based measurements
%   Kevin Berry         09/03/2009      Changed angle measurement to Unit
%   Kevin Berry         02/02/2010      Added the option to input station
%                                         locations in ECEF instead of 
%                                         using the NDOSL list at every 
%                                         time step
%   Russell Carpenter   02/11/2011      Added useAngles option

%% Get values from options
gsID         = getOdtbxOptions(options, 'gsID', [] );
gsList       = getOdtbxOptions(options, 'gsList', []);
gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
elMin        = getOdtbxOptions(options, 'gsElevationConstraint', 10)*pi/180; % convert from degs to rads
uselt        = getOdtbxOptions(options, 'useLightTime', false);
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
useUnit      = getOdtbxOptions(options, 'useUnit', false );
useAngles    = getOdtbxOptions(options, 'useAngles', false );
Sched        = getOdtbxOptions(options, 'Schedule',[]); %Tracking Schedule
numtypes     = useRange + useRangeRate + useDoppler+3*useUnit+2*useAngles;

if isempty(gsECEF)
    if( isempty(gsList) && ~isempty(gsID) ); gsList = createGroundStationList(); end
    gsECEF = zeros(3,length(gsID));
    for n=1:length(gsID)
            gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
    end
end
M            = size(gsECEF,2) * numtypes;
N            = length(t);
if size(t,1)==N, t=t'; end

if isnan(epoch); error('An epoch must be set in the options structure.'); end


%% call rrdot with all the options passing straight through
y    = nan(M,N);
H    = zeros(M,6,N);  
if uselt
    %x2 needs states before and after each state in x1 at time steps
    %comparable to the light time delay in order for the interpolation
    %within the lightTimeCorrection function to have a high level of
    %accuracy.
    c    = JATConstant('c')/1000;
    ltDT = sqrt(sum(x(1:3,:).^2))/c; %
    t2   =  unique([t-ltDT, t, t+ltDT]);

    % Get the ground station positions at times t2 (see subfunction below)
    x2 = getGSstates(gsECEF,epoch,t2);

    % Run rrdotlt for each ground station
    for n=1:size(gsECEF,2)
        % Check schedule for times when tracking is done
        if ~isempty(Sched)
            gSch = Sched(n==Sched(:,1),2:3);
            tind = [];
            for m=1:size(gSch,1)
                tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
            end
        else
            tind = 1:length(t);
        end
        t1  = t(tind);
        x1  = x(:,tind);
        
        if ~isempty(t1)

            [y1,H1,R,t2_lt,x2_lt] = rrdotlt(t1,x1,t2,x2(:,:,n),options);
            
            % apply the elevation constraint
            Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
            Ephem.SatCoords   = 'ECI';
            Ephem.Epoch       = epoch+t2_lt{1}/86400; %UTC
            Ephem.StationInfo = 'ECEF';
            Ephem.staPos      = gsECEF(:,n)*1000;
            [~,el]            = jatStaAzEl(Ephem);
            index0            = find( el < elMin );
            if length(x2_lt)==2 %then it was a 2way measurement
                Ephem.Epoch = epoch+t2_lt{2}/86400; %UTC
                [~,el]      = jatStaAzEl(Ephem);
                index1      = find( el < elMin );
                index0      = union(index0,index1);
            end
            y1(:,index0) = NaN;

            % combine with results from previous stations
            indstart                     = 1 + numtypes*(n-1);
            indstop                      = numtypes*n;
            y(indstart:indstop, tind)    = y1;
            H(indstart:indstop, :, tind) = H1;
        end
    end
    clear R;
else
    % Get the ground station positions at times t (see subfunction below)
    gx = getGSstates(gsECEF,epoch,t);

    % Run rrdot for each ground station
    for n=1:size(gsECEF,2)
        % Check schedule for times when tracking is done
        if ~isempty(Sched)
            gSch = Sched(n==Sched(:,1),2:3);
            tind = [];
            for m=1:size(gSch,1)
                tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
            end
        else
            tind = 1:length(t);
        end
        if isempty(tind),continue,end
        t1 = t(tind);
        x1 = x(:,tind);
        x2 = gx(:,tind,n);

        [y1,H1] = rrdotang(t1,x1,x2,options);

        % apply the elevation constraint
        Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
        Ephem.SatCoords   = 'ECI';
        Ephem.Epoch       = epoch+t1/86400;%UTC
        Ephem.StationInfo = 'ECEF';
        Ephem.staPos      = gsECEF(:,n)*1000;
        [~,el]            = jatStaAzEl(Ephem);        
        y1(:,el<elMin)    = NaN;

        % combine with results from previous stations
        indstart                     = 1 + numtypes*(n-1);
        indstop                      = numtypes*n;
        y(indstart:indstop, tind)    = y1;
        H(indstart:indstop, :, tind) = H1;
    end
end

%% Set the measurement covariance output
if nargout > 2,
    %get sigma out of the options
    sigmaDefault = ones(1,size(y,1))*1e-3;
    sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
    if( length(sigma)~=length(sigmaDefault) )
        disp('WARNING: In gsmeas, length(rSigma) does not match total number of measurements');
        disp('   Setting rSigma = default (1e-3) for all measurements');
        sigma = sigmaDefault;
    end
    sigma           = diag(sigma);
    R = repmat(sigma.^2,[1,1,N]);
end


%% Subfunction for getting ECI states of ground stations at times t
function gx = getGSstates(gsECEF,epoch,t)

D = jatDCM('ecef2eci', epoch+t/86400);
w = [0;0;JATConstant('wEarth')];

gx = zeros(6,length(t),size(gsECEF,2));
for nt = 1:length(t);
    M          = rotransf(-D(:,:,nt)*w,D(:,:,nt));
    gx(:,nt,:) = M(1:6,1:6)*[gsECEF;zeros(size(gsECEF))];
end