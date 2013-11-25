function [ drTropo ] = TropoModel_basic( t,r1,r2,options )
% Calculate the signal delay (in km) caused by a signal traveling 
% through the Troposphere. This function will only be accurate if r2 is a
% ground station.
%
% [ drTropo ] = TropoModel_basic( t,r1,r2,options )
%
% INPUTS
%   t         (1xN)	Times corresponding to r1 (secs)
%   r1        (3xN)	User spacecraft position (km) in ECI
%   r2        (3xN) Tracking spacecraft/ground station position (km) ECI
%   options   (1x1) Data structure
% OUTPUTS
%   drTropo   (1xN) Range error in meters due to Troposphere
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
%   Jason Schmidt       11/7/2013       original TropoModel_basic

    % If r2 is not a ground station, set drTropo to 0
    Re = JATConstant('rEarth')/1000; % Radius of the earth in km
    if all((sqrt((sum(r2.^2)))-Re) > 5)
        drTropo = zeros(1,numel(t));
        return
    end

    %r1 and r2 states must be in ECI.
    r      = r1 - r2 ;
    range  = sqrt(sum(r.^2));
    f                 = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    c                 = JATConstant('c'); %speed of light in meters
    epoch             = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SatCoords   = 'ECI';
    Ephem.StationInfo = 'ECEF';
    init2fixed        = jatDCM('eci2ecef', epoch+t/86400);
    drTropo           = zeros(1,numel(t)); %initialize before loop    
    for n = 1:numel(t)
        Ephem.Epoch          = epoch+t(n)/86400; %UTC
        Ephem.satPos         = r1(:,n)*1000;
        Ephem.staPos         = init2fixed(:,:,n)*r2(:,n)*1000;    
        [~,elAngle]    = jatStaAzEl(Ephem);
        % Check if the elevation is too negative.  Currently, this is -15
        % deg.  Consider making this a user option.
        min_elAngle = -15*pi/180;
        elAngle( elAngle < min_elAngle ) = min_elAngle; 
        drTropo(n)           = jatTropoModel(range(n)*1000, elAngle, c/f)/1000;
    end    
    
end

