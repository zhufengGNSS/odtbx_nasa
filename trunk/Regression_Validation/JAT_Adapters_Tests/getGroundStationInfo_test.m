function failed = getGroundStationInfo_test()
% Regression Test Case
% Function(s) getGroundStationInfo, createGroundStationList
%
% This function tests the ground station database functions using the DS16
% NDOSL data
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

failed  = 0;
tol     = 1e-10;

epoch           = datenum('June 1, 2011');
targetID        = 'DS16';
targetNASAID    = '1312';
targetLocation  = 'Goldstone, CA';
targetLatitude  = 35.3415393888889;
targetLongitude = 243.12635025;
targetAltitude  = 0.943977;
targetECEFPos   = [-2354.76332589963; -4646.78738368213; 3669.38335404832];
% A target's ECEF Velocity should be zero:
targetECEFVel   = [0.0; 0.0; 0.0];
targetECIPos    = [-3489.11010054197; 3864.47951575575; 3673.39438227322];
targetECIVel    = [-0.281804218157403; -0.254735171477489; 0.000319298187696472];

try
    gsList = createGroundStationList();
catch
    failed = 1;
    return
end

location        = getGroundStationInfo( gsList, targetID, 'location' );
nasaID          = getGroundStationInfo( gsList, targetID, 'nasa_id' );
latitude        = getGroundStationInfo( gsList, targetID, 'latitude' );
longitude       = getGroundStationInfo( gsList, targetID, 'longitude' );
altitude        = getGroundStationInfo( gsList, targetID, 'altitude' );
ecefPos         = getGroundStationInfo( gsList, targetID, 'ecefPos' );
ecefVel         = getGroundStationInfo( gsList, targetID, 'ecefVel' );
eciPos          = getGroundStationInfo( gsList, targetID, 'eciPos', epoch );
eciVel          = getGroundStationInfo( gsList, targetID, 'eciVel', epoch );

if( strcmpi(location,targetLocation)==false )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetLocation');
end
if( strcmpi(nasaID,targetNASAID)==false )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetNASAID');
end
if( abs(latitude-targetLatitude) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetLatitude');
end
if( abs(longitude-targetLongitude) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetLongitude');
end
if( abs(altitude-targetAltitude) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetAltitude');
end
if( max(abs(ecefPos-targetECEFPos)) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetECEFPos');
end
if( max(abs(ecefVel-targetECEFVel)) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetECEFVel');
end
if( max(abs(eciPos-targetECIPos)) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetECIPos');
end
if( max(abs(eciVel-targetECIVel)) > tol )
    failed = 1;
    disp('FAILED getGroundStationInfo_test: targetECIVel');
end
