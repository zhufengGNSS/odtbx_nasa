function d = jatDCM(type, time, todExpireDays)

% JATDCM  Returns the requested Direction Cosine Matrix from JAT.
%
%   d = jatDCM(type, time) returns a 3x3 direction cosine matrix, or a 3x3xN
%   array of DCMs, of the TYPE the user specifies at the requested times.  
%   Times are input in Matlab datenum format.
%
%   INPUTS
%   VARIABLE  SIZE      DESCRIPTION (Optional/Default)
%      type   (1x1)     String indicating type of DCM to return:
%         'eci2ecef'      ECI to ECEF transformation
%         'ecef2eci'      ECEF to ECI transformation
%         'tod'           J2000 to True of Date transformation
%         'precession'    Precession transformation of equatorial coordinates
%         'nutation'      Transformation from mean to true equator and equinox (Nutation matrix)
%         'gha'           Transformation from true equator and equinox to Earth equator and Greenwich meridian system
%         'polar'         Transformation from pseudo Earth-fixed to Earth-fixed coordinates (Polar Motion matrix)
%         'ecliptic'      Transformation of equatorial to ecliptical coordinates
%       time   (1xn)    UTC Times to compute DCMs
%      todExpireDays 1x1 In days.  For efficiency, eci-to-ecef
%                       calculations can reuse true of date (precession and
%                       nutation) calculations from previous time samples.
%                       This specifies for how long precession and nutation
%                       can be reused.  Note that reuse of precession and 
%                       nutation is within one call to jatDCM.  Subsequent
%                       calls to jatDCM will not reuse precession and
%                       nutation from previous calls no matter what
%                       expiration is set.  0 means always recalculate.
%                       This parameter is ignored for types besides
%                       eci2ecef and ecef2eci.
%
%   OUTPUTS 
%      d       (3x3xn)  Direction cosine matrix/matrices
%
%   Notes:
%     1) PolarMatrix always equals identity. This can/should be modified in
%        future versions of JAT to enable user to enter polar coefficients.
%
%     2) TODMatrix = NutationMatrix * PrecessionMatrix
%
%     3) ECI2ECEFMatrix = PolarMatrix * GHAMatrix * TODMatrix
%
%   keyword: Coordinate Transformations, JAT Adapter
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               		(MM/DD/YYYY)
%   Derek Surka          08/16/2007     Original
%   Kevin Berry          06/22/2009     Added time scale information to the
%                                       description in the help.

% Allowable types
validNames = {
    'eci2ecef'
    'ecef2eci'
    'tod'
    'precession'
    'nutation'
    'gha'
    'polar'
    'ecliptic'
    };

if( nargin == 0 )
    [temp,i] = sort(lower(validNames));
    disp('Available DCMs are:');
    disp(validNames(i));
    return
end

if( nargin < 2 )
    time = 0;
%     time = now;
end

if (nargin < 3)
    todExpireDays = 0;
end

[index,fullName] = getIndex(type,validNames);

% initialize output d
n = length(time);
mjd_utc = matlabTime2MJD(time);
switch fullName
    case 'eci2ecef'
        d = jat.matlabInterface.BatchCalcs.eciToEcefXform(mjd_utc, todExpireDays);
    case 'ecef2eci'
        d = jat.matlabInterface.BatchCalcs.ecefToEciXform(mjd_utc, todExpireDays);
    case 'tod'
        d = jat.matlabInterface.BatchCalcs.j2000ToTODXform(mjd_utc);
    case 'precession'
        d = jat.matlabInterface.BatchCalcs.getPrecession(mjd_utc);
    case 'nutation'
        d = jat.matlabInterface.BatchCalcs.getNutation(mjd_utc);
    case 'gha'
        d = jat.matlabInterface.BatchCalcs.getGHAMatrix(mjd_utc);
    case 'polar'
        d = jat.matlabInterface.BatchCalcs.getPoleMatrix(mjd_utc);
    case 'ecliptic'
        d = jat.matlabInterface.BatchCalcs.getEcliptic(mjd_utc);
end

    
end
