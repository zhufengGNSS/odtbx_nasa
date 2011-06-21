function y = getGroundStationInfo(gsList, stdn, param, time)

% GETGROUNDSTATIONINFO Returns the value of the specified ground station parameter.
%
%   y = getGroundStationInfo(gsList, stdn, param, time) returns the value of the
%   input parameter for the specified ground station. The list of available
%   parameters is below:
%
%   PARAMETER       DESCRIPTION                     TYPE/UNITS
%
%   NASA_ID         NASA Station ID                 string
%   location        city, state, country            string
%   latitude        latitude                        degs
%   longitude       longitude                       degs
%   altitude        height above ellipsoid          km
%   ecefPosition    ECEF position vector            km
%   ecefVelocity    ECEF velocity vector            km/s
%   eciPosition     ECI (J2000) position vector     km
%   eciVelocity     ECI (J2000) velocity vector     km/s
%
%   The list of groundstations is created by using the
%   createGroundStationList function.
%
%   INPUTS
%   VARIABLE        SIZE    		DESCRIPTION
%      gsList       Java object     List of ground stations from
%                                   createGroundStationList
%      stdn         string      	Station ID
%      param        (1x1)     		String specifying one of the parameters above
%      time         (1x1)           UTC in matlab datenum format
%
%   OUTPUTS
%      y            (1x1)			Output from function, Value of the field
%
%
%    keyword: groundstation, JAT Adapter
%    See also CREATEGROUNDSTATIONLIST
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
%               		(MM/DD/YYYY)
%   Derek Surka          08/27/2007     Original

gsIndex = findGS(gsList, stdn);

if( gsIndex < 0 )
    error(['There is no groundstation in the list matching ', stdn]);
end

validNames = {
    'nasa_id'
    'location'
    'latitude'
    'longitude'
    'altitude'
    'ecefposition'
    'ecefvelocity'
    'eciposition'
    'ecivelocity' };

[index,fullName] = getIndex(param,validNames);

switch fullName
    case 'nasa_id'
        y = gsList.get(gsIndex).getNASA();

    case 'location'
        y = gsList.get(gsIndex).getLocation();

    case 'latitude'
        y = gsList.get(gsIndex).getLatitude() * 180/pi;

    case 'longitude'
        y = gsList.get(gsIndex).getLongitude() * 180/pi;

    case 'altitude'
        y = gsList.get(gsIndex).getHAE() / 1000;

    case 'ecefposition'
        r = gsList.get(gsIndex).getECEFPosition();
        %y = [r.x(1); r.x(2); r.x(3)] / 1000;
	y = r / 1000;

    case 'ecefvelocity'
        v = gsList.get(gsIndex).getECEFVelocity();
        %y = [v.x(1); v.x(2); v.x(3)] / 1000;
	y = v / 1000;

    case 'eciposition'
        if( nargin < 4 )
            error('A date must be specified to retrieve ECI position');
        else
            eRef    = jat.spacetime.EarthRef(matlabTime2MJD(time));
            eci2ecf = eRef.ECI2ECEF();
            r = gsList.get(gsIndex).getECIPosition(eci2ecf);
            y = [r.x(1); r.x(2); r.x(3)] / 1000;
        end

    case 'ecivelocity'
        if( nargin < 4 )
            error('A date must be specified to retrieve ECI position');
        else
            eRef    = jat.spacetime.EarthRef(matlabTime2MJD(time));
            eci2ecf = eRef.ECI2ECEF();
            vECI = gsList.get(gsIndex).getECIVelocity(eci2ecf);
            y = [vECI.x(1); vECI.x(2); vECI.x(3)] / 1000;
            
            % Below does the same as what is done in JAT - must only be
            % sure that values of wEarth match
%             w = [0;0;JATConstant('wEarth')];
%             w = [0;0;7.292115E-05];
%             ecf2eci = jatDCM('ece',time);
%             Aecf2eci =rotransf(-w,ecf2eci);
%             rECF = gsList.get(gsIndex).getECEFPosition();
%             vECF = gsList.get(gsIndex).getECEFVelocity();
%             xECEF = [rECF.x(1); rECF.x(2); rECF.x(3); vECF.x(1); vECF.x(2); vECF.x(3)];
%             xECI = Aecf2eci(1:6,1:6)*xECEF;
%             y = xECI(4:6);
        end

end

end


%-----------------------------------------------------------
function y = findGS(gsList, str)

gsSize = gsList.size();

found = 0;
i     = 0;
y     = -1;

while( (i<gsSize) && (~found) )
    if( strcmpi( gsList.get(i).getSTDN(), str ) )
        found = 1;
        y     = i;
    else
        i = i+1;
    end
end
end




