function [pos, vel] = ephemDE405(planet, time, timescale)

% EPHEMDE405  Computes the position and velocity of the input planet at the given input time.
% The outputs are in the (J2000) International Celestial Reference Frame (ICRF).
%
%   [pos, vel] = ephemDE405(planet, time) calls JAT to compute the ICRF 
% position and velocity of the input planet at the given input Terrestrial
% Time (TT) in datenum format. 
%
%   [pos, vel] = ephemDE405(planet, time, timescale) calls JAT to compute
% the ICRF position and velocity of the input planet at the given input
% time in datenum format in the given time scale. All time scales are
% converted to TT before being sent to Jat.
%   
%   The latest JPL documentation about the DE405 ephemeris can be found at
%       ftp://ssd.jpl.nasa.gov/pub/eph/export/DE405/
%
%   INPUTS 
%   VARIABLE      SIZE    	DESCRIPTION (Optional/Default)
%      planet    (1X1)      Name of desired ephemeris. A list of valid
%                             names can be found by calling ephemDE405 
%                             without any input arguments.
%      time      (1X1)      Time in Matlab datenum time format.
%                             Default is time=0 (JAT gives an error)
%      timescale (1X1)      (Optional)String specifying input timescale. 
%                             A list of valid time scales can be found in
%                             convertTime.m
%                             Default = 'TT'
%   OUTPUTS 
%      pos       (3X1)      Position in km (ICRF)
%      vel       (3X1)      Velocity in km/s (ICRF)
%
%   keyword: JAT Adapter, DE405 
%   See also  
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
%               	   (MM/DD/YYYY)
%   Derek Surka         05/18/2007   	Original
%   Derek Surka         08/06/2007      Modified to work with new JAT DE405
%   Keith Speckman      03/22/2008      Modifed % DE405 output is of type VectorN section to
%                                         properly extract the java structure
%   Keith Speckman      08/04/2008      New JAT outputs in m - converted back to km
%   Kevin Berry         06/11/2009      Added timescale optional input

if( nargin == 0 )
    displayDE405Names;
    return
end

if( nargin < 2 )
    time = 0;
%     time = now;
end

if( nargin >= 3 )
    time = convertTime('TT', timescale, time);
end

de405Body = convertPlanetToBody( planet ); 
mjd = matlabTime2MJD(time);

% DMS must use following format because get_planet_posvel is not static
y = get_planet_posvel(jat.eph.DE405, de405Body, mjd);

% DE405 output is of type VectorN
outstruct = struct(y);
pos = outstruct.x(1:3) * 1e-3;
vel = outstruct.x(4:6) * 1e-3;

end

%-------------------------------------------------------------------------
function body = convertPlanetToBody( planet )
    validNames       = getDE405Names;
    [index,fullName] = getIndex(planet,validNames);
    bodyStr = ['body = jat.eph.DE405_Body.' fullName ';'];
    eval(bodyStr);
end

%-------------------------------------------------------------------------
function displayDE405Names()
    names = getDE405Names();
    disp('Valid DE405 body names are:');
    disp(names)
end

%-------------------------------------------------------------------------
function [names, bodies] = getDE405Names()
    bodies = jat.eph.DE405_Body.values();
    for i=1:length(bodies)
        names(i) = bodies(i).name();
    end
    names = cellstr(char(names));
end
