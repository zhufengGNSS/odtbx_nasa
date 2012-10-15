function failed = ephemDE405_test()
% Regression Test Case
% Function(s) ephemDE405
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
tol     = 1e-6;
disp('Tolerance has been dropped from 1e-10 to 1e-6 km')
epoch   = datenum('June 1, 2007');

load ephemDE405_testdata

bodies  = {
    'EARTH'
    'MERCURY'
    'VENUS'
    'EM_BARY'
    'MARS'
    'JUPITER'
    'SATURN'
    'URANUS'
    'NEPTUNE'
    'PLUTO'
    'GEOCENTRIC_MOON'
    'SUN'
    'MOON_LIB'
    'GEOCENTRIC_SUN'
    'SOLAR_SYSTEM_BARY'
    'MOON'
    };
    
n       = length(bodies); 

for i=1:n
    [p,v] = ephemDE405(bodies{i}, epoch);
    if( any( abs(p-targetPos(:,i)) > tol ) || any( abs(v-targetVel(:,i)) > tol ) )
        failed = 1;
	a=norm( abs(p-targetPos(:,i)))
    end
end

epoch_UTC = convertTime('UTC', 'TT', epoch);
for i=1:n
    [p,v] = ephemDE405(bodies{i}, epoch_UTC, 'UTC');
    if( any( abs(p-targetPos(:,i)) > tol ) || any( abs(v-targetVel(:,i)) > tol ) )
        failed = 1;
	a=norm( abs(p-targetPos(:,i)))
    end
end