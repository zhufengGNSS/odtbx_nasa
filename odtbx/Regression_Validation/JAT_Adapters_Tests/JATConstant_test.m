function failed = JATConstant_test()
% Regression Test Case
% Function(s) JATConstant
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

namevals  = {
    'arcsec2rad'       4.848136811095360e-006
    'au'               1.495978700000000e+011
    'c'                             299792458
    'days2sec'                          86400
    'deg2rad'          1.745329251994330e-002
    'epsilon'          2.343929111000000e+001
    'fEarth'           3.352810000000000e-003
    'g0'              -9.779999999999999e+000
    'j2Earth'          1.082630000000000e-003
    'L1Frequency'      1.575420000000000e+009
    'L1Wavelength'     1.902936727983649e-001
    'mjdJ2000'         7.304865000000000e+005
    'muEarth'          3.986004415000000e+014
    'muMoon'           4.902794900000000e+012
    'muSun'            1.327122000000000e+020
    'pSolar'           4.560000000000000e-006
    'rad2deg'          5.729577951308232e+001
    'rEarth'           6.378136300000000e+006
    'rSun'                          695990000
    'sec2days'         1.157407407407407e-005
    'wEarth'           7.292115000000000e-005
    };
    
n = length(namevals); 

for i=1:n
    val = JATConstant(namevals{i,1});
    if( abs(val-namevals{i,2}) > tol )
        failed = 1;
    end
end

namevals2  = {
    'earth'             6.37812e6
    'moon'              1.738e6;
    'sun'               6.9599e8;
    'mercury'           6.37812e6 * 0.382;
    'venus'             6.37812e6 * 0.949;
    'mars'              6.37812e6 * 0.532;
    'jupiter'           6.37812e6 * 11.209;
    'saturn'            6.37812e6 * 9.49;
    'uranus'            6.37812e6 * 4.007;
    'neptune'           6.37812e6 * 3.83;
    'pluto'             6.37812e6 * 0.18;
    };
    
n = length(namevals2); 

for i=1:n
    val = JATConstant('meanRadius',namevals2{i,1});
    if( abs(val-namevals2{i,2}) > tol )
        failed = 1;
    end
end