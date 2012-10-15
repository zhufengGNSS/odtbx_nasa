function y = JATConstant(constant,source)

% JATCONSTANT  Returns the value from JAT of the input constant from the specified source.
%
%   y = JATConstant(constant, source) returns the value from JAT of the
%   input constant from the optionally specified source. 
%
%   Calling JATCONSTANT without any input arguments displays the list of
%   available constants.
%
%   The list of available constants is below:
%
%   CONSTANT     DESCRIPTION                     UNITS       SOURCE
%
%   arcsec2Rad   Conversion factor from          arcsec/rad
%                arcseconds to radians
%   au           Astronomical Unit               m           IAU76 (default)
%   c            Speed of Light                  m/s         IAU76 (default)
%   days2Sec     Conversion factor for days      sec/day
%                  to seconds
%   deg2Rad      Conversion factor for degrees   rad/deg
%                  to radians                          
%   epsilon      Obliquity of the ecliptic J2000 deg
%   fEarth       Flattening factor of Earth      -           STKHPOP (default)
%                                                            WGS84
%   g0           Gravitational acceleration      m/s^2        
%                near Earth
%   L1Frequency  GPS L1 frequency                Hz
%   L1Wavelength GPS L1 wavelength               m
%   j2Earth      Earth's J2 Value                -           JGM3 (default)
%   mJDJ2000     Modified Julian Date of the     datenum
%                  J2000 epoch in TDB time
%   muEarth      Earth Gravity Constant (GM)     m^3/s^2     JGM3 (default)
%                                                            WGS84
%   muMoon       Moon Gravity Constant (GM)      m^3/s^2     STK (default)
%                                                            DE200
%                                                            JPLSSD
%   muSun        Sun Gravity Constant (GM)       m^3/s^2     STK (default)
%                                                            IAU76
%   muMercury    Mercury Gravitation Constant    m^3/s^2     Will Campbell Code
%                  (GM)
%   muVenus      Venus Gravitation Constant      m^3/s^2     Will Campbell Code
%                  (GM)
%   muMars       Mars Gravitation Constant       m^3/s^2     Will Campbell Code
%                  (GM)
%   muJupiter    Jupiter Gravitation Constant    m^3/s^2     Will Campbell Code
%                  (GM)
%   muSaturn     Saturn Gravitation Constant     m^3/s^2     Will Campbell Code
%                  (GM)
%   muUranus     Uranus Gravitation Constant     m^3/s^2     Will Campbell Code
%                  (GM)
%   muNeptune    Neptune Gravitation Constant    m^3/s^2     Will Campbell Code
%                  (GM)
%   muPluto      Pluto Gravitation Constant      m^3/s^2     Will Campbell Code
%                  (GM)
%   pSolar       Solar pressure at 1AU           N/m^2       IERS96 (default)             
%   rad2Deg      Conversion factor for radians   deg/rad
%                  to degrees
%   rEarth       Equatorial radius of Earth      m           STK (default)
%                                                            WGS84
%   rSun         Mean radius of Sun              m           STK (default)
%                                                            ScienceWorld
%   sec2Days     Conversion factor for seconds   day/sec
%                  to days
%   wEarth       Earth's rotation rate           rad/s       IERS96 (default)
%                                                            WGS84
%                                                            Vallado
%   meanRadius   Mean Radius of planetary bodies  m          Orbit Mechanics - Prussing, Conway 1993
%                  Enter 'EARTH', 'MOON', 'SUN'
%                  or planets for source - names
%                  match output of ephemDE405
%   
%   INPUTS 
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%     constant  (1X1)       String specifying one of the constants above.
%
%     source    (1X1)       Optional string that specifies the source of the data.
%
%   OUTPUTS 
%      y        (1X1)       The value of the constant from JAT
%
%   keyword: JAT Constant 
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
%   Derek Surka          06/04/2007     Original
%   Keith Speckman       05/28/2008     Added hard coded meanRadius option
%   Keith Speckman       10/10/2008     Added grav consts for all planets
%   Kevin Berry          06/22/2009     Added time scale information to the
%                                       mJDJ2000 description in the help.

if( nargin == 0)
    validNames = getJATConstantsNames();
    [~,i] = sort(lower(validNames));
    disp('Available constants are:');
    disp(validNames(i));
    return
elseif( nargin < 2 )
    source = 'default';
end

validNames          = getJATConstantsNames();
errorMessage        = ['JATConstant: There is no JAT constant corresponding to the name ' constant];
[~,fullName]    = getIndex(constant,validNames,errorMessage);

source = lower(source);

switch fullName
    case 'muEarth'
        switch source(1:3)
            case 'wgs'
                y = jat.cm.Constants.GM_WGS84;
            case {'jgm','def'}
                y = jat.cm.Constants.GM_Earth;
            otherwise
                error('Source is not recognized');
        end

    case 'rEarth'
         switch source(1:3)
            case 'wgs'
                y = jat.cm.Constants.rEarth_WGS84;
            case {'stk','def','jgm'}
                y = jat.cm.Constants.rEarth_STKJGM3;
            otherwise
                error('Source is not recognized');
        end
       
    case 'rSun'
        switch source(1:3)
            case 'sci'
                y = jat.cm.Constants.rSun_ScienceWorld;
            case {'stk','def'}
                y = jat.cm.Constants.rSun_STK;
            otherwise
                error('Source is not recognized');
        end

    case 'fEarth'
        switch source(1:3)
            case 'wgs'
                y = jat.cm.Constants.fEarth_WGS84;
            case {'stk','def'}
                y = jat.cm.Constants.fEarth_STKHPOP2;
            otherwise
                error('Source is not recognized');
        end

    case 'wEarth'
        switch source(1:3)
            case 'wgs'
                y = jat.cm.Constants.WE_WGS84;
            case 'val'
                y = jat.cm.Constants.wEarth_Vallado;
            case {'ier','def'}
                y = jat.cm.Constants.wEarth_IERS96;
            otherwise
                error('Source is not recognized');
        end

    case 'j2Earth'
        y = jat.cm.Constants.j2;
        
    case 'c'
        y = jat.cm.Constants.clight;
        
    case 'g0'
        y = jat.cm.Constants.g0;
        
    case 'deg2rad'
        y = jat.cm.Constants.deg2rad;
        
    case 'rad2deg'
        y = jat.cm.Constants.rad2deg;
        
    case 'arcsec2rad'
        y = jat.cm.Constants.arcsec2rad;
        
    case 'pSolar'
        y = jat.cm.Constants.P_Sol;
        
    case 'muMoon'
        switch source(1:3)
            case 'de2'
                y = jat.cm.Constants.muMoon_DE200;
            case 'jpl'
                y = jat.cm.Constants.muMoon_JPLSSD;
            case {'stk','def'}
                y = jat.cm.Constants.muMoon_STK;
            otherwise
                error('Source is not recognized');
        end

    case 'muSun'
        switch source(1:3)
            case 'iau'
                y = jat.cm.Constants.muSun_IAU76;
            case {'stk','def'}
                y = jat.cm.Constants.muSun_STK;
            otherwise
                error('Source is not recognized');
        end

    case 'muMercury'
        y = 22032e9;

    case 'muVenus'
        y = 324859e9;

    case 'muMars'
        y = 42828e9;

    case 'muJupiter'
        y = 126712768e9;

    case 'muSaturn'
        y = 37940626e9;

    case 'muUranus'
        y = 5794549e9;

    case 'muNeptune'
        y = 6836534e9;

    case 'muPluto'
        y = 982e9;

    case 'au'
        y = jat.cm.Constants.AU;
        
    case 'epsilon'
        y = jat.cm.Constants.eps;

    case 'days2sec'
        y = jat.spacetime.TimeUtils.days2sec;
    
    case 'sec2days'
        y = jat.spacetime.TimeUtils.sec2days;
            
    case 'mjdJ2000'
        mjd = jat.spacetime.TimeUtils.MJD_J2000;
        y   = mJD2MatlabTime(mjd);
        
    case 'L1Wavelength'
        y = jat.gps.GPS_Utils.lambda;
        
    case 'L1Frequency'
        y = jat.gps.GPS_Utils.L1_freq;

    case 'meanRadius'
        switch source
            case 'earth'
                y = 6.37812e6;
            case 'moon'
                y = 1.738e6;
            case 'sun'
                y = 6.9599e8;
            case 'mercury'
                y = 6.37812e6 * 0.382;
            case 'venus'
                y = 6.37812e6 * 0.949;
            case 'mars'
                y = 6.37812e6 * 0.532;
            case 'jupiter'
                y = 6.37812e6 * 11.209;
            case 'saturn'
                y = 6.37812e6 * 9.49;
            case 'uranus'
                y = 6.37812e6 * 4.007;
            case 'neptune'
                y = 6.37812e6 * 3.83;
            case 'pluto'
                y = 6.37812e6 * 0.18;
            otherwise
                error('Body not available');
        end

    otherwise % This should never be executed
        error('There is no JAT constant corresponding to that name')
end
        
end

% List of all valid constants
%----------------------------
function validNames = getJATConstantsNames()

validNames = {
    'arcsec2rad'
    'au'
    'c'
    'days2sec'
    'deg2rad'
    'epsilon'
    'fEarth'
    'g0'
    'j2Earth'
    'mjdJ2000'
    'muEarth'
    'muMoon'
    'muSun'
    'muMercury'
    'muVenus'
    'muMars'
    'muJupiter'
    'muSaturn'
    'muUranus'
    'muNeptune'
    'muPluto'
    'pSolar'
    'rad2deg'
    'rEarth'
    'rSun'
    'sec2days'
    'wEarth'
    'L1Wavelength' 
    'L1Frequency'
    'meanRadius' };

end
