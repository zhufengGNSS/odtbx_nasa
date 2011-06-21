%
% kep2cart_test
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

% Brent Wm. Barbee

% Begin function.
function [failed] = kep2cart_test

    % Default to 'pass'.
    failed = 0;

    % Set accuracy tolerance for pass/fail.
    TOL = 2e-5; % = 2 cm (for position) or 2 cm/s (for velocity)
    
    % Gravitational parameter.
    GM = 3.986004415e5; % km^3/s^2
    
    % Test case values to cover the 4 possible orbit types:
    %
    % Case 1: Circular Equatorial
    KOE(1).sma  = 14378.137; % km
    KOE(1).ecc  = 0.0;
    KOE(1).incl = 0.0*(pi/180.0); % rad
    KOE(1).raan = 0.0*(pi/180.0); % rad
    KOE(1).argp = 0.0*(pi/180.0); % rad
    KOE(1).tran = 25.0*(pi/180.0); % rad
    
    % Case 2: Circular Inclined
    KOE(2).sma  = 6778.137; % km
    KOE(2).ecc  = 0.0;
    KOE(2).incl = 62.0*(pi/180.0); % rad
    KOE(2).raan = 360.0*(pi/180.0); % rad
    KOE(2).argp = 0.0*(pi/180.0); % rad
    KOE(2).tran = 127.0*(pi/180.0); % rad
    
    % Case 3: Elliptical Inclined
    KOE(3).sma  = 20378.137; % km
    KOE(3).ecc  = 0.65;
    KOE(3).incl = 70.0*(pi/180.0); % rad
    KOE(3).raan = 40.0*(pi/180.0); % rad
    KOE(3).argp = 121.0*(pi/180.0); % rad
    KOE(3).tran = 258.0*(pi/180.0); % rad
    
    % Case 4: Elliptical Equatorial
    KOE(4).sma  = 11978.137; % km
    KOE(4).ecc  = 0.43;
    KOE(4).incl = 0.0*(pi/180.0); % rad
    KOE(4).raan = 0.0*(pi/180.0); % rad
    KOE(4).argp = 360.0*(pi/180.0); % rad
    KOE(4).tran = 181.0*(pi/180.0); % rad

   
    % Call the kep2cart function.
    X = kep2cart(KOE, GM);
    
    
    % Specify what the results should be (obtained from independent
    % validated code - BWB). Units are km and km/s.
    Xtrue(1,:) = [13031.0175143697 -4079.18467803767 8881.95135138262 -17124.1596122153];
    Xtrue(2,:) = [6076.46326050249 2541.37205609786 9430.77519139708 -298.903317872381];
    Xtrue(3,:) = [0.0 4779.62568419486 4162.93561538469 0.0];
    Xtrue(4,:) = [-2.2251835382434 -6.12438286743884 -4.71684695624137 0.11151264614313];
    Xtrue(5,:) = [4.77192150243214 -2.16663638202305 -2.37093913899615 -3.64105801763251];
    Xtrue(6,:) = [0.0 -4.07485038445297 3.34006991068724 0.0];
    
    % Difference the values returned by kep2cart with what the values
    % should be.
    diff_X = X - Xtrue;
    disp(diff_X);
    
    if(max(max(abs(diff_X))) > TOL)
    
        % If tolerance requirement is not met, set to 'failed'.
        failed = 1;
        
    end

% End of function.
end






























