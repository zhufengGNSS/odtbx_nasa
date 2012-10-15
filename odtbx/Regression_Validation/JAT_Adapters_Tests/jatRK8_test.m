
%
% jatRK8_test
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

%
% Brent Wm. Barbee
% brent.barbee@emergentspace.com
% Emergent Space Technologies, Inc.
%

% Begin function.
function [failed] = jatRK8_test

    % Default to 'pass'.
    failed = 0;
    
    % Tolerance for considering orbit propagation data to be equivalent.
    % In this case the units are meters, meters/second, and seconds since
    % position, velocity, and time are all tested against this same
    % tolerance value.
    JGM2_TOL = 1e-07;
    JGM3_TOL = 1e-12;

    % The initial conditions specified below are for an elliptical inclined
    % Earth-centered orbit with the following orbital parameters:
    %    
    % sma  = 20378.137; % km
    % ecc  =     0.65;
    % incl =    70.0;   % deg
    % raan =    40.0;   % deg
    % argp =   121.0;   % deg
    % tran =   258.0;   % deg
    %
    x    = 8881.95135138262;    % km
    y    = 9430.77519139708;    % km
    z    = 4162.93561538469;    % km
    xdot =   -4.71684695624137; % km/s
    ydot =   -2.37093913899615; % km/s
    zdot =    3.34006991068724; % km/s

    % Specify spacecraft physical parameters and force model options.
    coefD     = 2.22;    % Spacecraft drag coefficient.
    cr        = 1.1;     % Spacecraft coefficient of reflectivity.
    mass      = 500.0;   % Spacecraft mass, kg.
    cArea     = 10.0;    % Spacecraft cross-sectional area, m^2.
    mjd_utc   = 53100.0; % Modified Julian Date (UTC) of initial epoch.
    JGMOrder  = 18;      % Order for JGM gravity model.
    JGMDegree = 18;      % Degree for JGM gravity model.
    
    % Specify simulation time parameters.
    t0 =     0.0; % initial time, seconds
    dt =    60.0; % step size, seconds
    tf = 28980.0; % final time, seconds (~ 1 orbit period)
    
    % Set the time span.
    tspan = (t0:dt:tf);
    
    % Assemble the initial state vector in units of meters and meters/sec.
    R0 = 1000.0*[x y z xdot ydot zdot];
    
    % Assemble data structure of spacecraft parameters force model options.
    options = setJatRK8Options('stepSize', dt, 'coefD', coefD, 'cr', cr, 'mass', mass, 'cArea', cArea, 'mjd_utc', mjd_utc, 'JGMOrder', JGMOrder, 'JGMDegree', JGMDegree);
    
    % Call the numerical integrator with JGM2 EOMs.
    [tjgm2, yjgm2] = jatRK8('JatUniverseJGM2', tspan, R0, options);
    
    % Call the numerical integrator with the JGM3 EOMs.
    [tjgm3, yjgm3] = jatRK8('JatUniverseJGM3', tspan, R0, options);
    
    % Load the regression test data. The variables stored in this .mat file
    % are: tjgm2_store, yjgm2_store, tjgm3_store, and yjgm3_store.
    load jatRK8_reg_dat.mat
    
    % Compute the differences between the current JGM2 integration results
    % and the stored regression test data.
    diff_tjgm2 = tjgm2 - tjgm2_store;
    diff_jgm2  = yjgm2 - yjgm2_store;
    
    % Compute the differences between the current JGM3 integration results
    % and the stored regression test data.
    diff_tjgm3 = tjgm3 - tjgm3_store;
    diff_jgm3 = yjgm3 - yjgm3_store;
    
    % Test the JGM2 results.
    if(max(max(abs(diff_jgm2))) > JGM2_TOL)
    
        % If tolerance requirement is not met, set to 'failed'.
        failed = 1;
        
        fprintf(2, '\n*** Warning: State mismatch for JGM2 results.\n');
        
    end
    
    if(max(abs(diff_tjgm2)) > JGM2_TOL)
    
        % If tolerance requirement is not met, set to 'failed'.
        failed = 1;
        
        fprintf(2, '\n*** Warning:Time value mismatch for JGM2 results.\n');
        
    end

    % Test the JGM3 results.
    if(max(max(abs(diff_jgm3))) > JGM3_TOL)
    
        % If tolerance requirement is not met, set to 'failed'.
        failed = 1;
        
        fprintf(2, '\n*** Warning: State mismatch for JGM3 results.\n');
        
    end
    
    if(max(abs(diff_tjgm3)) > JGM3_TOL)
    
        % If tolerance requirement is not met, set to 'failed'.
        failed = 1;
        
        fprintf(2, '\n*** Warning:Time value mismatch for JGM3 results.\n');
        
    end

% End of function.
end
