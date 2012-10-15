
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
% Author: Brent Wm. Barbee
%         brent.barbee@emergentspace.com
%         Emergent Space Technologies, Inc.
%
%         May 2009
%

% Begin function.
function [failed] = jatForces_test

    failed = 0;

    load jatForces_reg_data.mat;

    % Initialize options structure.
    jOptions = odtbxOptions('force');

    % Set options values.
    jOptions = setOdtbxOptions(jOptions, 'epoch', JATConstant('MJDJ2000') );
    jOptions = setOdtbxOptions(jOptions, 'cD', 2.2);
    jOptions = setOdtbxOptions(jOptions, 'cR', 0.7);
    jOptions = setOdtbxOptions(jOptions, 'mass', 1000);
    jOptions = setOdtbxOptions(jOptions, 'dragArea', 10, 'srpArea', 8);
    jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM3');
    jOptions = setOdtbxOptions(jOptions, 'gravDegree', 4, 'gravOrder', 4);
    jOptions = setOdtbxOptions(jOptions, 'useSolarGravity', true);
    jOptions = setOdtbxOptions(jOptions, 'useLunarGravity', true);
    jOptions = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', true);
    jOptions = setOdtbxOptions(jOptions, 'useAtmosphericDrag', true);
    jOptions = setOdtbxOptions(jOptions, 'atmosphereModel', 'HP');
    jOptions = setOdtbxOptions(jOptions, 'nParameterForHPModel', 2);
    jOptions = setOdtbxOptions(jOptions, 'f107Daily', 150);
    jOptions = setOdtbxOptions(jOptions, 'f107Average', 150);
    jOptions = setOdtbxOptions(jOptions, 'ap', 15);

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
    
    % Assemble the initial state vector in units of km and km/sec.
    R0_km = [x; y; z; xdot; ydot; zdot];
    
    % Assemble the initial state vector in units of meters and meters/sec.
    R0_m = 1000.0*R0_km;
    
    % Specify simulation time parameters.
    t0 =     0.0; % initial time, seconds
    dt =    60.0; % step size, seconds
    tf = 28980.0; % final time, seconds (~ 1 orbit period)
    
    % Assemble the time span array.
    tspan = (t0:dt:tf);
    
    % jatWorld is an ODToolboxJATModel object that stores relevant force
    % and spacecraft information.
    jatWorld = createJATWorld(jOptions);
    
    % Set integrator options.
    odeopts    = odeset('reltol', 1e-9, 'abstol', 1e-9, 'initialstep', dt); 
    odeOptions = setOdtbxOptions('OdeSolver', @ode45, 'OdeSolvOpts', odeopts);
    
    %
    % jatForces
    %
    
    [t_m, R_m] = integ(@jatForces, tspan, R0_m, odeOptions, jatWorld);
    
    % Transpose the results matrix to match standard convention.
    R_m = R_m';
    
    % Convert back to km for comparisons.
    R_m = (1/1000.0)*R_m;
    
    %
    % jatForces_km
    %
    
    [t_km, R_km] = integ(@jatForces_km, tspan, R0_km, odeOptions, jatWorld);
    
    % Transpose the results matrix to match standard convention.
    R_km = R_km';
    
    
    [Rdot_m, A_m, Q_m] = jatForces(0, R0_m, jatWorld);
    
    [Rdot_km, A_km, Q_km] = jatForces_km(0, R0_km, jatWorld);

    
    % Uncomment the following code to save a new set of regression test
    % data. Re-comment the code afterwards (and be sure to move the saved
    % .mat file to the proper directory).
    %t_m_store     = t_m;
    %R_m_store     = R_m;
    %t_km_store    = t_km;
    %R_km_store    = R_km;
    %Rdot_m_store  = Rdot_m;
    %A_m_store     = A_m;
    %Q_m_store     = Q_m;
    %Rdot_km_store = Rdot_km;
    %A_km_store    = A_km;
    %Q_km_store    = Q_km;
    %
    %save jatForces_reg_data.mat t_m_store R_m_store t_km_store R_km_store Rdot_m_store A_m_store Q_m_store Rdot_km_store A_km_store Q_km_store
    % 
    %end
    
    % deriv check based on the same initial state and time:
    if(max(max(abs(R_km(1,:) - R_m(1,:)))) > 1e-15)
        warning('jatForces_test: km and m results mismatch.');
        failed = 1;
    end
    
    % deriv check over the entire two trajectories
    if(max(max(abs(R_km - R_m))) > 3e-5)
        warning('jatForces_test: km and m results mismatch.');
        failed = 1;
    end
    
    if(max(max(abs(R_m_store - R_m))) > 1e-8)
        warning('jatForces_test: meters state results do not match stored data.');
        failed = 1;
    end
    
    if(max(abs(t_m_store - t_m)) > 1e-15)
        warning('jatForces_test: meters time results do not match stored data.');
        failed = 1;
    end
    
    if(max(max(abs(R_km_store - R_km))) > 2e-9)
        warning('jatForces_test: km state results do not match stored data.');
        failed = 1;
    end
    
    if(max(abs(t_km_store - t_km)) > 1e-15)
        warning('jatForces_test: km time results do not match stored data.');
        failed = 1;
    end
    
    if(max(abs(Rdot_m_store - Rdot_m)) > 1e-15)
        warning('jatForces_test: meters derivs results do not match stored data.');
        failed = 1;
    end
    
    if(max(abs(Rdot_km_store - Rdot_km)) > 1e-15)
        warning('jatForces_test: km derivs results do not match stored data.');
        failed = 1;
    end
    
    if(max(max(abs(A_m_store - A_m))) > 1e-15)
        warning('jatForces_test: meters STM results do not match stored data.');
        failed = 1;
    end
    
    if(max(max(abs(Q_m_store - Q_m))) > 1e-24)
        warning('jatForces_test: meters Q matrix results do not match stored data.');
        failed = 1;
    end
    
    if(max(max(abs(A_km_store - A_km))) > 1e-15)
        warning('jatForces_test: km STM results do not match stored data.');
        failed = 1;
    end
    
    if(max(max(abs(Q_km_store - Q_km))) > 1e-24)
        warning('jatForces_test: km Q matrix results do not match stored data.');
        failed = 1;
    end
    
    if failed
        % plot the errors:
        figure;
        hold on;
        plot(t_m,R_m(:,1:3) - R_m_store(:,1:3));
        plot(t_m,R_km(:,1:3) - R_m_store(:,1:3));
        hold off;
        xlabel('time (sec)');
        ylabel('Err (km/s)');
        legend('x','y','z','x_km','y_km','z_km');
        title('jatForces and jatForces_km Regression Position State Deriv Differences');
        figure;
        hold on;
        plot(t_m,R_m(:,4:6) - R_m_store(:,4:6));
        plot(t_m,R_km(:,4:6) - R_m_store(:,4:6));
        hold off;
        xlabel('time (sec)');
        ylabel('Err (km/s^2)');
        legend('x','y','z','x_km','y_km','z_km');
        title('jatForces() Regression Velocity State Deriv Differences');
        
    end
    
% End of function.
end

























