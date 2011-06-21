function editratio_test
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

    % Specify the initial reference state.
    %
    % sma  = 12345.6  km
    % ecc  =     0.46
    % incl =     0.00 rad
    % raan =     0.00 rad
    % argp =     0.00 rad
    % tran =     0.00 rad
    %
    xo(1,1) = 6666.624;            % x,  km
    xo(2,1) =    0.000;            % y,  km
    xo(3,1) =    0.000;            % vx, km/s
    xo(4,1) =    9.34312877843682; % vy, km/s

    % Specify the a priori covariance matrix.
    Po = diag([1e-2 1e-4 1e-4 1e-7].^2); % km^2 & km^2/s^2

    % Specify the dynamics models for the truth and the estimator. 
    dynfun.tru = @pr2bp; 
    dynfun.est = @pr2bp;

    % Specify the inputs to the dynamics models
    mu = 3.986004415e5; % km^3/sec^2
    dynarg.tru = mu;
    dynarg.est = mu;

    % Set up the solve-for and consider mapping.
    S = eye(4); % Solve-for map - solve for all 4 states
    C = [];     % Consider map - no consider states
    
    % Use the same a priori for truth and estimator:
    Xnot.Xo    = xo;   % True initial state
    Xnot.Xbaro = S*xo; % Estimator initial state
    Pnot.Po    = Po;   % True initial covariance
    
    if isempty(Po)
        Pnot.Pbaro = [];
    else
        Pnot.Pbaro = S*Po*S'; % Estimator initial covariance
    end
    
    % Specify the measurement models for the truth and the estimator. 
    datfun.tru = @rrdot2D;
    datfun.est = @rrdot2D;
    
    % Specify the inputs to the measurement models, which in this case is the 
    % measurement noise.  Here, the estimator will assume 3x more measurement 
    % noise than the truth:
    sig = diag([1e-3,1e-6]); % km & km/s
    datarg.tru = sig;
    datarg.est = 3*sig;

    % Specify the time vector at which the measurements are taken. 
    tspan = (0:60:(4*3600));

    % Set some options for the sequential filter. 
    opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 1);

    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);

    % Set measurement editing options.
    opts = setOdtbxOptions(opts, 'EditRatio', [0.5 1], 'EditFlag', [1 1]);

    % Run (and time) the sequential estimator.
    tic
    [t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw, Sig_sa, eflag] = estseq(...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    toc

    colstr{1} = 'ro';
    colstr{2} = 'go';
    colstr{3} = 'bo';
    colstr{4} = 'mo';
    colstr{5} = 'ko';
    colstr{6} = 'co';

    k = 1;
    figure(1);
    m = length(eflag);
    for(a=1:m)

        n = size(eflag{a},1);

        for(b=1:n)

            plot(eflag{a}(b,:), colstr{k});

            k = k + 1;

            if(k > 6)
                fprintf(2, '\n *** Number of available color settings exceeded! ***\n');
                return;
            end

            if(a == 1)
                grid on;
            end

            if((a == 1) && ((n > 1) || (m > 1)))
                hold on;
            end

        end

    end

    % Run (and time) the batch estimator.
    %tic
    %[t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw] = estbat(...
    %    dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    %toc
    

    % Plot estimator error, measurment residuals, and variance sandpiles in
    % terms of meters
    %plot_results(t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw, datarg, S, 1e3);


    %========================================%

    %tstep = 30.0; % seconds

    %sig = [];

    %GM = 3.986004415e5; % km^3/s^2

    %KOE.sma  = 12345.6;  % km
    %KOE.ecc  =     0.46;
    %KOE.incl =     0.00; % rad
    %KOE.raan =     0.00; % rad
    %KOE.argp =     0.00; % rad
    %KOE.tran =     0.00; % rad

    % Compute orbit period in units of seconds.
    %P = 2.0*pi*sqrt((KOE.sma^3)/GM);

    % Construct the array of time values for propagation.
    %t = (0.0:tstep:2*P);

    % Propagate the orbit.
    %KOEf = kepprop2b(KOE, t, GM);

    % Convert the orbit time history to cartesian coordinates.
    %X = kep2cart(KOEf, GM);

    %figure;
    %plot(X(1,:), X(2,:), 'b');
    %grid on;
    %xlabel('X, km');
    %ylabel('Y, km');
    %axis equal;

    %figure;
    %plot3(X(1,:), X(2,:), X(3,:), 'b');
    %grid on;
    %xlabel('X, km');
    %ylabel('Y, km');
    %zlabel('Z, km');
    %axis equal;

    %[y, H, R] = range2D(t, X, sig);

end
