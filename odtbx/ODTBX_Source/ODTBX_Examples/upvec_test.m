function upvec_test
    % estseq example.
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
    opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 3);

    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 691);

    % Set measurement editing options.
    opts = setOdtbxOptions(opts, 'EditRatio', [9 9], 'EditFlag', [1 1]);

    % Set the update vectorized flag to No.
    opts = setOdtbxOptions(opts, 'UpdateVectorized', 0);

    % Run (and time) the sequential estimator with update vectorized set to No.
    tic
    [t0, xhat0, P0, e0, dy0, Pa0, Pv0, Pw0, Phata0, Phatv0, Phatw0, Sig_sa0, eflag0, Pdy0, Pdyt0] = estseq(...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    toc

    % Set the update vectorized flag to Yes.
    opts = setOdtbxOptions(opts, 'UpdateVectorized', 1);

    % Run (and time) the sequential estimator with update vectorized set to Yes.
    tic
    [t1, xhat1, P1, e1, dy1, Pa1, Pv1, Pw1, Phata1, Phatv1, Phatw1, Sig_sa1, eflag1, Pdy1, Pdyt1] = estseq(...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    toc


    % Compute differences.

    dt = t1{1} - t0{1};

    if(any(dt))
        error('\n***Time mismatch!\n');
    end
    
    dxhat = xhat0{1} - xhat1{1};

    dP = P0{1} - P1{1};

    de = e0{1} - e1{1};

    ddy = dy0{1} - dy1{1};

    def = eflag0{1} - eflag1{1};

    dPa_max = max(max(max(abs(Pa0 - Pa1))));
    dPa_min = min(min(min(abs(Pa0 - Pa1))));

    dPv_max = max(max(max(abs(Pv0 - Pv1))));
    dPv_min = min(min(min(abs(Pv0 - Pv1))));

    dPw_max = max(max(max(abs(Pw0 - Pw1))));
    dPw_min = min(min(min(abs(Pw0 - Pw1))));

    dPhata_max = max(max(max(abs(Phata0 - Phata1))));
    dPhata_min = min(min(min(abs(Phata0 - Phata1))));

    dPhatv_max = max(max(max(abs(Phatv0 - Phatv1))));
    dPhatv_min = min(min(min(abs(Phatv0 - Phatv1))));

    dPhatw_max = max(max(max(abs(Phatw0 - Phatw1))));
    dPhatw_min = min(min(min(abs(Phatw0 - Phatw1))));

    dSigsa_max = max(max(max(abs(Sig_sa0 - Sig_sa1))));
    dSigsa_min = min(min(min(abs(Sig_sa0 - Sig_sa1))));

    dPdy_max = max(max(max(abs(Pdy0{1} - Pdy1{1}))));
    dPdy_min = min(min(min(abs(Pdy0{1} - Pdy1{1}))));

    dPdyt_max = max(max(max(abs(Pdyt0 - Pdyt1))));
    dPdyt_min = min(min(min(abs(Pdyt0 - Pdyt1))));

    fprintf(1, '\n');
    fprintf(1, 'Min dPa --> %e\n', dPa_min);
    fprintf(1, 'Max dPa --> %e\n', dPa_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPv --> %e\n', dPv_min);
    fprintf(1, 'Max dPv --> %e\n', dPv_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPw --> %e\n', dPw_min);
    fprintf(1, 'Max dPw --> %e\n', dPw_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPhata --> %e\n', dPhata_min);
    fprintf(1, 'Max dPhata --> %e\n', dPhata_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPhatv --> %e\n', dPhatv_min);
    fprintf(1, 'Max dPhatv --> %e\n', dPhatv_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPhatw --> %e\n', dPhatw_min);
    fprintf(1, 'Max dPhatw --> %e\n', dPhatw_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dSigsa --> %e\n', dSigsa_min);
    fprintf(1, 'Max dSigsa --> %e\n', dSigsa_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPdy --> %e\n', dPdy_min);
    fprintf(1, 'Max dPdy --> %e\n', dPdy_max);
    fprintf(1, '\n');
    fprintf(1, 'Min dPdyt --> %e\n', dPdyt_min);
    fprintf(1, 'Max dPdty --> %e\n', dPdyt_max);
    fprintf(1, '\n');


    figure;
    subplot(4,1,1);
    plot(t0{1}, dxhat(1,:), 'r');
    ylabel('dxhat_1');
    grid on;
    subplot(4,1,2);
    plot(t0{1}, dxhat(2,:), 'g');
    ylabel('dxhat_2');
    grid on;
    subplot(4,1,3);
    plot(t0{1}, dxhat(3,:), 'b');
    ylabel('dxhat_3');
    grid on;
    subplot(4,1,4);
    plot(t0{1}, dxhat(4,:), 'm');
    ylabel('dxhat_4');
    grid on;
    xlabel('Time');

    figure;
    subplot(5,1,1);
    plot(t0{1}, dP(1,:), 'r');
    ylabel('dP_1');
    grid on;
    subplot(5,1,2);
    plot(t0{1}, dP(2,:), 'g');
    ylabel('dP_2');
    grid on;
    subplot(5,1,3);
    plot(t0{1}, dP(3,:), 'b');
    ylabel('dP_3');
    grid on;
    subplot(5,1,4);
    plot(t0{1}, dP(4,:), 'm');
    ylabel('dP_4');
    grid on;
    subplot(5,1,5);
    plot(t0{1}, dP(5,:), 'm');
    ylabel('dP_5');
    grid on;
    xlabel('Time');

    figure;
    subplot(5,1,1);
    plot(t0{1}, dP(6,:), 'r');
    ylabel('dP_6');
    grid on;
    subplot(5,1,2);
    plot(t0{1}, dP(7,:), 'g');
    ylabel('dP_7');
    grid on;
    subplot(5,1,3);
    plot(t0{1}, dP(8,:), 'b');
    ylabel('dP_8');
    grid on;
    subplot(5,1,4);
    plot(t0{1}, dP(9,:), 'm');
    ylabel('dP_9');
    grid on;
    subplot(5,1,5);
    plot(t0{1}, dP(10,:), 'm');
    ylabel('dP_{10}');
    grid on;
    xlabel('Time');

    figure;
    subplot(4,1,1);
    plot(t0{1}, de(1,:), 'r');
    ylabel('de_1');
    grid on;
    subplot(4,1,2);
    plot(t0{1}, de(2,:), 'g');
    ylabel('de_2');
    grid on;
    subplot(4,1,3);
    plot(t0{1}, de(3,:), 'b');
    ylabel('de_3');
    grid on;
    subplot(4,1,4);
    plot(t0{1}, de(4,:), 'm');
    ylabel('de_4');
    grid on;
    xlabel('Time');

    figure;
    subplot(2,1,1);
    plot(t0{1}, ddy(1,:), 'r*');
    ylabel('ddy_1');
    grid on;
    subplot(2,1,2);
    plot(t0{1}, ddy(2,:), 'g*');
    ylabel('ddy_2');
    grid on;
    xlabel('Time');

    figure;
    subplot(2,1,1);
    plot(t0{1}, def(1,:), 'r*');
    ylabel('def_1');
    grid on;
    subplot(2,1,2);
    plot(t0{1}, def(2,:), 'g*');
    ylabel('def_2');
    grid on;
    xlabel('Time');

    %colstr{1} = 'ro';
    %colstr{2} = 'go';
    %colstr{3} = 'bo';
    %colstr{4} = 'mo';
    %colstr{5} = 'ko';
    %colstr{6} = 'co';

    %k = 1;
    %figure(1);
    %m = length(eflag);
    %for(a=1:m)
    %
    %    n = size(eflag{a},1);
    %
    %    for(b=1:n)
    %
    %        plot(eflag{a}(b,:), colstr{k});
    %
    %        k = k + 1;
    %
    %        if(k > 6)
    %            fprintf(2, '\n *** Number of available color settings exceeded! ***\n');
    %            return;
    %        end
    %
    %        if(a == 1)
    %            grid on;
    %        end
    %
    %        if((a == 1) && ((n > 1) || (m > 1)))
    %            hold on;
    %        end
    %
    %    end
    %
    %end

    % Run (and time) the batch estimator.
    %tic
    %[t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw] = estbat(...
    %    dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    %toc
    

    % Plot estimator error, measurment residuals, and variance sandpiles in
    % terms of meters
    %plot_results(t, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw, Pdy, Pdyt, S, 1e3);

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
