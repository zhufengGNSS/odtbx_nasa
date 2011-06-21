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
function [failed] = montecarloseed_test

    % Default to pass.
    failed = 0;

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

    % Specify the time vector at which the measurements are taken for the sequential estimator cases. 
    tspan = (0:60:(4*3600));

    % Specify the time vector at which the measurements are taken for the batch estimator cases. 
    tspan_bat = (0:120:(3*3600));

    %*********************************************************************************************%
    % Block of code for generating and storing test data sets.
    %
    %==========%
    %= CASE 1 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    % opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    % opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    %
    % Run the sequential estimator.
    % [t_c1, xhat_c1, P_c1, e_c1, dy_c1, Pa_c1, Pv_c1, Pw_c1, Phata_c1, Phatv_c1, Phatw_c1, ...
    % Sig_sa_c1, eflag_c1, Pdy_c1, Pdyt_c1] = estseq(dynfun, datfun, tspan, Xnot, Pnot, ...
    %                                                opts, dynarg, datarg, S, C);
    %
    % Store the test case data.
    % save montecarloseed_test_data_case1.mat t_c1 xhat_c1 P_c1 e_c1 dy_c1 Pa_c1 Pv_c1 Pw_c1 ...
    %                                        Phata_c1 Phatv_c1 Phatw_c1 Sig_sa_c1 eflag_c1 ...
    %                                        Pdy_c1 Pdyt_c1
    %
    %==========%
    %= CASE 2 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    % opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    % opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    %
    % Run the sequential estimator.
    % [t_c2, xhat_c2, P_c2, e_c2, dy_c2, Pa_c2, Pv_c2, Pw_c2, Phata_c2, Phatv_c2, Phatw_c2, ...
    % Sig_sa_c2, eflag_c2, Pdy_c2, Pdyt_c2] = estseq(dynfun, datfun, tspan, Xnot, Pnot, ...
    %                                                opts, dynarg, datarg, S, C);
    %
    % Store the test case data.
    % save montecarloseed_test_data_case2.mat t_c2 xhat_c2 P_c2 e_c2 dy_c2 Pa_c2 Pv_c2 Pw_c2 ...
    %                                        Phata_c2 Phatv_c2 Phatw_c2 Sig_sa_c2 eflag_c2 ...
    %                                        Pdy_c2 Pdyt_c2
    %
    %==========%
    %= CASE 3 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    % opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    % opts = setOdtbxOptions(opts, 'MonteCarloSeed', [10 20 30]);
    %
    % Run the sequential estimator.
    % [t_c3, xhat_c3, P_c3, e_c3, dy_c3, Pa_c3, Pv_c3, Pw_c3, Phata_c3, Phatv_c3, Phatw_c3, ...
    % Sig_sa_c3, eflag_c3, Pdy_c3, Pdyt_c3] = estseq(dynfun, datfun, tspan, Xnot, Pnot, ...
    %                                                opts, dynarg, datarg, S, C);
    %
    % Store the test case data.
    % save montecarloseed_test_data_case3.mat t_c3 xhat_c3 P_c3 e_c3 dy_c3 Pa_c3 Pv_c3 Pw_c3 ...
    %                                        Phata_c3 Phatv_c3 Phatw_c3 Sig_sa_c3 eflag_c3 ...
    %                                        Pdy_c3 Pdyt_c3
    %
    %==========%
    %= CASE 4 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    %opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    %opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    %
    % Run the batch estimator.
    %[t_c4, xhat_c4, P_c4, e_c4, dy_c4, Pa_c4, Pv_c4, Pw_c4, Phata_c4, Phatv_c4, Phatw_c4, ...
    % Sig_sa_c4, Pdy_c4, Pdyt_c4] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
    %                                      dynarg, datarg, S, C);
    %
    % Store the test case data.
    %save montecarloseed_test_data_case4.mat t_c4 xhat_c4 P_c4 e_c4 dy_c4 Pa_c4 Pv_c4 Pw_c4 ...
    %                                        Phata_c4 Phatv_c4 Phatw_c4 Sig_sa_c4 Pdy_c4 ...
    %                                        Pdyt_c4
    %
    %==========%
    %= CASE 5 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    %opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    %opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    %
    % Run the batch estimator.
    %[t_c5, xhat_c5, P_c5, e_c5, dy_c5, Pa_c5, Pv_c5, Pw_c5, Phata_c5, Phatv_c5, Phatw_c5, ...
    % Sig_sa_c5, Pdy_c5, Pdyt_c5] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
    %                                      dynarg, datarg, S, C);
    %
    % Store the test case data.
    %save montecarloseed_test_data_case5.mat t_c5 xhat_c5 P_c5 e_c5 dy_c5 Pa_c5 Pv_c5 Pw_c5 ...
    %                                        Phata_c5 Phatv_c5 Phatw_c5 Sig_sa_c5 Pdy_c5 ...
    %                                        Pdyt_c5
    %
    %==========%
    %= CASE 6 =%
    %==========%
    % Set the number of Monte Carlo cases and the number of update iterations. 
    %opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    %
    % Set the Monte Carlo Seed.
    %opts = setOdtbxOptions(opts, 'MonteCarloSeed', [10 20 30]);
    %
    % Run the batch estimator.
    %[t_c6, xhat_c6, P_c6, e_c6, dy_c6, Pa_c6, Pv_c6, Pw_c6, Phata_c6, Phatv_c6, Phatw_c6, ...
    % Sig_sa_c6, Pdy_c6, Pdyt_c6] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
    %                                      dynarg, datarg, S, C);
    %
    % Store the test case data.
    %save montecarloseed_test_data_case6.mat t_c6 xhat_c6 P_c6 e_c6 dy_c6 Pa_c6 Pv_c6 Pw_c6 ...
    %                                        Phata_c6 Phatv_c6 Phatw_c6 Sig_sa_c6 Pdy_c6 ...
    %                                        Pdyt_c6
    %
    %*********************************************************************************************%


    %
    % Comparisons:
    %

    est{1} = @estseq;
    est{2} = @estsrif;
    
    %
    % Case 1:
    %
    
    % Load the comparison data.
    load montecarloseed_test_data_case1.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations. 
    opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    
    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    for iest = 1:length(est)
        
        % Run the sequential estimator.
        [t_c1n, xhat_c1n, P_c1n, e_c1n, dy_c1n, Pa_c1n, Pv_c1n, Pw_c1n, Phata_c1n, Phatv_c1n, Phatw_c1n, ...
            Sig_sa_c1n, eflag_c1n, Pdy_c1n, Pdyt_c1n] = feval(est{iest},dynfun, datfun, tspan, Xnot, Pnot, ...
            opts, dynarg, datarg, S, C);
        
        
        for a=1:length(t_c1n)
            dt = t_c1n{a} - t_c1{a};
            if(any(dt))
                failed = 1;
                fprintf(2, '\n *** montecarloseed_test, %s case 1: t test failed.\n',char(est{iest}));
            end
        end
        
        dx_max = 5e-9;
        
        for a=1:length(xhat_c1n)
            dx = abs(xhat_c1n{a} - xhat_c1{a});
            if(max(max(dx)) > dx_max)
                failed = 1;
                disp(max(max(dx)))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: xhat test failed.\n',char(est{iest}));
            end
        end
        
        dP_max = 5e-9;
        
        for a=1:length(P_c1n)
            dP = abs(P_c1n{a} - P_c1{a});
            if(max(max(dP)) > dP_max)
                failed = 1;
                disp(max(max(dP)))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: P test failed.\n',char(est{iest}));
            end
        end
        
        de_max = 5e-9;
        
        for a=1:length(e_c1n)
            de = abs(e_c1n{a} - e_c1{a});
            if(max(max(de)) > de_max)
                failed = 1;
                disp(max(max(de)))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: e test failed.\n',char(est{iest}));
            end
        end
        
        ddy_max = 5e-9;
        
        for a=1:length(dy_c1n)
            ddy = abs(dy_c1n{a} - dy_c1{a});
            if(max(max(ddy)) > ddy_max)
                failed = 1;
                disp(max(max(ddy)))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: dy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pa_c1n - Pa_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Pa test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pv_c1n - Pv_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Pv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pw_c1n - Pw_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Pw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phata_c1n - Phata_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Phata test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatv_c1n - Phatv_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Phatv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatw_c1n - Phatw_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Phatw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Sig_sa_c1n - Sig_sa_c1))))) > 1e-6)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Sig_sa test failed.\n',char(est{iest}));
        end
        
        def_max = 5e-9;
        
        for a=1:length(eflag_c1n)
            def = abs(eflag_c1n{a} - eflag_c1{a});
            if(max(max(def)) > def_max)
                failed = 1;
                disp(max(max(def)))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: eflag test failed.\n',char(est{iest}));
            end
        end
        
        dPdy_max = 5e-9;
        
        for a=1:length(Pdy_c1n) 
            dPdy = abs(Pdy_c1n{a} - Pdy_c1{a});
            if(max(max(max(dPdy))) > dPdy_max)
                failed = 1;
                disp(max(max(max(dPdy))))
                fprintf(2, '\n *** montecarloseed_test, %s case 1: Pdy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pdyt_c1n - Pdyt_c1))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 1: Pdyt test failed.\n',char(est{iest}));
        end
    end
        
    %
    % Case 2:
    %
    
    % Load the comparison data.
    load montecarloseed_test_data_case2.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations.
    opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);
    
    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    for iest = 1:length(est)
        
        % Run the sequential estimator.
        [t_c2n, xhat_c2n, P_c2n, e_c2n, dy_c2n, Pa_c2n, Pv_c2n, Pw_c2n, Phata_c2n, Phatv_c2n, Phatw_c2n, ...
            Sig_sa_c2n, eflag_c2n, Pdy_c2n, Pdyt_c2n] = feval(est{iest},dynfun, datfun, tspan, Xnot, Pnot, ...
            opts, dynarg, datarg, S, C);
        
        
        for a=1:length(t_c2n)
            dt = t_c2n{a} - t_c2{a};
            if(any(dt))
                failed = 1;
                fprintf(2, '\n *** montecarloseed_test, %s case 2: t test failed.\n',char(est{iest}));
            end
        end
        
        dx_max = 5e-9;
        
        for a=1:length(xhat_c2n)
            dx = abs(xhat_c2n{a} - xhat_c2{a});
            if(max(max(dx)) > dx_max)
                failed = 1;
                disp(max(max(dx)))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: xhat test failed.\n',char(est{iest}));
            end
        end
        
        dP_max = 5e-9;
        
        for a=1:length(P_c2n)
            dP = abs(P_c2n{a} - P_c2{a});
            if(max(max(dP)) > dP_max)
                failed = 1;
                disp(max(max(dP)))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: P test failed.\n',char(est{iest}));
            end
        end
        
        de_max = 5e-9;
        
        for a=1:length(e_c2n)
            de = abs(e_c2n{a} - e_c2{a});
            if(max(max(de)) > de_max)
                failed = 1;
                disp(max(max(de)))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: e test failed.\n',char(est{iest}));
            end
        end
        
        ddy_max = 5e-9;
        
        for a=1:length(dy_c2n)
            ddy = abs(dy_c2n{a} - dy_c2{a});
            if(max(max(ddy)) > ddy_max)
                failed = 1;
                disp(max(max(ddy)))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: dy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pa_c2n - Pa_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Pa test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pv_c2n - Pv_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Pv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pw_c2n - Pw_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Pw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phata_c2n - Phata_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Phata test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatv_c2n - Phatv_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Phatv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatw_c2n - Phatw_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Phatw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Sig_sa_c2n - Sig_sa_c2))))) > 1e-6)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Sig_sa test failed.\n',char(est{iest}));
        end
        
        def_max = 5e-9;
        
        for a=1:length(eflag_c2n)
            def = abs(eflag_c2n{a} - eflag_c2{a});
            if(max(max(def)) > def_max)
                failed = 1;
                disp(max(max(def)))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: eflag test failed.\n',char(est{iest}));
            end
        end
        
        dPdy_max = 5e-9;
        
        for a=1:length(Pdy_c2n)
            dPdy = abs(Pdy_c2n{a} - Pdy_c2{a});
            if(max(max(max(dPdy))) > dPdy_max)
                failed = 1;
                disp(max(max(max(dPdy))))
                fprintf(2, '\n *** montecarloseed_test, %s case 2: Pdy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pdyt_c2n - Pdyt_c2))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 2: Pdyt test failed.\n',char(est{iest}));
        end
    end
        
    %
    % Case 3:
    %
    
    % Load the comparison data.
    load montecarloseed_test_data_case3.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations.
    opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', [10 20 30]);
    
    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    for iest = 1:length(est)
        
        % Run the sequential estimator.
        [t_c3n, xhat_c3n, P_c3n, e_c3n, dy_c3n, Pa_c3n, Pv_c3n, Pw_c3n, Phata_c3n, Phatv_c3n, Phatw_c3n, ...
            Sig_sa_c3n, eflag_c3n, Pdy_c3n, Pdyt_c3n] = feval(est{iest},dynfun, datfun, tspan, Xnot, Pnot, ...
            opts, dynarg, datarg, S, C);
        
        
        for a=1:length(t_c3n) 
            dt = t_c3n{a} - t_c3{a};
            if(any(dt))
                failed = 1;
                fprintf(2, '\n *** montecarloseed_test, %s case 3: t test failed.\n',char(est{iest}));
            end
        end
        
        dx_max = 5e-9;
        
        for a=1:length(xhat_c3n)
            dx = abs(xhat_c3n{a} - xhat_c3{a});
            if(max(max(dx)) > dx_max)
                failed = 1;
                disp(max(max(dx)))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: xhat test failed.\n',char(est{iest}));
            end
        end
        
        dP_max = 5e-9;
        
        for a=1:length(P_c3n)
            dP = abs(P_c3n{a} - P_c3{a});
            if(max(max(dP)) > dP_max)
                failed = 1;
                disp(max(max(dP)))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: P test failed.\n',char(est{iest}));
            end
        end
        
        de_max = 5e-9;
        
        for a=1:length(e_c3n)
            de = abs(e_c3n{a} - e_c3{a});
            if(max(max(de)) > de_max)
                failed = 1;
                disp(max(max(de)))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: e test failed.\n',char(est{iest}));
            end
        end
        
        ddy_max = 5e-9;
        
        for a=1:length(dy_c3n)
            ddy = abs(dy_c3n{a} - dy_c3{a});
            if(max(max(ddy)) > ddy_max)
                failed = 1;
                disp(max(max(ddy)))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: dy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pa_c3n - Pa_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Pa test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pv_c3n - Pv_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Pv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Pw_c3n - Pw_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Pw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phata_c3n - Phata_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Phata test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatv_c3n - Phatv_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Phatv test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Phatw_c3n - Phatw_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Phatw test failed.\n',char(est{iest}));
        end
        
        if((max(max(max(abs(Sig_sa_c3n - Sig_sa_c3))))) > 1e-6)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Sig_sa test failed.\n',char(est{iest}));
        end
        
        def_max = 5e-9;
        
        for a=1:length(eflag_c3n)
            def = abs(eflag_c3n{a} - eflag_c3{a});
            if(max(max(def)) > def_max)
                failed = 1;
                disp(max(max(def)))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: eflag test failed.\n',char(est{iest}));
            end
        end
        
        dPdy_max = 5e-9;
        
        for a=1:length(Pdy_c3n)
            dPdy = abs(Pdy_c3n{a} - Pdy_c3{a});
            if(max(max(max(dPdy))) > dPdy_max)
                failed = 1;
                disp(max(max(max(dPdy))))
                fprintf(2, '\n *** montecarloseed_test, %s case 3: Pdy test failed.\n',char(est{iest}));
            end
        end
        
        if((max(max(max(abs(Pdyt_c3n - Pdyt_c3))))) > 5e-9)
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, %s case 3: Pdyt test failed.\n',char(est{iest}));
        end
    end

    %
    % Case 4:
    %

    % Load the comparison data.
    load montecarloseed_test_data_case4.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations. 
    opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);

    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    % Run the batch estimator.
    [t_c4n, xhat_c4n, P_c4n, e_c4n, dy_c4n, Pa_c4n, Pv_c4n, Pw_c4n, Phata_c4n, Phatv_c4n, Phatw_c4n, ...
     Sig_sa_c4n, Pdy_c4n, Pdyt_c4n] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
                                          dynarg, datarg, S, C);

    for a=1:length(t_c4n)
        dt = t_c4n{a} - t_c4{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, case 4: t test failed.\n');
        end
    end
   
    dx_max = 2e-8;
    
    for a=1:length(xhat_c4n)
        dx = abs(xhat_c4n{a} - xhat_c4{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** montecarloseed_test, case 4: xhat test failed.\n');
        end
    end

    dP_max = 5e-9;
    
    for a=1:length(P_c4n)
       dP = abs(P_c4n{a} - P_c4{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** montecarloseed_test, case 4: P test failed.\n');
        end
    end

    de_max = 2e-8;
    
    for a=1:length(e_c4n)
        de = abs(e_c4n{a} - e_c4{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** montecarloseed_test, case 4: e test failed.\n');
        end
    end

    ddy_max = 2e-8;
    
    for a=1:length(dy_c4n)
        ddy = abs(dy_c4n{a} - dy_c4{a});
        if(max(max(ddy)) > ddy_max)
            failed = 1;
            disp(max(max(ddy)))
            fprintf(2, '\n *** montecarloseed_test, case 4: dy test failed.\n');
        end
    end

    if((max(max(max(abs(Pa_c4n - Pa_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Pa test failed.\n');
    end
    
    if((max(max(max(abs(Pv_c4n - Pv_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Pv test failed.\n');
    end
    
    if((max(max(max(abs(Pw_c4n - Pw_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Pw test failed.\n');
    end
    
    if((max(max(max(abs(Phata_c4n - Phata_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Phata test failed.\n');
    end
    
    if((max(max(max(abs(Phatv_c4n - Phatv_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Phatv test failed.\n');
    end
    
    if((max(max(max(abs(Phatw_c4n - Phatw_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Phatw test failed.\n');
    end

    if((max(max(max(abs(Sig_sa_c4n - Sig_sa_c4))))) > 1e-6)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Sig_sa test failed.\n');
    end

    dPdy_max = 5e-9;
    
    for a=1:length(Pdy_c4n)
        dPdy = abs(Pdy_c4n{a} - Pdy_c4{a});
        if(max(max(max(dPdy))) > dPdy_max)
            failed = 1;
            disp(max(max(max(dPdy))))
            fprintf(2, '\n *** montecarloseed_test, case 4: Pdy test failed.\n');
        end
    end

    if((max(max(max(abs(Pdyt_c4n - Pdyt_c4))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 4: Pdyt test failed.\n');
    end


    %
    % Case 5:
    %

    % Load the comparison data.
    load montecarloseed_test_data_case5.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations. 
    opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', 10);

    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    % Run the batch estimator.
    [t_c5n, xhat_c5n, P_c5n, e_c5n, dy_c5n, Pa_c5n, Pv_c5n, Pw_c5n, Phata_c5n, Phatv_c5n, Phatw_c5n, ...
     Sig_sa_c5n, Pdy_c5n, Pdyt_c5n] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
                                          dynarg, datarg, S, C);

    for a=1:length(t_c5n)
        dt = t_c5n{a} - t_c5{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, case 5: t test failed.\n');
        end
    end
   
    dx_max = 1e-5;
    
    for a=1:length(xhat_c5n)
        dx = abs(xhat_c5n{a} - xhat_c5{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** montecarloseed_test, case 5: xhat test failed.\n');
        end
    end

    dP_max = 5e-9;
    
    for a=1:length(P_c5n)
       dP = abs(P_c5n{a} - P_c5{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** montecarloseed_test, case 5: P test failed.\n');
        end
    end

    de_max = 1e-5;
    
    for a=1:length(e_c5n) 
        de = abs(e_c5n{a} - e_c5{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** montecarloseed_test, case 5: e test failed.\n');
        end
    end

    ddy_max = 1e-5;
    
    for a=1:length(dy_c5n)
        ddy = abs(dy_c5n{a} - dy_c5{a});
        if(max(max(ddy)) > ddy_max)
            failed = 1;
            disp(max(max(ddy)))
            fprintf(2, '\n *** montecarloseed_test, case 5: dy test failed.\n');
        end
    end

    if((max(max(max(abs(Pa_c5n - Pa_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Pa test failed.\n');
    end
    
    if((max(max(max(abs(Pv_c5n - Pv_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Pv test failed.\n');
    end
    
    if((max(max(max(abs(Pw_c5n - Pw_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Pw test failed.\n');
    end
    
    if((max(max(max(abs(Phata_c5n - Phata_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Phata test failed.\n');
    end
    
    if((max(max(max(abs(Phatv_c5n - Phatv_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Phatv test failed.\n');
    end
    
    if((max(max(max(abs(Phatw_c5n - Phatw_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Phatw test failed.\n');
    end

    if((max(max(max(abs(Sig_sa_c5n - Sig_sa_c5))))) > 1e-6)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Sig_sa test failed.\n');
    end

    dPdy_max = 5e-9;
    
    for a=1:length(Pdy_c5n)
        dPdy = abs(Pdy_c5n{a} - Pdy_c5{a});
        if(max(max(max(dPdy))) > dPdy_max)
            failed = 1;
            disp(max(max(max(dPdy))))
            fprintf(2, '\n *** montecarloseed_test, case 5: Pdy test failed.\n');
        end
    end

    if((max(max(max(abs(Pdyt_c5n - Pdyt_c5))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 5: Pdyt test failed.\n');
    end


    %
    % Case 6:
    %

    % Load the comparison data.
    load montecarloseed_test_data_case6.mat;
    
    % Execute the current processing.
    % Set the number of Monte Carlo cases and the number of update iterations. 
    opts = setOdtbxOptions('MonteCarloCases', 3, 'UpdateIterations', 1);
    
    % Set the Monte Carlo Seed.
    opts = setOdtbxOptions(opts, 'MonteCarloSeed', [10 20 30]);

    % Set the Edit Flag
    opts = setOdtbxOptions(opts,'EditFlag',[2 2]);
    
    % Run the batch estimator.
    [t_c6n, xhat_c6n, P_c6n, e_c6n, dy_c6n, Pa_c6n, Pv_c6n, Pw_c6n, Phata_c6n, Phatv_c6n, Phatw_c6n, ...
     Sig_sa_c6n, Pdy_c6n, Pdyt_c6n] = estbat(dynfun, datfun, tspan_bat, Xnot, Pnot, opts, ...
                                          dynarg, datarg, S, C);

    for a=1:length(t_c6n)
        dt = t_c6n{a} - t_c6{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** montecarloseed_test, case 6: t test failed.\n');
        end
    end
   
    dx_max = 2e-8;
    
    for a=1:length(xhat_c6n)
        dx = abs(xhat_c6n{a} - xhat_c6{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** montecarloseed_test, case 6: xhat test failed.\n');
        end
    end

    dP_max = 5e-9;
    
    for a=1:length(P_c6n)
       dP = abs(P_c6n{a} - P_c6{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** montecarloseed_test, case 6: P test failed.\n');
        end
    end

    de_max = 2e-8;
    
    for a=1:length(e_c6n)
        de = abs(e_c6n{a} - e_c6{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** montecarloseed_test, case 6: e test failed.\n');
        end
    end

    ddy_max = 2e-8;
    
    for a=1:length(dy_c6n)
        ddy = abs(dy_c6n{a} - dy_c6{a});
        if(max(max(ddy)) > ddy_max)
            failed = 1;
            disp(max(max(ddy)))
            fprintf(2, '\n *** montecarloseed_test, case 6: dy test failed.\n');
        end
    end

    if((max(max(max(abs(Pa_c6n - Pa_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Pa test failed.\n');
    end
    
    if((max(max(max(abs(Pv_c6n - Pv_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Pv test failed.\n');
    end
    
    if((max(max(max(abs(Pw_c6n - Pw_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Pw test failed.\n');
    end
    
    if((max(max(max(abs(Phata_c6n - Phata_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Phata test failed.\n');
    end
    
    if((max(max(max(abs(Phatv_c6n - Phatv_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Phatv test failed.\n');
    end
    
    if((max(max(max(abs(Phatw_c6n - Phatw_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Phatw test failed.\n');
    end

    if((max(max(max(abs(Sig_sa_c6n - Sig_sa_c6))))) > 1e-6)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Sig_sa test failed.\n');
    end

    dPdy_max = 5e-9;
    
    for a=1:length(Pdy_c6n)
        dPdy = abs(Pdy_c6n{a} - Pdy_c6{a});
        if(max(max(max(dPdy))) > dPdy_max)
            failed = 1;
            disp(max(max(max(dPdy))))
            fprintf(2, '\n *** montecarloseed_test, case 6: Pdy test failed.\n');
        end
    end

    if((max(max(max(abs(Pdyt_c6n - Pdyt_c6))))) > 5e-9)
        failed = 1;
        fprintf(2, '\n *** montecarloseed_test, case 6: Pdyt test failed.\n');
    end


% End of function.
end
