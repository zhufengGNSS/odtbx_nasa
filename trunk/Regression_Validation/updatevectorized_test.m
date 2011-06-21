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

function [failed] = updatevectorized_test

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

% Specify the time vector at which the measurements are taken.
tspan = (0:60:(4*3600));

% Set some options for the sequential filter.
opts = setOdtbxOptions('MonteCarloCases', 1, 'UpdateIterations', 3);

% Set the Monte Carlo Seed.
opts = setOdtbxOptions(opts, 'MonteCarloSeed', 691);

% Set measurement editing options.
opts = setOdtbxOptions(opts, 'EditRatio', [9 9], 'EditFlag', [1 1]);

%*********************************************************************************************
% Block of code for generating regression test data sets.
%
%=============================================%
%=  Case 1: UpdateVectorized flag set to NO. =%
%=============================================%
%
% Set the update vectorized flag to No.
%opts = setOdtbxOptions(opts, 'UpdateVectorized', 0);
%
% Run the sequential estimator with update vectorized set to No.
%[t_c1, xhat_c1, P_c1, e_c1, dy_c1, Pa_c1, Pv_c1, Pw_c1, Phata_c1, Phatv_c1, ...
% Phatw_c1, Sig_sa_c1, eflag_c1, Pdy_c1, Pdyt_c1] = estseq(dynfun, datfun, ...
% tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
%
% Store the data for test case 1.
%save updatevectorized_test_data_case1.mat t_c1 xhat_c1 P_c1 e_c1 dy_c1 Pa_c1 Pv_c1 Pw_c1 ...
%                                          Phata_c1 Phatv_c1 Phatw_c1 Sig_sa_c1 eflag_c1 ...
%                                          Pdy_c1 Pdyt_c1
%
%=============================================%
%= Case 2: UpdateVectorized flag set to YES. =%
%=============================================%
%
% Set the update vectorized flag to Yes.
%opts = setOdtbxOptions(opts, 'UpdateVectorized', 1);
%
% Run the sequential estimator with update vectorized set to Yes.
%[t_c2, xhat_c2, P_c2, e_c2, dy_c2, Pa_c2, Pv_c2, Pw_c2, Phata_c2, Phatv_c2, ...
% Phatw_c2, Sig_sa_c2, eflag_c2, Pdy_c2, Pdyt_c2] = estseq(dynfun, datfun, ...
% tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
%
% Store the data for test case 2.
%save updatevectorized_test_data_case2.mat t_c2 xhat_c2 P_c2 e_c2 dy_c2 Pa_c2 Pv_c2 Pw_c2 ...
%                                          Phata_c2 Phatv_c2 Phatw_c2 Sig_sa_c2 eflag_c2 ...
%                                          Pdy_c2 Pdyt_c2
%
%*********************************************************************************************

%
% Comparisons:
%

est{1} = @estseq;
est{2} = @estsrif;

%
% Case 1:
%

% Load the comparison data.
load updatevectorized_test_data_case1.mat;

% Execute the current processing.
% Set the update vectorized flag to No.
opts = setOdtbxOptions(opts, 'UpdateVectorized', 0);

for iest = 1:length(est)
    
    % Run the sequential estimator.
    [t_c1n, xhat_c1n, P_c1n, e_c1n, dy_c1n, Pa_c1n, Pv_c1n, Pw_c1n, Phata_c1n, Phatv_c1n, Phatw_c1n, ...
        Sig_sa_c1n, eflag_c1n, Pdy_c1n, Pdyt_c1n] = feval(est{iest},dynfun, datfun, tspan, Xnot, Pnot, ...
        opts, dynarg, datarg, S, C);
    
    for a=1:length(t_c1n)
        dt = t_c1n{a} - t_c1{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** updatevectorized_test, %s case 1: t test failed.\n',char(est{iest}));
        end
    end
    
    dx_max = 2e-9;
    
    for a=1:length(xhat_c1n)
        dx = abs(xhat_c1n{a} - xhat_c1{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: xhat test failed.\n',char(est{iest}));
        end
    end
    
    dP_max = 2e-9;
    
    for a=1:length(P_c1n)
        dP = abs(P_c1n{a} - P_c1{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: P test failed.\n',char(est{iest}));
        end
    end
    
    de_max = 2e-9;
    
    for a=1:length(e_c1n)
        de = abs(e_c1n{a} - e_c1{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: e test failed.\n',char(est{iest}));
        end
    end
    
    ddy_max = 2e-9;
    
    for a=1:length(dy_c1n)
        ddy = abs(dy_c1n{a} - dy_c1{a});
        if(max(max(ddy)) > ddy_max)
            failed = 1;
            disp(max(max(ddy)))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: dy test failed.\n',char(est{iest}));
        end
    end
    
    if((max(max(max(abs(Pa_c1n - Pa_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Pa test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Pv_c1n - Pv_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Pv test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Pw_c1n - Pw_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Pw test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phata_c1n - Phata_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Phata test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phatv_c1n - Phatv_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Phatv test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phatw_c1n - Phatw_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Phatw test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Sig_sa_c1n - Sig_sa_c1))))) > 1e-6)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Sig_sa test failed.\n',char(est{iest}));
    end
    
    def_max = 2e-9;
    
    for a=1:length(eflag_c1n)
        def = abs(eflag_c1n{a} - eflag_c1{a});
        if(max(max(def)) > def_max)
            failed = 1;
            disp(max(max(def)))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: eflag test failed.\n',char(est{iest}));
        end
    end
    
    dPdy_max = 2e-9;
    
    for a=1:length(Pdy_c1n)
        dPdy = abs(Pdy_c1n{a} - Pdy_c1{a});
        if(max(max(max(dPdy))) > dPdy_max)
            failed = 1;
            disp(max(max(max(dPdy))))
            fprintf(2, '\n *** updatevectorized_test, %s case 1: Pdy test failed.\n',char(est{iest}));
        end
    end
    
    if((max(max(max(abs(Pdyt_c1n - Pdyt_c1))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 1: Pdyt test failed.\n',char(est{iest}));
    end
    
end

%
% Case 2:
%

% Load the comparison data.
load updatevectorized_test_data_case2.mat;

% Execute the current processing.
% Set the update vectorized flag to Yes.
opts = setOdtbxOptions(opts, 'UpdateVectorized', 1);

for iest = 1:length(est)
    
    % Run the sequential estimator.
    [t_c2n, xhat_c2n, P_c2n, e_c2n, dy_c2n, Pa_c2n, Pv_c2n, Pw_c2n, Phata_c2n, Phatv_c2n, Phatw_c2n, ...
        Sig_sa_c2n, eflag_c2n, Pdy_c2n, Pdyt_c2n] = feval(est{iest},dynfun, datfun, tspan, Xnot, Pnot, ...
        opts, dynarg, datarg, S, C);
    
    
    for a=1:length(t_c2n)
        dt = t_c2n{a} - t_c2{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** updatevectorized_test, %s case 2: t test failed.\n',char(est{iest}));
        end
    end
    
    dx_max = 2e-9;
    
    for a=1:length(xhat_c2n)
        dx = abs(xhat_c2n{a} - xhat_c2{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: xhat test failed.\n',char(est{iest}));
        end
    end
    
    dP_max = 2e-9;
    
    for a=1:length(P_c2n)
        dP = abs(P_c2n{a} - P_c2{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: P test failed.\n',char(est{iest}));
        end
    end
    
    de_max = 2e-9;
    
    for a=1:length(e_c2n)
        de = abs(e_c2n{a} - e_c2{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: e test failed.\n',char(est{iest}));
        end
    end
    
    ddy_max = 2e-9;
    
    for a=1:length(dy_c2n)
        ddy = abs(dy_c2n{a} - dy_c2{a});
        if(max(max(ddy)) > ddy_max)
            failed = 1;
            disp(max(max(ddy)))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: dy test failed.\n',char(est{iest}));
        end
    end
    
    if((max(max(max(abs(Pa_c2n - Pa_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Pa test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Pv_c2n - Pv_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Pv test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Pw_c2n - Pw_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Pw test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phata_c2n - Phata_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Phata test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phatv_c2n - Phatv_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Phatv test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Phatw_c2n - Phatw_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Phatw test failed.\n',char(est{iest}));
    end
    
    if((max(max(max(abs(Sig_sa_c2n - Sig_sa_c2))))) > 1e-6)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Sig_sa test failed.\n',char(est{iest}));
    end
    
    def_max = 2e-9;
    
    for a=1:length(eflag_c2n)
        def = abs(eflag_c2n{a} - eflag_c2{a});
        if(max(max(def)) > def_max)
            failed = 1;
            disp(max(max(def)))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: eflag test failed.\n',char(est{iest}));
        end
    end
    
    dPdy_max = 2e-9;
    
    for a=1:length(Pdy_c2n)
        dPdy = abs(Pdy_c2n{a} - Pdy_c2{a});
        if(max(max(max(dPdy))) > dPdy_max)
            failed = 1;
            disp(max(max(max(dPdy))))
            fprintf(2, '\n *** updatevectorized_test, %s case 2: Pdy test failed.\n',char(est{iest}));
        end
    end
    
    if((max(max(max(abs(Pdyt_c2n - Pdyt_c2))))) > 2e-9)
        failed = 1;
        fprintf(2, '\n *** updatevectorized_test, %s case 2: Pdyt test failed.\n',char(est{iest}));
    end
end

end
