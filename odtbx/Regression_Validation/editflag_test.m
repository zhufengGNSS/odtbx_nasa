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
% Sun Hur-Diaz (Added vector editing test)

% Begin function.
function [failed] = editflag_test

% Default to pass.
failed = 0;
warning('off', 'ODTBX:COVSMPL:seedReset');

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

optsnom = opts;

%
% Block of code for generating and storing test data sets.
%
% Set measurement editing options.
%opts = setOdtbxOptions(opts, 'EditRatio', [0.5 1], 'EditFlag', [1 1]);
%
% Run the sequential estimator.
%[t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw, Sig_sa, eflag] = estseq(...
%    dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
%
%eflag_s_case_1 = eflag;
%
% Set measurement editing options.
%opts = setOdtbxOptions(opts, 'EditRatio', [0.5 1], 'EditFlag', [2 2]);
%
% Run the sequential estimator.
%[t, xhat, P, e, dy, Pa, Pv, Pw, Phata, Phatv, Phatw, Sig_sa, eflag] = estseq(...
%    dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
%
%eflag_s_case_2 = eflag;
%
%save eflag_test_data.mat eflag_s_case_1 eflag_s_case_2
%

est{1} = @estseq;
est{2} = @estsrif;

% Load the comparison data.
load eflag_test_data.mat;

for iest = 1:length(est)
    
    % Case 1.
    % Set measurement editing options.
    opts = setOdtbxOptions(optsnom, 'EditRatio', [0.5 1], 'EditFlag', [1 1]);
    
    % Run the sequential estimator.
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, eflag_case_1] = feval(est{iest},...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    
    % Case 2.
    % Set measurement editing options.
    opts = setOdtbxOptions(optsnom, 'EditRatio', [0.5 1], 'EditFlag', [2 2]);
    
    % Run the sequential estimator.
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, eflag_case_2] = feval(est{iest},...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
   
    % Case 3 - Vector Editing
    % Set measurement editing options.
    opts = setOdtbxOptions(optsnom, 'EditRatio', [0.5], 'EditVector', 1);
    
    % Run the sequential estimator.
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, eflag_case_3] = feval(est{iest},...
        dynfun, datfun, tspan, Xnot, Pnot, opts, dynarg, datarg, S, C);
    
    %
    % Comparisons.
    %
    if(size(eflag_case_1) ~= size(eflag_s_case_1))
        failed = 1;
        fprintf('Editflag test: %s case 1 eflag size comparison failed.\n',char(est{iest}));
    end
    
    if(size(eflag_case_2) ~= size(eflag_s_case_2))
        failed = 1;
        fprintf('Editflag test: %s case 2 eflag size comparison failed.\n',char(est{iest}));
    end
    
    n = size(eflag_case_1{1});
    
    if(n ~= size(eflag_s_case_1{1}))
        failed = 1;
        fprintf('Editflag test: %s case 1 eflag array size comparison failed.\n',char(est{iest}));
    end
    
    n2 = size(eflag_case_1{1}, 2);
    
    for a=1:n2 
        for b=1:2
            if(isnan(eflag_case_1{1}(b,a)))
                if(~isnan(eflag_s_case_1{1}(b,a)))
                    failed = 1;
                    fprintf('Editflag test: %s case 1 eflag array value mismatch (NaNs).\n',char(est{iest}));
                end
            else
                if(eflag_case_1{1}(b,a) ~= eflag_s_case_1{1}(b,a))
                    failed = 1;
                    fprintf('Editflag test: %s case 1 eflag array value mismatch.\n',char(est{iest}));
                end
            end
        end
    end
    
    n = size(eflag_case_2{1});
    
    if(n ~= size(eflag_s_case_2{1}))
        failed = 1;
        fprintf('Editflag test: %s case 2 eflag array size comparison failed.\n',char(est{iest}));
    end
    
    n2 = size(eflag_case_2{1}, 2);
    
    for a=1:n2 
        for b=1:2 
            if(isnan(eflag_case_2{1}(b,a)))
                if(~isnan(eflag_s_case_2{1}(b,a)))
                    failed = 1;
                    fprintf('Editflag test: %s case 2 eflag array value mismatch (NaNs).\n',char(est{iest}));
                end
            else
                if(eflag_case_2{1}(b,a) ~= eflag_s_case_2{1}(b,a))
                    failed = 1;
                    fprintf('Editflag test: %s case 2 eflag array value mismatch.\n',char(est{iest}));
                end
            end
        end
    end
    
    n = size(eflag_case_3{1});
    
    if(n ~= size(eflag_s_case_3{1}))
        failed = 1;
        fprintf('Editflag test: %s case 3 eflag array size comparison failed.\n',char(est{iest}));
    end
    
    n2 = size(eflag_case_3{1}, 2);
    
    for a=1:n2 
        for b=1:2 
            if(isnan(eflag_case_3{1}(b,a)))
                if(~isnan(eflag_s_case_3{1}(b,a)))
                    failed = 1;
                    fprintf('Editflag test: %s case 3 eflag array value mismatch (NaNs).\n',char(est{iest}));
                end
            else
                if(eflag_case_3{1}(b,a) ~= eflag_s_case_3{1}(b,a))
                    failed = 1;
                    fprintf('Editflag test: %s case 3 eflag array value mismatch.\n',char(est{iest}));
                end
            end
        end
    end

end

% End of function.
end
