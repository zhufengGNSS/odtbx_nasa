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

function [failed] = estbat_no_proc_noise_test

    % Default to passed.
    failed = 0;

    %
    % This code is taken from estbat_example.m.
    %

    % The next two variables are inputs to the dynamics and measurement models:
    mu = 3.986004415e5; % km^3'sec^2
    sig = diag(1000*[1e-4,1e-7]); % km & km/s

    % Run the estimator over a 5 minute pass at 1 second intervals:
    tspan = 0:1:300;

    % The next two variables define the initial reference state and a priori
    % covariance matrix:
    xo = [6878;0.00;0.00;0.00;0.00;8.339]; %km & km/s
    Po = diag([1e-2 1e-4 1e-1 1e-4 1e-7 1e-5].^2); % km^2 & km^2/s^2

    % Set some options, including the option to not use process noise:
    opts = setOdtbxOptions('MonteCarloCases',30,'UpdateIterations',6, 'UseProcNoise', 0);

    % Set up the solve-for and consider mapping, and the true and formal
    % dynamics and measurement models:
    I6 = eye(6);
    
    S = I6; % Solve-for map; this says solve-for all 6 states
    
    C = zeros(0,0,length(tspan)); % Consider map; this says no consider states
    
    % The next four lines say to use the same models for truth and estimator
    % (but note that the process noise in r2bp is ignored by the batch):
    dynfun.tru = @r2bp;
    datfun.tru = @rrdot3D;
    dynfun.est = @r2bp;
    datfun.est = @rrdot3D;
    
    % The next four lines say to use the same a prior for truth and estimator:
    Pnot.Po = Po;
    Pnot.Pbaro = S*Po*S';
    Xnot.Xo = xo;
    Xnot.Xbaro = S*xo;
    
    % The next two lines say the estimator will assume 3x more measurement
    % noise than the truth:
    datarg.tru = sig;
    datarg.est = 3*sig;
    
    % Now run the batch estimator:
    [t, xhat, Phat, e, y, Pa, Pv, Pw, Phata, Phatv, Phatw]=estbat(...
        dynfun, datfun, tspan, Xnot, Pnot, opts, mu, datarg, S, C);
    
    %
    % Code block to store data for regression testing:
    %
    %t_b = t;
    %xhat_b = xhat;
    %Phat_b = Phat;
    %e_b = e;
    %y_b = y;
    %Pa_b = Pa;
    %Pv_b = Pv;
    %Pw_b = Pw;
    %Phata_b = Phata;
    %Phatv_b = Phatv;
    %Phatw_b = Phatw;
    %
    %save estbat_no_proc_noise_test.mat t_b xhat_b Phat_b e_b y_b Pa_b Pv_b Pw_b Phata_b Phatv_b Phatw_b
    
    % Load data for regression testing:
    load estbat_no_proc_noise_data.mat;

    
    for(a=1:length(t))
        dt = t{a} - t_b{a};
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** estbat_no_proc_noise_test: t test failed.\n');
        end
    end
    
    dx_max = 0.6;
    
    for(a=1:length(xhat))
        dx = abs(xhat{a} - xhat_b{a});
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** estbat_no_proc_noise_test: xhat test failed.\n');
        end
    end
    
    dP_max = 3.5e-7;
    
    for(a=1:length(Phat))
        dP = abs(Phat{a} - Phat_b{a});
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** estbat_no_proc_noise_test: Phat test failed.\n');
        end
    end
    
    de_max = 0.85;
    
    for(a=1:length(e))
        de = abs(e{a} - e_b{a});
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** estbat_no_proc_noise_test: e test failed.\n');
        end
    end
    
    dy_max = 0.75;
    
    for(a=1:length(y))
        dy = abs(y{a} - y_b{a});
        if(max(max(dy)) > dy_max)
            failed = 1;
            disp(max(max(dy)))
            fprintf(2, '\n *** estbat_no_proc_noise_test: y test failed.\n');
        end
    end
    
    if((max(max(max(abs(Pa - Pa_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Pa test failed.\n');
    end
    
    if((max(max(max(abs(Pv - Pv_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Pv test failed.\n');
    end
    
    if((max(max(max(abs(Pw - Pw_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Pw test failed.\n');
    end
    
    if((max(max(max(abs(Phata - Phata_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Phata test failed.\n');
    end
    
    if((max(max(max(abs(Phatv - Phatv_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Phatv test failed.\n');
    end
    
    if((max(max(max(abs(Phatw - Phatw_b))))) > eps)
        failed = 1;
        fprintf(2, '\n *** estbat_no_proc_noise_test: Phatw test failed.\n');
    end
    
end







