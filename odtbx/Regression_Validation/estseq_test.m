function failed = estseq_test
%
% estseq_test Regression testing and demo for estseq.m
% See also estseq.m
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
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur             08/29/2012      Extract regression test from
%                                           original estseq.m
%   Ravi Mathur             05/02/2013      Fully extracted regression test
%                                           from original estseq.m

totaltests = 4;
disp('estseq_test: regression testing estseq...')

% Run all the test cases
for k = 1:totaltests,
    disp(['Case ',num2str(k),'...'])
    fail(k,:) = run_test(k);
end

% If system supports parallel processing run all
% the test cases again in parallel
testparallel = 0;
try % See if system supports parallel operation
    if matlabpool('size') == 0,
        matlabpool('open');
        testparallel = 1;
    else
        disp('Skipping parallel test.');
        disp('Self test already running in parallel environment.');
    end
catch %#ok<CTCH>
    disp('Skipping parallel test.');
    disp('No support for parallel processing.');
end
if testparallel,
    for k = 1:totaltests,
        disp(['Parallel Case ',num2str(k),'...'])
        fail(k+totaltests,:) = run_test(k); %#ok<AGROW>
    end
    matlabpool('close');
end

failed = any(any(fail));

end %function

% Runs the actual estbat tests. The contents of this function have been
% extracted from estbat.
function fail = run_test(testnum)
    % Set functions based on the test being performed
    switch testnum
        case 1
            %% Case 1
            % This is a scalar constant state case.  Both the true and the assumed
            % dynamics are the same (i.e., S=1, C=0) and are given below:
            %
            % $$ \begin{array}{ll} \dot{x} = \dot{\hat{x}} = 0, & df/dx = d\hat{f}/d\hat{x} = 0 \end{array}$$
            %
            % The process noise PSD of the truth and the assumed are
            %
            % $$ \begin{array}{ll} Q = 1, & \hat{Q} = 0 \end{array} $$
            %
            % Both the true and the assumed measurement models are the same and are given below:
            %
            % $$ \begin{array}{ll} Y = x, & H = 1 \\
            %                      \hat{Y} = \hat{x}, & \hat{H} = 1\end{array}$$
            %
            % The measurement noise covariance of the truth and the assumed are the
            % same
            %
            % $$ \begin{array}{ll} R = 1, & \hat{R} = 1 \end{array} $$
            %
            % The initial state and covariance for the truth and assumed are
            %
            % $$ \begin{array}{ll} x_o = 0, & P_o = 10\\
            %                      \bar{x}_o = 0, & \bar{P}_o = 10 \end{array} $$
            %
            % The measurement updates occur at t = [1:1:5]. There are 125 Monte Carlo
            % simulations.  In the Monte Carlo simulations, there are 31 interior time
            % steps for the integrator between measurement updates, and each
            % measurement update is iterated 3 times.
            %
            % Since the only difference between the true and the assumed models is in
            % the process noise, we expect the variances due to _a priori_ error and
            % the measurement noise to be the same between the truth and the assumed
            % while the variance due to process noise of the truth to be higher than
            % the assumed.  With no assumed process noise, we expect the total assumed
            % error covariance to decrease with each measurement update.  On the other
            % hand, the true error covarince increases over time since the gains based
            % on the assumed decrease over time and not enough error reduction is
            % performed to cancel the growth during the time update.

            dynfun.tru = @rwdyn;
            dynfun.est = @rwdyn;
            datfun.tru = @rwdat;
            datfun.est = @rwdat;
            load estseq_test1;

        case 2 % Consider covariance
            %% Case 2
            % This is an 8-state case of simple 3-D kinematic equations of motion with
            % the measurement bias (b) and its rate (r) included:
            %
            % $$ \vec{x} = [\begin{array}{cccccccc}x & y & z & u & v & w & b & r\end{array}]^T $$
            %
            % where the true dynamics and its true process noise are given by
            %
            % $$
            % \frac{d}{dt}\left[\begin{array}{c}x\\y\\z\\u\\v\\w\\b\\r\end{array}\right
            % ]
            %                = \left[\begin{array}{cccccccc}
            %                    0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
            %                    \end{array}\right]\vec{x} $$
            %
            % $$
            % Q = \left[\begin{array}{cccccccc}
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 1e-6 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 1e-6 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 1e-6 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 1e-6
            %                    \end{array}\right] $$
            %
            % The true measurement noise model is given by
            %
            % $$ Y = \left[\begin{array}{cccccccc}
            %        1 & 0 & 0 & 0 & 0 & 0 & -1 & 0 \\
            %        0 & 1 & 0 & 0 & 0 & 0 & -1 & 0 \\
            %        0 & 0 & 1 & 0 & 0 & 0 & -1 & 0
            %        \end{array}\right]\vec{x} $$
            %
            % $$ R = \left[\begin{array}{ccc}
            %        1 & 0 & 0 \\
            %        0 & 1 & 0 \\
            %        0 & 0 & 1
            %        \end{array}\right] $$
            %
            % In the estimator, the dynamics model and its process noise consist only
            % of the solve-for states corresponding to the 3-D position and velocity.
            %
            % $$ \vec{x} = [\begin{array}{cccccccc}x & y & z & u & v & w\end{array}]^T $$
            %
            % $$
            % \frac{d}{dt}\left[\begin{array}{c}x\\y\\z\\u\\v\\w\end{array}\right
            % ]
            %                = \left[\begin{array}{cccccc}
            %                    0 & 0 & 0 & 1 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 1 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 1 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0
            %                    \end{array}\right]\vec{x} $$
            %
            % $$
            % Q = \left[\begin{array}{cccccc}
            %                    0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 0 \\
            %                    0 & 0 & 0 & 1e-2 & 0 & 0 \\
            %                    0 & 0 & 0 & 0 & 1e-2 & 0 \\
            %                    0 & 0 & 0 & 0 & 0 & 1e-2
            %                    \end{array}\right] $$
            %
            % The estimator measurement model also consists of only the solve-for
            % states:
            %
            % $$ Y = \left[\begin{array}{cccccc}
            %        1 & 0 & 0 & 0 & 0 & 0 \\
            %        0 & 1 & 0 & 0 & 0 & 0 \\
            %        0 & 0 & 1 & 0 & 0 & 0
            %        \end{array}\right]\vec{x} $$
            %
            % $$ R = \left[\begin{array}{ccc}
            %        1e4 & 0 & 0 \\
            %        0 & 1e4 & 0 \\
            %        0 & 0 & 1e4
            %        \end{array}\right] $$
            %
            % The partition matrices are
            %
            % $$ S = \left[\begin{array}{cccccccc}
            %              1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
            %              0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
            %              0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
            %              0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
            %              0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
            %              0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\end{array}\right] $$
            %
            % $$ C = \left[\begin{array}{cccccccc}
            %              0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
            %              0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\end{array}\right] $$
            %
            % The initial states and covariances are
            %
            % $$ x_o = [\begin{array}{cccccccc}0 & 0 & 0 & 0 & 0 & 0 & 0 &
            % 0\end{array}]^T$$
            %
            % $$ P_o = \left[\begin{array}{cc}
            %                I_{6x6} & O_{6x2}\\
            %                O_{2x6} & 2I_{2x2}\end{array}\right] $$
            %
            % $$ \bar{x}_o = [\begin{array}{cccccc}0 & 0 & 0 & 0 & 0 & 0\end{array}]^T$$
            %
            % $$ \bar{P}_o = 2I_{6x6} $$
            %
            % The measurement updates occur at t = [1:1:30]. There are 12 Monte Carlo
            % simulations.  In the Monte Carlo simulations, there are 3 interior time
            % steps for the integrator between measurement updates, and each
            % measurement update is iterated 3 times.
            
            dynfun.tru = @irwbdyn;
            dynfun.est = @irwdyn;
            datfun.tru = @irwbdat;
            datfun.est = @irwdat;
            load estseq_test2;

        case 3 % Schmidt Kalman filter version of case 2
            %% Case 3
            % Case 3 is a Schmidt-Kalman filter implementation of Case 2.  Both the true
            % and the estimator models are identical corresponding to the true model of
            % Case 2 above.
            
            dynfun.tru = @irwbdyn;
            dynfun.est = @irwbdyn;% Modified from test2
            datfun.tru = @irwbdat;
            datfun.est = @irwbdat;% Modified from test2
            load estseq_test3;

        case 4 % Estimation using JAT forces
            %% Case 4
            % This regression test case needs to be documented.
            dynfun.tru = @jatForces_km;
            dynfun.est = dynfun.tru;
            datfun.tru = @gsmeas;
            datfun.est = datfun.tru;
            load estseq_test4;
            
            % Create JAT structures and initialize integrator
            dynarg.tru = createJATWorld(options); 
            dynarg.est = dynarg.tru;
            gsList  = createGroundStationList('DBS_NDOSL_WGS84_Mod_Example.txt');
            options = setOdtbxOptions(options,'gsList',gsList);
            datarg.tru  = options;
            datarg.est  = setOdtbxOptions(datarg.tru,'rSigma',3*[1e-3 1e-6]);
            
            options = setOdtbxOptions(options,'OdeSolver',@ode113,'OdeSolvOpts',...
                odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10));
            
% The estseq_test4.mat file was creating using the following data, in
% addition to the truth test results data. This comment can be deleted
% pending results verification.
%             tspan = 0:10:300;
%             S = eye(6);            % Solve-for map - solve for all 6 states
%             Xo = [6878;0.00;0.00;0.00;0.00;8.339];       % km & km/sec                     
%             Xbaro = S*Xo;
%             Po = diag([1e-2 1e-4 1e-1 1e-4 1e-7 1e-5].^2); % km^2 & km^2/s^2
%             Pbaro = S*Po*S';
%             C=[];
% 
%             options = odtbxOptions;
%             options = setOdtbxOptions(options, 'epoch', JATConstant('MJDJ2000') );
%             options = setOdtbxOptions(options, 'cD', 2.2);
%             options = setOdtbxOptions(options, 'cR', 0.7);
%             options = setOdtbxOptions(options, 'mass', 1000);
%             options = setOdtbxOptions(options, 'draga', 20, 'srpArea', 20);
%             options = setOdtbxOptions(options, 'earthGravityModel', '2body');
%             options = setOdtbxOptions(options, 'gravDeg', 2, 'gravOrder', 2);
%             options = setOdtbxOptions(options, 'useSolarGravity', false);
%             options = setOdtbxOptions(options, 'useLunarGravity', false);
%             options = setOdtbxOptions(options, 'useSolarRadiationPressure', false);
%             options = setOdtbxOptions(options, 'useAtmosphericDrag', false);
%             options = setOdtbxOptions(options, 'atmosphereModel', 'HP');
%             options = setOdtbxOptions(options, 'nParameterForHPModel', 2);
%             options = setOdtbxOptions(options, 'f107Daily', 150);
%             options = setOdtbxOptions(options, 'f107Average', 150);
%             options = setOdtbxOptions(options, 'ap', 15);
% 
%             options = setOdtbxOptions(options,'epoch',datenum(2006,12,31,23,59,38.3));
%             options = setOdtbxOptions(options,'useRange',true);
%             options = setOdtbxOptions(options,'rangeType','2way');
%             options = setOdtbxOptions(options,'useRangerate',true);
%             options = setOdtbxOptions(options,'useDoppler',false);
%             options = setOdtbxOptions(options,'rSigma',[1e-3 1e-6]);
%             options = setOdtbxOptions(options,'useTroposphere',false);
%             options = setOdtbxOptions(options,'gsID',{'ZZOD'});
%             options = setOdtbxOptions(options,'gsElevationConstraint',0);
% 
%             options = setOdtbxOptions(options, 'UpdateIterations',2,'MonteCarloCases',2,...
%                                                'MonteCarloSeed', 1, 'refint',3,'EditFlag',[2 2]);
    end
    
    %% Call estseq with loaded data
    [~,~,Phat,e,~,Pa,Pv,Pw,Phata,Phatv,Phatw,Sig_sa,~,~,~,Pm,~,~,S,C] = ...
        estseq(dynfun, datfun, tspan, Xo, Po, options, dynarg, datarg, S, C);

    %% Get needed Monte Carlo run information
    ncases = getOdtbxOptions(options,'MonteCarloCases',1);
    ischmidt = getOdtbxOptions(options,'SchmidtKalman',0);
    [~,~,lentr] = size(Sig_sa_test);
    ns = size(S(:,:,1),1);
    
    %% Check results against accepted values
    [dPhat{1:ncases}] = deal(NaN(size(Phat)));
    [de{1:ncases}] = deal(NaN(size(e)));
    Pfail = zeros(lentr,ncases);
    efail = zeros(lentr,ncases,ns);
    for k = ncases:-1:1,
        dPhat{k} = Phat_test{k} - Phat{k}; %#ok<USENS>
        %sPhat{k} = Phat_test{k} + Phat{k}; %#ok<AGROW>
        de{k} = e_test{k} - e{k}; %#ok<USENS>
        % Is each error sample within prob of 1e-9 of its corresponding test
        % value?  Is each Phat sample "close enough," in terms of the
        % matrix 2-norm (largest singular value) to its test value?
        try
            chistat = chi2inv(1 - 1.0e-9,ns);
            vectest = true;
        catch %#ok<CTCH>
            chistat = 37.325;
            vectest = false;
        end
        for i = lentr:-1:1,
            SPS = S(:,:,i)*unscrunch(P_test(:,i))*S(:,:,i)';  
            if ischmidt==1
                dek=S(:,:,i)*de{k}(:,i); % Only check the solved-for errors
            else
                dek=de{k}(:,i);
            end
            if vectest,
                efail(i,k,1) = ...
                    dek'/SPS*dek > chistat;
                    %dek'*inv(SPS)*dek > chistat; 
            else
                for j = ns:-1:1,
                    efail(i,k,j) = dek(j)^2/SPS(j,j) > 37.325; 
                end
            end
            Pfail(i,k) = ...
                norm(unscrunch(dPhat{k}(:,i))) > ...
                0.1*norm(unscrunch(Phat_test{k}(:,i))); 
        end
    end
    P = scrunch(Pa + Pv + Pw + Pm);
    fail = logical([...
        any(any(any(abs(P_test - P) > 4e-10))), ...
        any(any(any(abs(Pa_test - Pa) > 1e-10))), ...
        any(any(any(abs(Pv_test - Pv) > 1e-10))),...
        any(any(any(abs(Pw_test - Pw) > 4e-10))), ...
        any(any(any(abs(Phatw_test - Phatw) > 1e-10))), ...
        any(any(any(abs(Phatv_test - Phatv) > 1e-10))), ...
        any(any(any(abs(Phata_test - Phata) > 1e-10))), ...
        any(any(any(abs(Sig_sa_test - Sig_sa) > 2e-10))), ...
        any(any(any(efail))), ...
        any(any(Pfail))]);
end

%% Self-test user functions

%% Dynamics function for case 1: true and estimate
function [Xdot,A,Q] = rwdyn(t,X,q)
el = length(t);
A = zeros(1,1,el);
Xdot = zeros(size(X));
Q = q*ones(1,1,el);
end

%% Dynamics function for case 2: true
function [Xdot,A,Q] = irwbdyn(t,X,q)
el = length(t);
A([1:3 7:8],[4:6 8]) = eye(5,4);
Xdot = A*X;
A = repmat(A,[1 1 el]);
Q([4:6 8],[4:6 8]) = q*eye(4);
Q = repmat(Q,[1 1 el]);
end

%% Dynamics function for case 2: estimate, solve for only
function [Xdot,A,Q] = irwdyn(t,X,q)
el = length(t);
A(:,4:6) = eye(6,3);
Xdot = A*X;
A = repmat(A,[1 1 el]);
Q(4:6,4:6) = q*eye(3);
Q = repmat(Q,[1 1 el]);
end

%% Measurement function for case 1: true and estimate
function [Y,H,R] = rwdat(t,X,r)
Y = X;
H = ones(1,1,length(t));
R = r*ones(1,1,length(t));
end

%% Measurement function for case 2: true
function [Y,H,R] = irwbdat(t,X,r)
el = length(t);
H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
Y = H*X;
H = repmat(H,[1 1 el]);
R = r*repmat(eye(3),[1 1 el]);
end

%% Measurement function for case 2: estimate, solve for only
function [Y,H,R] = irwdat(t,X,r)
el = length(t);
H = [eye(3) zeros(3,3)];
Y = H*X;
H = repmat(H,[1 1 el]);
R = r*repmat(eye(3),[1 1 el]);
end