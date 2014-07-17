function failed = estsrif_test
%
% estsrif_test Regression testing and demo for estsrif.m
% See also estsrif.m
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
%
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
%   Ravi Mathur             08/28/2012      Extract regression test from
%                                           original estsrif.m
%
%
%
% Note that the brunt of the regression testing still occurs in
% estsrif.m. It is deeply embedded and will require considerable effort
% to extract all of it into estsrif_test.m.

totaltests = 3;
disp('estsrif_test: regression testing estsrif...')
warning('off', 'ODTBX:COVSMPL:seedReset');

% Run all the test cases
for k = 1:totaltests,
    disp(['Case ',num2str(k),'...'])
    fail(k,:) = run_test(k); %#ok<AGROW>
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
end % function

% Runs the actual estsrif tests. The contents of this function have been
% extracted from estsrif.
function fail = run_test(testnum)

    switch testnum
        case 1
            dynfun.tru = @rwdyn;
            dynfun.est = @rwdyn;
            datfun.tru = @rwdat;
            datfun.est = @rwdat;
            load estseq_test1;
            options = setOdtbxOptions(options,'UseSmoother',1);

        case 2 % Consider covariance
            dynfun.tru = @irwbdyn;
            dynfun.est = @irwdyn;
            datfun.tru = @irwbdat;
            datfun.est = @irwdat;
            load estseq_test2;
            options = setOdtbxOptions(options,'UseSmoother',1);

        case 3 % Estimation using JAT forces
            dynfun.tru = @jatForces_km;
            dynfun.est = dynfun.tru;
            datfun.tru = @gsmeas;
            datfun.est = datfun.tru;
            load estseq_test4;
            
            % Create JAT structures and initialize integrator
            gsList  = createGroundStationList('DBS_NDOSL_WGS84_Mod_Example.txt');
            options = setOdtbxOptions(options,'gsList',gsList);
            options = setOdtbxOptions(options,'OdeSolver',@ode113,'OdeSolvOpts',...
                odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10));
            options = setOdtbxOptions(options,'UseSmoother',1);
            dynarg.tru = createJATWorld(options);
            dynarg.est = dynarg.tru;
            datarg.tru  = options;
            datarg.est  = setOdtbxOptions(datarg.tru,'rSigma',3*[1e-3 1e-6]);
    end

    %% Call estsrif with loaded data
    [~,~,Phat,e,~,Pa,Pv,Pw,Phata,Phatv,Phatw,Sig_sa,~,~,~,~,Pstar,estar,Pm,~,~,S,C] = ...
        estsrif(dynfun, datfun, tspan, Xo, Po, options, dynarg, datarg, S, C);
    
    %% Get needed Monte Carlo run information
    ncases = getOdtbxOptions(options,'MonteCarloCases',1);
    [~,~,lentr] = size(Sig_sa_test);
    ns = size(S(:,:,1),1);    
    
    %% Check results against accepted values
    [dPhat{1:ncases}] = deal(NaN(size(Phat)));
    [dPstar{1:ncases}] = deal(NaN(size(Pstar)));
    [de{1:ncases}] = deal(NaN(size(e)));
    [destar{1:ncases}] = deal(NaN(size(estar)));
    Pfail = zeros(lentr,ncases);
    Pstarfail = zeros(lentr,ncases);
    efail = zeros(lentr,ncases,ns);
    estarfail = zeros(lentr,ncases,ns);
    
    for k = ncases:-1:1,
        dPhat{k} = Phat_test{k} - Phat{k}; %#ok<USENS>
        dPstar{k} = Pstar_test{k} - Pstar{k};
        de{k} = e_test{k} - e{k}; %#ok<USENS>
        destar{k} = estar_test{k} - estar{k}; %#ok<USENS>
        % Is each error sample within prob of 1e-9 of its corresponding test
        % value?  Is each Phat sample "close enough," in terms of the
        % matrix 2-norm (largest singular value) to its test value?
        try
            chistat = chi2inv(1 - 1e-9,ns);
            vectest = true;
        catch %#ok<CTCH>
            chistat = 37.325;
            vectest = false;
        end
        for i = lentr:-1:1,
            SPS = S(:,:,i)*unscrunch(P_test(:,i))*S(:,:,i)';  
            
            dek=de{k}(:,i);
            dekstar=destar{k}(:,i);
            
            if vectest,
                efail(i,k,1) = ...
                    dek'/SPS*dek > chistat;
                estarfail(i,k,1) = ...
                    dekstar'/SPS*dekstar > chistat;
            else
                for j = ns:-1:1,
                    efail(i,k,j) = dek(j)^2/SPS(j,j) > 37.325; 
                    estarfail(i,k,j) = dekstar(j)^2/SPS(j,j) > 37.325;
                end
            end
            Pfail(i,k) = ...
                norm(unscrunch(dPhat{k}(:,i))) > ...
                0.1*norm(unscrunch(Phat_test{k}(:,i))); 
            Pstarfail(i,k) = ...
                norm(unscrunch(dPstar{k}(:,i))) > ...
                0.1*norm(unscrunch(Pstar_test{k}(:,i)));
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
        any(any(any(estarfail))), ...
        any(any(Pfail)), ...
        any(any(Pstarfail))]);
end % function

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
