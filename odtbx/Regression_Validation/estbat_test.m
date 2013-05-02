function failed = estbat_test
%
% estbat_test Regression testing and demo for estbat.m
% See also estbat.m
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
%   Ravi Mathur             08/29/2012      Extract regression test from
%                                           original estbat.m
%   Ravi Mathur             05/01/2013      Fully extracted regression test
%                                           from original estbat.m

totaltests = 3;
disp('estbat_test: regression testing estbat...')

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
        fail(k+totaltests,:) = run_test(k);
    end
    matlabpool('close');
end

failed = any(any(fail));

end % function estbat_test

% Runs the actual estbat tests. The contents of this function have been
% extracted from estbat.
function fail = run_test(testnum)
    % Set functions based on the test being performed
    switch testnum
        case 1
            dynfun.tru = @rwdyn;
            dynfun.est = @rwdyn;
            datfun.tru = @rwdat;
            datfun.est = @rwdat;
            load estbat_test1;

        case 2
            dynfun.tru = @irwbdyn;
            dynfun.est = @irwdyn;          
            datfun.tru = @irwbdat;
            datfun.est = @irwdat;
            load estbat_test2;

        case 3
            dynfun.tru = @sogmbdyn;
            dynfun.est = @rwdyn;          
            datfun.tru = @sogmbdat;
            datfun.est = @rwdat;
            load estbat_test3;
    end
    
    % Call estbat with loaded data
    ncases = getOdtbxOptions(options.tru,'MonteCarloCases',1);
    lent = length(tspan);
    [~,~,Phat,e,~,Pa,Pv,Pw,Phata,Phatv,~,Sigma_a] = ...
        estbat(dynfun, datfun, tspan, Xo, Po, options, dynarg, datarg, S, C);
    
    % Compare estbat output to validated output
    for k = ncases:-1:1,
        % TODO: use chi2 test like in estseq instead of "9\sigma"
        dPhat{k} = Phat_test{k} - Phat{k}; %#ok<USENS>
        sPhat{k} = Phat_test{k} + Phat{k};
        de{k} = e_test{k} - e{k}; %#ok<USENS>
        % Is each error sample within 9\sigma of its corresponding test
        % value?  Is each Phat sample "close enough," in terms of the
        % matrix 2-norm (largest singular value) to its test value?
        for i = lent:-1:1,
            efail(i,k) = ...
                de{k}(:,i)'*inv(unscrunch(sPhat{k}(:,i)))*de{k}(:,i) ...
                > 81;  %#ok<MINV>
            Pfail(i,k) = ...
                norm(unscrunch(dPhat{k}(:,i))) > ...
                0.1*norm(unscrunch(Phat_test{k}(:,i)));
        end
    end
    P = Pa + Pv + Pw;
    fail = logical([...
        any(any(any(abs(P_test - P) > 1e-10))), ...
        any(any(any(abs(Pa_test - Pa) > 1e-10))), ...
        any(any(any(abs(Pv_test - Pv) > 1e-10))),...
        any(any(any(abs(Pw_test - Pw) > 1e-10))), ...
        any(any(any(abs(Phatv_test - Phatv) > 1e-10))), ...
        any(any(any(abs(Phata_test - Phata) > 1e-10))), ...
        any(any(any(abs(Sigma_a_test - Sigma_a) > 1e-10))), ...
        any(any(efail)), ...
        any(any(Pfail))]);
end

% Self-test user functions
function [Xdot,A,Q] = rwdyn(t,X,q) % Test 1 dynfun
el = length(t);
A = zeros(1,1,el);
Xdot = A*X;
Q = q*ones(1,1,el);
end

function [Y,H,R] = rwdat(t,X,r) % Test 1 datfun
Y = X;
H = ones(1,1,length(t));
R = r*ones(1,1,length(t));
end

function [Xdot,A,Q] = irwbdyn(t,X,q) % Test 2 dynfun.tru
el = length(t);
A([1:3 7:8],[4:6 8]) = eye(5,4);
Xdot = A*X;
A = repmat(A,[1 1 el]);
Q([4:6 8],[4:6 8]) = q*eye(4);
Q = repmat(Q,[1 1 el]);
end

function [Xdot,A,Q] = irwdyn(t,X,q) % Test 2 dynfun.est
el = length(t);
A(:,4:6) = eye(6,3);
Xdot = A*X;
A = repmat(A,[1 1 el]);
Q(4:6,4:6) = q*eye(3);
Q = repmat(Q,[1 1 el]);
end

function [Xdot,A,Q] = sogmbdyn(t,X,c) % Test 3 dynfun.tru
el = length(t);
A(:,2:3) = eye(3,2);
A(3,2:3) = [-c.w_n^2, -2*c.zeta*c.w_n];
Xdot = A*X;
A = repmat(A,[1 1 el]);
Q(3,3) = c.q;
Q = repmat(Q,[1 1 el]);
end

function [Y,H,R] = sogmbdat(t,X,r) % Test 3 datfun.tru
Y = X(1,:) + X(2,:);
H = [1 1 0];
H = repmat(H,[1,1,length(t)]);
R = r*ones(1,1,length(t));
end

function [Y,H,R] = irwbdat(t,X,r) % Test 2 datfun.tru
el = length(t);
H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
Y = H*X;
H = repmat(H,[1 1 el]);
R = r*repmat(eye(3),[1 1 el]);
end

function [Y,H,R] = irwdat(t,X,r) % Test 2 datfun.est
el = length(t);
H = [eye(3) zeros(3,3)];
Y = H*X;
H = repmat(H,[1 1 el]);
R = r*repmat(eye(3),[1 1 el]);
end