function failed = estspf_test
%
% estspf_test Regression testing and demo for estspf.m
% See also estspf.m
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

% Author: Sun Hur-Diaz
%         sun@emergentspace.com
%         Emergent Space Technologies, Inc.
%
%         May 2009
%
% Ravi Mathur         08/28/2012      Renamed to conform to new
%                                     regression testing framework
% Ravi Mathur         05/03/2013      Fully extracted regression test
%                                     from original estspf.m

totaltests = 4;
disp('estspf_test: regression testing estspf...')

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

end % function

% Runs the actual estspf tests. The contents of this function have been
% extracted from estspf.
function fail = run_test(testnum)
    % Set functions based on the test being performed
    switch testnum
        case 1
            %% Case 1: Needs documentation
            propfun.tru = @rwprop;
            propfun.est = @rwprop;
            measfun.tru = @rwmeas;
            measfun.est = @rwmeas;
            load estspf_test1;
            tol = 1e-10;

        case 2
            %% Case 2: Consider covariance
            propfun.tru = @irwbprop;
            propfun.est = @irwprop;
            measfun.tru = @irwbmeas;
            measfun.est = @irwmeas;
            load estspf_test2;
            tol = 1e-10;

        case 3
            %% Case 3: Nonlinear example of the planar 2-body problem
            propfun.tru = @pr2bpd;
            propfun.est = @pr2bpd;
            measfun.tru = @range2Dd;
            measfun.est = @range2Dd;
            load estspf_test3;
            randn('state',1);
            tol = 5e-9;

        case 4
            %% Case 4: Rerun test 3 with measurement editing
            % Note that case 3 MUST be run before case 4, due to the
            % explicit dependence on the random number generator state.
            propfun.tru = @pr2bpd;
            propfun.est = @pr2bpd;
            measfun.tru = @range2Dd;
            measfun.est = @range2Dd;
            load estspf_test4;
            tol = 5e-9;
    end
    
    %% Call estspf with loaded data
    [~,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW,EFLAG,PDY]=...
        estspf(propfun,measfun,tspan,Xo,Po,options,proparg,measarg,S,C);
    
    %% Check results against accepted values
    fail = logical([...
        any(any(any(abs(cell2mat(P_test) - cell2mat(P)) > tol))), ...
        any(any(any(abs(PA_test - PA) > tol))), ...
        any(any(any(abs(PV_test - PV) > tol))),...
        any(any(any(abs(PW_test - PW) > tol))), ...
        any(any(any(abs(PHATW_test - PHATW) > tol))), ...
        any(any(any(abs(PHATV_test - PHATV) > tol))), ...
        any(any(any(abs(PHATA_test - PHATA) > tol))), ...
        any(any(any(abs(cell2mat(X_test) - cell2mat(X)) > tol))), ...
        any(any(any(abs(cell2mat(E_test) - cell2mat(E)) > tol))), ...
        any(any(any(abs(cell2mat(DY_test) - cell2mat(DY)) > tol))) ...
        ]);
    
    % Test case 4 has extra checks
    if testnum == 4,
        fail = logical([fail, ...
            any(any(any(abs(cell2mat(EFLAG_test) - cell2mat(EFLAG)) > tol))), ...
            any(any(any(abs(cell2mat(PDY_test) - cell2mat(PDY)) > tol))) ...
            ]);
    else
        fail = [fail, 0, 0]; % Add dummy values to match array size with test 4
    end
end % function
                                      
function [Xnext,rtQd] = pr2bpd(tspan,Xcurr,mu)
% PR2BPD Discrete version of PR2BP: planar restricted two-body problem.
%   [Xnext,rtQd] = pr2bpd(tspan,Xcurr,mu)

[lx,nx]=size(Xcurr);
Xnext=zeros(lx,nx);
for i=1:nx
    if i==1
        [t,x,Phi,S] = integ(@pr2bp,tspan,Xcurr(:,i),[],mu);
        lt=length(t);
        Xnext(:,i)=x(:,lt);
        rtQd=chol(S(:,:,lt))';
    else
        [t,x] = integ(@pr2bp,tspan,Xcurr(:,i),[],mu);
        lt=length(t);
        Xnext(:,i)=x(:,lt);
    end
end

end % function

function [y,rtR] = range2Dd(t,x,sig)
% RANGE2DD Discrete version of range2D: Range measurement model for 2-D
% inertial state space.
%   [y,rtR] = range2Dd(t,x,sig)

[lx,nx]=size(x);
y=zeros(1,nx);
for i=1:nx
    if i==1
        [y(:,i),H,R] = range2D(t,x(:,i),sig);
        rtR=chol(R(:,:,1))';
    else
        y(:,i) = range2D(t,x(:,i),sig);
    end
end

end % function

%% Functions used by the tests

% Propagation function for case 1: true and estimate
function [Xnext,rtQd] = rwprop(tspan,Xcurr,q)
dt = tspan(2)-tspan(1);
Xnext = Xcurr;
rtQd = sqrt(q*dt);
end % function

% Propagation function for case 2: true
function [Xnext,rtQd] = irwbprop(tspan,Xcurr,q)
dt = tspan(2)-tspan(1);
n = size(Xcurr,1);
O = zeros(n,n);
A([1:3 7:8],[4:6 8]) = eye(5,4);
Xnext = expm(A*dt)*Xcurr;
Q([4:6 8],[4:6 8]) = q*eye(4);
Psi = expm([-A Q; O A']*dt);
Phi = Psi(n+1:2*n,n+1:2*n)';
rtQd = chol(Phi*Psi(1:n,n+1:2*n))';
% The following is the discrete formulation from the continuous model
% [lx,nx]=size(Xcurr);
% Xnext=zeros(lx,nx);
% for i=1:nx
%     if i==1
%         [t,x,Phi,S] = integ(@irwbdyn,tspan,Xcurr(:,i),[],q);
%         lt=length(t);
%         Xnext(:,i)=x(:,lt);
%         rtQd=chol(S(:,:,lt))';
%     else
%         [t,x] = integ(@irwbdyn,tspan,Xcurr(:,i),[],q);
%         lt=length(t);
%         Xnext(:,i)=x(:,lt);
%     end
% end
end % function

% Propagation function for case 2: estimate, solve for only
function [Xnext,rtQd] = irwprop(tspan,Xcurr,q)
dt = tspan(2)-tspan(1);
n = size(Xcurr,1);
O = zeros(n,n);
A(:,4:6) = eye(6,3);
Xnext = expm(A*dt)*Xcurr;
Q(4:6,4:6) = q*eye(3);
Psi = expm([-A Q; O A']*dt);
Phi = Psi(n+1:2*n,n+1:2*n)';
rtQd = chol(Phi*Psi(1:n,n+1:2*n))';
% The following is the discrete formulation from the continuous model
% [lx,nx]=size(Xcurr);
% Xnext=zeros(lx,nx);
% for i=1:nx
%     if i==1
%         [t,x,Phi,S] = integ(@irwdyn,tspan,Xcurr(:,i),[],q);
%         lt=length(t);
%         Xnext(:,i)=x(:,lt);
%         rtQd=chol(S(:,:,lt))';
%     else
%         [t,x] = integ(@irwdyn,tspan,Xcurr(:,i),[],q);
%         lt=length(t);
%         Xnext(:,i)=x(:,lt);
%     end
% end
end % function

% Measurement function for case 1: true and estimate
function [Y,rtR] = rwmeas(t,X,sig) %#ok<INUSL>
Y = X;
rtR = sig;
end % function

% Measurement function for case 2: true
function [Y,rtR] = irwbmeas(t,X,sig) %#ok<INUSL>
H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
Y = H*X;
rtR = sig*eye(3);
% The following is the discrete formulation from the continuous model
% [lx,nx]=size(X);
% Y=zeros(3,nx);
% for i=1:nx
%     if i==1
%         [Y(:,i),H,R] = irwbdat(t,X(:,i),sig^2);
%         rtR=chol(R(:,:,1))';
%     else
%         Y(:,i) = irwbdat(t,X(:,i),sig^2);
%     end
% end
end % function

% Measurement function for case 2: estimate, solve for only
function [Y,rtR] = irwmeas(t,X,sig) %#ok<INUSL>
H = [eye(3) zeros(3,3)];
Y = H*X;
rtR = sig*eye(3);
% The following is the discrete formulation from the continuous model
% [lx,nx]=size(X);
% Y=zeros(3,nx);
% for i=1:nx
%     if i==1
%         [Y(:,i),H,R] = irwdat(t,X(:,i),sig^2);
%         rtR=chol(R(:,:,1))';
%     else
%         Y(:,i) = irwdat(t,X(:,i),sig^2);
%     end
% end
end % function

%% Functions used by the discrete formulation from the continuous model

% % Continuous functions of Case 2 from estseq.m
% % ----------------------------------------------
% % Dynamics function for case 2: true
% function [Xdot,A,Q] = irwbdyn(t,X,q) %#ok<DEFNU>
% el = length(t);
% A([1:3 7:8],[4:6 8]) = eye(5,4);
% Xdot = A*X;
% A = repmat(A,[1 1 el]);
% Q([4:6 8],[4:6 8]) = q*eye(4);
% Q = repmat(Q,[1 1 el]);
% end % function
% 
% % Dynamics function for case 2: estimate, solve for only
% function [Xdot,A,Q] = irwdyn(t,X,q) %#ok<DEFNU>
% el = length(t);
% A(:,4:6) = eye(6,3);
% Xdot = A*X;
% A = repmat(A,[1 1 el]);
% Q(4:6,4:6) = q*eye(3);
% Q = repmat(Q,[1 1 el]);
% end % function
% 
% % Measurement function for case 2: true
% function [Y,H,R] = irwbdat(t,X,r) %#ok<DEFNU>
% el = length(t);
% H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
% Y = H*X;
% H = repmat(H,[1 1 el]);
% R = r*repmat(eye(3),[1 1 el]);
% end % function
% 
% % Measurement function for case 2: estimate, solve for only
% function [Y,H,R] = irwdat(t,X,r) %#ok<DEFNU>
% el = length(t);
% H = [eye(3) zeros(3,3)];
% Y = H*X;
% H = repmat(H,[1 1 el]);
% R = r*repmat(eye(3),[1 1 el]);
% end % function