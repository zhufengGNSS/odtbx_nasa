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

% Begin function.
function failed = estspf_regression

failed = 0;
tol = 5e-9;

% Run the self tests first
disp('Running ESTSPF selftests...')

failed=estspf;
if failed
    warning('estspf_regression: Selftests failure.');
end

disp('Running ESTSPF regression test case 1...')

% Now check the nonlinear example of the planar 2-body problem

load estspf_reg_data.mat

% Set up initial state
X0=[6578 0 0 sqrt(3.986e5/6578)]';

% Set up initial covariance matrix
P0=eye(4);P0(1,1)=100;P0(2,2)=100;P0(3,3)=.01;P0(4,4)=.01;
% Argument in propagation function
MU=3.986e5;
% Time span
tspan=0:10:180;

ns=length(X0);
ncases=2;
options=setOdtbxOptions('MonteCarloCases',ncases);
options=setOdtbxOptions(options,'EditFlag',2);  % Accept all measurements

randn('state',1);
[~,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW]=...
                        estspf(@pr2bpd,@range2Dd,tspan,X0,P0,options,MU,.01);
                    
if any(any(any(abs(cell2mat(P) - cell2mat(Ps)) > tol)))
    warning('estspf_regression: Phat mismatch.');
    failed = 1;
end
if any(any(any(abs(PA - PAs) > tol)))
    warning('estspf_regression: PA mismatch.');
    failed = 1;
end
if  any(any(any(abs(PV - PVs) > tol)))
    warning('estspf_regression: PV mismatch.');
    failed = 1;
end
if any(any(any(abs(PW - PWs) > tol)))
    warning('estspf_regression: PW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATW - PHATWs) > tol)))
    warning('estspf_regression: PHATW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATV - PHATVs) > tol)))
    warning('estspf_regression: PHATV mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATA - PHATAs) > tol)))
    warning('estspf_regression: PHATA mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(X) - cell2mat(Xs)) > tol)))
    warning('estspf_regression: X mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(E) - cell2mat(Es)) > tol)))
    warning('estspf_regression: E mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(DY) - cell2mat(DYs)) > tol)))
    warning('estspf_regression: DY mismatch.');
    failed = 1;
end

% Rerun test with measurement editing
options=setOdtbxOptions(options,'EditFlag',1);
options=setOdtbxOptions(options,'EditRatio',3);

disp('Running ESTSPF regression test case 2 - measurement editing...')

[~,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW,EFLAG,PDY]=...
                        estspf(@pr2bpd,@range2Dd,tspan,X0,P0,options,MU,.01);
disp('Completed')
                    
test2 = load('estspf_reg_data2.mat');

if any(any(any(abs(cell2mat(P) - cell2mat(test2.P)) > tol)))
    warning('estspf_regression: Phat mismatch.');
    failed = 1;
end
if any(any(any(abs(PA - test2.PA) > tol)))
    warning('estspf_regression: PA mismatch.');
    failed = 1;
end
if  any(any(any(abs(PV - test2.PV) > tol)))
    warning('estspf_regression: PV mismatch.');
    failed = 1;
end
if any(any(any(abs(PW - test2.PW) > tol)))
    warning('estspf_regression: PW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATW - test2.PHATW) > tol)))
    warning('estspf_regression: PHATW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATV - test2.PHATV) > tol)))
    warning('estspf_regression: PHATV mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATA - test2.PHATA) > tol)))
    warning('estspf_regression: PHATA mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(X) - cell2mat(test2.X)) > tol)))
    warning('estspf_regression: X mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(E) - cell2mat(test2.E)) > tol)))
    warning('estspf_regression: E mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(DY) - cell2mat(test2.DY)) > tol)))
    warning('estspf_regression: DY mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(EFLAG) - cell2mat(test2.EFLAG)) > tol)))
    warning('estspf_regression: DY mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(PDY) - cell2mat(test2.PDY)) > tol)))
    warning('estspf_regression: DY mismatch.');
    failed = 1;
end



% End of function.

                                      
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

% % Propagation function for case 2: true
% function [Xnext,rtQd] = irwbprop(tspan,Xcurr,q)
% % dt = tspan(2)-tspan(1);
% % n = size(Xcurr,1);
% % O = zeros(n,n);
% % A([1:3 7:8],[4:6 8]) = eye(5,4);
% % Xnext = expm(A*dt)*Xcurr;
% % Q([4:6 8],[4:6 8]) = q*eye(4);
% % Psi = expm([-A Q; O A']*dt);
% % Phi = Psi(n+1:2*n,n+1:2*n)';
% % rtQd = chol(Phi*Psi(1:n,n+1:2*n))';
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
% 
% % Propagation function for case 2: estimate, solve for only
% function [Xnext,rtQd] = irwprop(tspan,Xcurr,q)
% % dt = tspan(2)-tspan(1);
% % n = size(Xcurr,1);
% % O = zeros(n,n);
% % A(:,4:6) = eye(6,3);
% % Xnext = expm(A*dt)*Xcurr;
% % Q(4:6,4:6) = q*eye(3);
% % Psi = expm([-A Q; O A']*dt);
% % Phi = Psi(n+1:2*n,n+1:2*n)';
% % rtQd = chol(Phi*Psi(1:n,n+1:2*n))';
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
% 
% % Measurement function for case 1: true and estimate
% function [Y,rtR] = rwmeas(t,X,sig) %#ok<INUSL>
% Y = X;
% rtR = sig;
% 
% % Measurement function for case 2: true
% function [Y,rtR] = irwbmeas(t,X,sig) %#ok<INUSL>
% % H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
% % Y = H*X;
% % rtR = sig*eye(3);
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
% 
% % Measurement function for case 2: estimate, solve for only
% function [Y,rtR] = irwmeas(t,X,sig) %#ok<INUSL>
% % H = [eye(3) zeros(3,3)];
% % Y = H*X;
% % rtR = sig*eye(3);
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
% 
% % Continuous function from estseq.m
% 
% % Dynamics function for case 2: true
% function [Xdot,A,Q] = irwbdyn(t,X,q)
% el = length(t);
% A([1:3 7:8],[4:6 8]) = eye(5,4);
% Xdot = A*X;
% A = repmat(A,[1 1 el]);
% Q([4:6 8],[4:6 8]) = q*eye(4);
% Q = repmat(Q,[1 1 el]);
% 
% % Dynamics function for case 2: estimate, solve for only
% function [Xdot,A,Q] = irwdyn(t,X,q)
% el = length(t);
% A(:,4:6) = eye(6,3);
% Xdot = A*X;
% A = repmat(A,[1 1 el]);
% Q(4:6,4:6) = q*eye(3);
% Q = repmat(Q,[1 1 el]);
% 
% % Measurement function for case 2: true
% function [Y,H,R] = irwbdat(t,X,r)
% el = length(t);
% H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
% Y = H*X;
% H = repmat(H,[1 1 el]);
% R = r*repmat(eye(3),[1 1 el]);
% 
% % Measurement function for case 2: estimate, solve for only
% function [Y,H,R] = irwdat(t,X,r)
% el = length(t);
% H = [eye(3) zeros(3,3)];
% Y = H*X;
% H = repmat(H,[1 1 el]);
% R = r*repmat(eye(3),[1 1 el]);
