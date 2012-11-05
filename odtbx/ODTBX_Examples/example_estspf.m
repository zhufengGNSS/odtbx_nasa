% This file is an example for the ODTBX Sigma Point Filter Estimator, estspf.m, with Monte Carlo runs.
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

% Begin function.
function example_estspf
X0=[6578 0 0 sqrt(3.986e5/6578)]';
P0=eye(4);P0(1,1)=100;P0(2,2)=100;P0(3,3)=.01;P0(4,4)=.01;
MU=3.986e5;
tspan=0:10:180;

ns=length(X0);
ncases=30;
options=setOdtbxOptions('MonteCarloCases',ncases);
options=setOdtbxOptions(options,'EditFlag',2);

randn('state',1);
[T,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW] = estseq(@pr2bp,@range2D,tspan,X0,P0,options,MU,.01);
randn('state',1);
[Ts,Xs,Ps,Es,DYs,PAs,PVs,PWs,PHATAs,PHATVs,PHATWs]=...
                                          estspf(@pr2bpd,@range2Dd,tspan,X0,P0,options,MU,.01);

ksig=3;
% Plot estseq data
lent = length(T{1});
ss = std(reshape([E{:}],ns,lent,ncases),0,3);
for i = 1:ns,
    figure(i)
    subplot(2,1,1)
    h=plot(T{1},ksig*[1;-1]*squeeze(sqrt(PHATA(i,i,:)+PHATV(i,i,:)+PHATW(i,i,:)))','g--',...
        T{1},ksig*[1;-1]*squeeze(sqrt(PA(i,i,:)+PV(i,i,:)+PW(i,i,:)))','c--');
    set(h,'linewidth',4);
    hold on
    h=plot(T{1},ksig*[1;-1]*ss(i,:),'r--');
    set(h,'linewidth',4);
    ch = get(gca,'children');
    legend(ch(1:2:5),['Empirical: \pm', num2str(ksig),'-\sigma'], ...
        ['"True:" \pm', num2str(ksig),'-\sigma'],...
        ['Formal: \pm', num2str(ksig),'-\sigma'])
    for j = 1:ncases,
        plot(T{j},E{j}(i,:))
    end
    hold off
    title(['ESTSEQ: ',num2str(ncases),'-case Monte Carlo: X(',num2str(i),')'])
    xlabel('t (sec)')
    ylabel(['Error and ',num2str(ksig),'-\sigma'])
end

% Plot estspf data
lent = length(Ts{1});
ss = std(reshape([Es{:}],ns,lent,ncases),0,3);
for i = 1:ns,
    figure(i)
    subplot(2,1,2)
    h=plot(Ts{1},ksig*[1;-1]*squeeze(sqrt(PHATAs(i,i,:)+PHATVs(i,i,:)+PHATWs(i,i,:)))','g--',...
        Ts{1},ksig*[1;-1]*squeeze(sqrt(PAs(i,i,:)+PVs(i,i,:)+PWs(i,i,:)))','c--');
    set(h,'linewidth',4);
    hold on
    h= plot(Ts{1},ksig*[1;-1]*ss(i,:),'r--');
    set(h,'linewidth',4);
    ch = get(gca,'children');
    legend(ch(1:2:5),['Empirical: \pm', num2str(ksig),'-\sigma'], ...
        ['"True:" \pm', num2str(ksig),'-\sigma'],...
        ['Formal: \pm', num2str(ksig),'-\sigma'])
    for j = 1:ncases,
        plot(Ts{j},Es{j}(i,:))
    end
    hold off
    title(['ESTSPF: ',num2str(ncases),'-case Monte Carlo: X(',num2str(i),')'])
    xlabel('t (sec)')
    ylabel(['Error and ',num2str(ksig),'-\sigma'])
end

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