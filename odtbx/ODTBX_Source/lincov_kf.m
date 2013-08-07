function varargout = lincov_kf(varargin)
% LINCOV_KF  Linear covariance analysis for Kalman Filter.
% Perform a general covariance analysis linearized about the reference.
% Assume that design values and true values may differ for the initial
% covariance (Pa), measurement noise covariance (Pv), and process noise 
% covariance (Pw).  Assume that the dynamics and measurement partials, and 
% the process and measurement noise covariances, may be time-varying. 
% Assume that linear, possibly time-varying, transformations partition the 
% state space into a "solve-for" subspace which is to be estimated from the 
% measurements, and a "consider" subspace which will not be estimated. 
% Compute the contributions due to _a priori_ uncertainty, measurement 
% noise, and process noise.  Compute the sensitivity to _a priori_ errors, 
% both solve-for and consider.
%
%  This function is called by ESTSEQ and ESTSRIF.
%
%   [P,PA,PV,PW,PM,PHATA,PHATV,PHATW,PHATM,SIG_A,PDYT] = LINCOV_KF(...
%   TSPAN,TINT,TITER,NITER,S,C,PO,PBARO,XREF,HREF,YREF,R,XSREF,HSREF,...
%   YBAR,RHAT,DYNFUN,DYNARG,DEMOMODE,ISCHMIDT,OPTIONS)
%     The inputs and outputs are all specific to the estimators estseq and
%     estsrif. See those functions for more information.
%
%   [...] = LINCOV_KF(...,PAO,PVO,PWO,PMO,PHATAO,PHATVO,PHATWO,PHATMO,SIG_AO)
%     Includes the "restart record" information as needed by the estimators
%     estseq and estsrif.
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


tspan=varargin{1};
tint=varargin{2};
titer=varargin{3};
niter=varargin{4};
S=varargin{5};
C=varargin{6};
Po=varargin{7};
Pbaro=varargin{8};
Xref=varargin{9};
Href=varargin{10};
Yref=varargin{11};
R=varargin{12};
Xsref=varargin{13};
Hsref=varargin{14};
Ybar=varargin{15};
Rhat=varargin{16};
dynfun=varargin{17};
dynarg=varargin{18};
demomode=varargin{19};
ischmidt=varargin{20};

% The input ODTBX options structure can be specified singly, or separately
% for truth and estimated states.
if all(isfield(varargin{21}, {'tru','est'})),
    options = varargin{21};
else
    options.tru = varargin{21};
    options.est = options.tru;
end


if nargin > 21
    Pao = varargin{22};
    Pvo = varargin{23};
    Pwo = varargin{24};
    Pmo = varargin{25};
    Phatao = varargin{26};
    Phatvo = varargin{27};
    Phatwo = varargin{28};
    Phatmo = varargin{29};
    Sig_ao = varargin{30};
    restart = 1;
else
    restart = 0;
end

%%
% *Solve-For and Consider Mapping*
%
% The mapping of the state-space into solve-for and consider subspaces is
% defined according to
%
% $$ s(t) = S(t) x(t), \quad c(t) = C(t) x(t) $$
%
% $$ M(t) = \Bigl[ S(t);\, C(t) \Bigr], \quad
% M^{-1}(t) = \Bigl[ \tilde{S}(t),\, \tilde{C}(t) \Bigr]$$
%
% $$ x(t) = \tilde{S}(t) s(t) + \tilde{C}(t) c(t) $$

lentr = length(titer);
m = size(Yref,1);
lents = length(tspan);
% Indices within tint that point back to tspan, i.e., tint(ispan)=tspan
[~,ispan] = ismember(tspan,tint,'legacy');
% Find indices within titer that point back to tint
[~,iint] = ismember(tint,titer,'legacy');

ns = size(S(:,:,1),1);
nc = size(C(:,:,1),1);
n = ns + nc;
Stilde=zeros(n,ns,lentr);
Ctilde=zeros(n,nc,lentr);
for i = lentr:-1:1,
    Minv = inv([S(:,:,i);C(:,:,i)]);
    Stilde(:,:,i) = Minv(:,1:ns); 
    Ctilde(:,:,i) = Minv(:,ns+1:n); 
end

%%
% *Sequential Kalman Filter Measurement Update*
%
% $$\quad \hat{P}(t_i^-) = \hat{P}_a(t_i^-)+\hat{P}_v(t_i^-)+\hat{P}_w(t_i^-)$$
%
% $$ K_i = \hat{P}(t_i^-) H_s(t_i)'
% \Bigl( H_s(t_i) \hat{P}(t_i^-) H_s(t_i)' + \hat{R}(t_i) \Bigr)^{-1} $$
%
% $$ \hat{P}_a(t_i^+) = \Bigl(I - K_i H_s(t_i)\Bigr) \hat{P}_a(t_i^-)
% \Bigl(I - K_i H_s(t_i)\Bigr)', \quad \hat{P}_a(t_1^-) = \hat{P}_o^- $$
%
% $$ P_a(t_i^+) = \Bigl(I - \tilde{S}(t_i) K_i H(t_i)\Bigr) P_a(t_i^-)
% \Bigl(I - \tilde{S}(t_i) K_i H(t_i)\Bigr)', \quad P_a(t_1^-) = P_o^- $$
%
% $$ \hat{P}_v(t_i^+) = \Bigl(I - K_i H_s(t_i)\Bigr) \hat{P}_v(t_i^-)
% \Bigl(I - K_i H_s(t_i)\Bigr)' + K_i \hat{R}(t_i) K_i', \quad \hat{P}_v(t_1^-) = 0 $$
%
% $$\begin{array}{l} 
% P_v(t_i^+)=\Bigl(I - \tilde{S}(t_i) K_i H(t_i)\Bigr) P_v(t_i^-)
%                  \Bigl(I - \tilde{S}(t_i)K_i H(t_i)\Bigr)'
%             + \tilde{S}(t_i)K_i R(t_i) K_i'\tilde{S}'(t_i),\\ 
% \quad P_v(t_1^-) = 0 \end{array}$$ 
%
% $$ \hat{P}_w(t_i^+) = \Bigl(I - K_i H_s(t_i)\Bigr) \hat{P}_w(t_i^-)
% \Bigl(I - K_i H_s(t_i)\Bigr)',\quad \hat{P}_w(t_1^-) = 0 $$
%
% $$ P_w(t_i^+) = \Bigl(I - \tilde{S}(t_i) K_i H(t_i)\Bigr) P_w(t_i^-)
% \Bigl(I - \tilde{S}(t_i) K_i H(t_i)\Bigr)',\quad P_w(t_1^-) = 0 $$
%
% *Sequential Kalman Filter Time Update*
%
% $$ \hat{P}_a(t_i^-) = \Phi_{ss}(t_i,t_{i-1}) \hat{P}_a(t_{i-1}^+)
% \Phi_{ss}'(t_i,t_{i-1})$$
%
% $$ P_a(t_i^-) = \Phi(t_i,t_{i-1}) P_a(t_{i-1}^+) \Phi'(t_i,t_{i-1})$$
%
% $$ \hat{P}_v(t_i^-) = \Phi_{ss}(t_i,t_{i-1}) \hat{P}_v(t_{i-1}^+)
% \Phi_{ss}'(t_i,t_{i-1})$$
%
% $$ P_v(t_i^-) = \Phi(t_i,t_{i-1}) P_v(t_{i-1}^+) \Phi'(t_i,t_{i-1})$$
%
% $$ \hat{P}_w(t_i^-) = \Phi_{ss}(t_i,t_{i-1}) \hat{P}_w(t_{i-1}^+)
% \Phi_{ss}'(t_i,t_{i-1}) + \hat{Q}_d(t_i,t_{i-1})$$
%
% $$ P_w(t_i^-) = \Phi(t_i,t_{i-1}) P_w(t_{i-1}^+) \Phi'(t_i,t_{i-1})
% + Q_d(t_i,t_{i-1})$$
%
% *Sensitivity Matrix Time and Measurement Update*
%
% The sensitivity matrix shows the linear sensitivity of the solution to
% mismodeling of the distributions of the _a priori_ parameters (both
% solve-fors and considers):
%
% $$ \Sigma_a(t_i) = \Bigl(I - \tilde{S}(t_i) K_i
% H(t_i) \Bigr) \Phi(t_i,t_{i-1}) \Sigma_a(t_{i-1}), \quad
% \Sigma_a(t_o) = \Bigl[\tilde{S}(t_o), \tilde{C}(t_o) \Bigr] $$
%
% In the particular case for which R = r*I and Rhat = rhat*I, then
%
% $$\Delta P_v = S P_v S' - \hat{P}_v $$
%
% will contain terms of the form K K' (r - rhat). In this case, if one
% chooses (r-rhat) = 1, delta P_v will represent the sensitivity to
% measurement noise mismodeling.  Similarly, when Q = q*I, and Qhat =
% qhat*I, (q-qhat) will factor out of the process noise partition, so that
% if one chooses q-qhat = 1, delta P_w will represent the sensitivity to
% process noise mismodeling.

% Note that it is possible to compute the sensitivities to each particular
% measurement and process noise sample, but it is somewhat complicated for
% the sequential estimator, and in any case, does not appear to be useful.

% Pre-allocate arrays for consider covariance analysis
% True (total) covariance
Pa    = NaN([n,n,lentr]);       % a-priori
Pv    = NaN(size(Pa));          % measurement noise
Pw    = NaN(size(Pa));          % process noise
Pm    = NaN(size(Pa));          % maneuver noise
P     = NaN(size(Pa));          % total
Pdyt  = NaN(m,m,lentr);         % measurement innovations
% Assumed covariance of solved for states only if ischmidt=0 otherwise
% it is of all the states
Phata = NaN([size(Pbaro),lentr]);
Phatv = NaN(size(Phata));
Phatw = NaN(size(Phata));
Phatm = NaN(size(Phata));
Phat = NaN(size(Phata));
% Covariance error of solved for states only
dPa   = NaN([ns,ns,lentr]);
dPv   = NaN(size(dPa));
dPw   = NaN(size(dPa));
dPm   = NaN(size(dPa));
% Sensitivity matrix of all states to a-priori
Sig_a = NaN(size(Pa));

% Loop forward over the data span, performing the sequential covariance
% analysis.  Since this is based on a linearization about a fixed reference,
% iterations of the measurement update will have no effect here; however,
% we want everything to synch up with any possible iterations of the EKF in
% the simululation below, so we'll just copy the update niter times.
for i = 1:lents,

    % ---------------------------------
    % Time update
    % ---------------------------------

    if i == 1, % Assign a-priori data only to a priori partitions
        thisint = 1;

        if restart
            % True covariance
            Pa(:,:,1) = Pao;
            Pv(:,:,1) = Pvo;
            Pw(:,:,1) = Pwo;
            Pm(:,:,1) = Pmo;
            P(:,:,1) = Pao+Pvo+Pwo+Pmo;
            % Assumed covariance
            Phata(:,:,1) = Phatao;
            Phatv(:,:,1) = Phatvo;
            Phatw(:,:,1) = Phatwo;
            Phatm(:,:,1) = Phatmo;
            Phat(:,:,1) = Phatao+Phatvo+Phatwo+Phatmo;
            Sig_a(:,:,1) = Sig_ao;
        else
            % True covariance
            Pa(:,:,1) = Po;
            Pv(:,:,1) = zeros(size(Po));
            Pw(:,:,1) = zeros(size(Po));
            Pm(:,:,1) = zeros(size(Po));
            P(:,:,1) = Pa(:,:,1);
            % Assumed covariance
            Phata(:,:,1) = Pbaro;
            Phatv(:,:,1) = zeros(size(Pbaro));
            Phatw(:,:,1) = zeros(size(Pbaro));
            Phatm(:,:,1) = zeros(size(Pbaro));
            Phat(:,:,1) = Phata(:,:,1);
            Sig_a(:,:,1) = [Stilde(:,:,1), Ctilde(:,:,1)];
        end
    else % Assign process noise only to process noise partitions

        thisint = iint(ispan(i-1)):iint(ispan(i))-niter;

        % True covariance
        [~,~,Phi,sdum] = integ(dynfun.tru,titer(thisint),Xref(:,ispan(i-1)),options.tru,dynarg.tru);
        if length(thisint) == 2 % This is because for time vector of length 2, ode outputs >2
            Phi(:,:,2) = Phi(:,:,end);
            Phi = Phi(:,:,1:2);
            sdum(:,:,2) = sdum(:,:,end);
            sdum = sdum(:,:,1:2);
        end
        for k = 2:length(thisint)
            sdum(:,:,k) = (sdum(:,:,k) + sdum(:,:,k)')/2;
            Pw(:,:,thisint(k)) = Phi(:,:,k)*Pw(:,:,thisint(1))*Phi(:,:,k)' + sdum(:,:,k);
            Pw(:,:,thisint(k)) = (Pw(:,:,thisint(k)) + Pw(:,:,thisint(k))')/2;
            Pa(:,:,thisint(k)) = Phi(:,:,k)*Pa(:,:,thisint(1))*Phi(:,:,k)';
            Pv(:,:,thisint(k)) = Phi(:,:,k)*Pv(:,:,thisint(1))*Phi(:,:,k)';
            Pm(:,:,thisint(k)) = Phi(:,:,k)*Pm(:,:,thisint(1))*Phi(:,:,k)';
            Sig_a(:,:,thisint(k)) = Phi(:,:,k)*Sig_a(:,:,thisint(1));
        end

        % Assumed covariance
        % If ischmidt=1, the Phiss corresponds to the full state transition
        % matrix, and not just the solved-for states.  Similarly for Qdhat.
         [~,~,Phiss,sdum] = integ(dynfun.est,titer(thisint),Xsref(:,ispan(i-1)),options.est,dynarg.est);
         if length(thisint) == 2 % This is because for time vector of length 2, ode outputs >2
             Phiss(:,:,2) = Phiss(:,:,end);
             Phiss = Phiss(:,:,1:2);
             sdum(:,:,2) = sdum(:,:,end);
             sdum = sdum(:,:,1:2);
         end
         for k = 2:length(thisint)
             sdum(:,:,k) = (sdum(:,:,k) + sdum(:,:,k)')/2;
             Phatw(:,:,thisint(k)) = Phiss(:,:,k)*Phatw(:,:,thisint(1))*Phiss(:,:,k)' + sdum(:,:,k);
             Phatw(:,:,thisint(k)) = (Phatw(:,:,thisint(k)) + Phatw(:,:,thisint(k))')/2;
             Phata(:,:,thisint(k)) = ...
                 Phiss(:,:,k)*Phata(:,:,thisint(1))*Phiss(:,:,k)';
             Phatv(:,:,thisint(k)) = ...
                 Phiss(:,:,k)*Phatv(:,:,thisint(1))*Phiss(:,:,k)';
             Phatm(:,:,thisint(k)) = ...
                 Phiss(:,:,k)*Phatm(:,:,thisint(1))*Phiss(:,:,k)';
         end

    end

    % Covariance differences and formal covariance over prop interval
    for j = thisint,
        if ischmidt == 1 % Both the true and assumed are the same size
            dPa(:,:,j) = S(:,:,j)*(Pa(:,:,j) - Phata(:,:,j))*S(:,:,j)';
            dPv(:,:,j) = S(:,:,j)*(Pv(:,:,j) - Phatv(:,:,j))*S(:,:,j)';
            dPw(:,:,j) = S(:,:,j)*(Pw(:,:,j) - Phatw(:,:,j))*S(:,:,j)';
            dPm(:,:,j) = S(:,:,j)*(Pm(:,:,j) - Phatm(:,:,j))*S(:,:,j)';
        else
            dPa(:,:,j) = S(:,:,j)*Pa(:,:,j)*S(:,:,j)' - Phata(:,:,j);
            dPv(:,:,j) = S(:,:,j)*Pv(:,:,j)*S(:,:,j)' - Phatv(:,:,j);
            dPw(:,:,j) = S(:,:,j)*Pw(:,:,j)*S(:,:,j)' - Phatw(:,:,j);
            dPm(:,:,j) = S(:,:,j)*Pm(:,:,j)*S(:,:,j)' - Phatm(:,:,j);
        end
        Phat(:,:,j) = Phata(:,:,j) + Phatv(:,:,j) + Phatw(:,:,j)...
            + Phatm(:,:,j);
        P(:,:,j) = Pa(:,:,j) + Pv(:,:,j) + Pw(:,:,j) + Pm(:,:,j);
    end

    % ---------------------------------
    % Kalman filter measurement update
    % ---------------------------------

    k = thisint(end)+1;

    inan = isnan(Yref(:,i)) | isnan(Ybar(:,i));
    Hsrefnn = Hsref(~inan,:,i);
    Hrefnn = Href(~inan,:,i);
    Rhatnn = Rhat(~inan,~inan,i);
    Rnn = R(~inan,~inan,i);
    
    Pdyt(:,:,k-1) = NaN(size(R(:,:,i)));
    Pdytnn = (Hrefnn*P(:,:,k-1)*Hrefnn' + Rnn);
    
    Pdyt(~inan,~inan,k-1) = Pdytnn;
    
    % Compute the gains
    K = Phat(:,:,k-1)*Hsrefnn'/...
        (Hsrefnn*Phat(:,:,k-1)*Hsrefnn' + Rhatnn);

    if ischmidt == 1

        % We apply the gains only to the solve-for
        K = S(:,:,i)*K;
                    
        % Calling kalmup may have extra overhead but could allow it for
        % iterative Kalman filter - need to flesh this out
        %         [xtmp,Ptmp,efltmp,dytmp,Pdytemp,K] = kalmup(datfun.est,...
        %             tspan(i),Xsref(:,i),Phat(:,:,k-1),Ybar(:,i),options,[],...
        %             [],datarg.est,[],S(:,:,i),C(:,:,i))

        % Update the total asssumed covariance
        ImSKH = eye(ns+nc) - Stilde(:,:,k-1)*K*Hsrefnn;
        Phata(:,:,k) = ImSKH*Phata(:,:,k-1)*ImSKH';
        Phatv(:,:,k) = ImSKH*Phatv(:,:,k-1)*ImSKH' + ...
            Stilde(:,:,k-1)*K*Rhatnn*K'*Stilde(:,:,k-1)';
        Phatw(:,:,k) = ImSKH*Phatw(:,:,k-1)*ImSKH';
        Phatm(:,:,k) = ImSKH*Phatm(:,:,k-1)*ImSKH';

    else

        % Update the asssumed covariance
        ImKH = eye(ns) - K*Hsrefnn;
        Phata(:,:,k) = ImKH*Phata(:,:,k-1)*ImKH';
        Phatv(:,:,k) = ImKH*Phatv(:,:,k-1)*ImKH' + K*Rhatnn*K';
        Phatw(:,:,k) = ImKH*Phatw(:,:,k-1)*ImKH';
        Phatm(:,:,k) = ImKH*Phatm(:,:,k-1)*ImKH';

    end

    % Update the true covariance. Assign measurement noise only to
    % measurement noise partitions of the total covariance
    ImSKH = eye(ns+nc) - Stilde(:,:,k-1)*K*Hrefnn;
    Pa(:,:,k) = ImSKH*Pa(:,:,k-1)*ImSKH';
    Pv(:,:,k) = ImSKH*Pv(:,:,k-1)*ImSKH' ...
        + Stilde(:,:,k-1)*K*Rnn*K'*Stilde(:,:,k-1)';
    Pw(:,:,k) = ImSKH*Pw(:,:,k-1)*ImSKH';
    Pm(:,:,k) = ImSKH*Pm(:,:,k-1)*ImSKH';
    P(:,:,k) = Pa(:,:,k) + Pv(:,:,k) + Pw(:,:,k) + Pm(:,:,k);

    % Update the sensitivity matrix to apriori
    Sig_a(:,:,k) = ImSKH*Sig_a(:,:,k-1);

    if ischmidt == 1

        % Post-update covariance differences and formal covariance
        dPa(:,:,k) = S(:,:,j)*(Pa(:,:,k) - Phata(:,:,k))*S(:,:,j)';
        dPv(:,:,k) = S(:,:,j)*(Pv(:,:,k) - Phatv(:,:,k))*S(:,:,j)';
        dPw(:,:,k) = S(:,:,j)*(Pw(:,:,k) - Phatw(:,:,k))*S(:,:,j)';
        dPm(:,:,k) = S(:,:,j)*(Pm(:,:,k) - Phatm(:,:,k))*S(:,:,j)';

    else

        % Post-update covariance differences and formal covariance
        dPa(:,:,k) = S(:,:,k)*Pa(:,:,k)*S(:,:,k)' - Phata(:,:,k);
        dPv(:,:,k) = S(:,:,k)*Pv(:,:,k)*S(:,:,k)' - Phatv(:,:,k);
        dPw(:,:,k) = S(:,:,k)*Pw(:,:,k)*S(:,:,k)' - Phatw(:,:,k);
        dPm(:,:,k) = S(:,:,k)*Pm(:,:,k)*S(:,:,k)' - Phatm(:,:,k);

    end

    Phat(:,:,k) = Phata(:,:,k) + Phatv(:,:,k) + Phatw(:,:,k)...
        + Phatm(:,:,k);

    % Now copy the data if necessary to fill in for any iterations below.
    % Zero-order hold
    for j = (k+1):iint(ispan(i)),
        Pa(:,:,j) = Pa(:,:,k);
        Pv(:,:,j) = Pv(:,:,k);
        Pw(:,:,j) = Pw(:,:,k);
        Pm(:,:,j) = Pm(:,:,k);
        Phata(:,:,j) = Phata(:,:,k);
        Phatv(:,:,j) = Phatv(:,:,k);
        Phatw(:,:,j) = Phatw(:,:,k);
        Phatm(:,:,j) = Phatm(:,:,k);
        dPa(:,:,j) = dPa(:,:,k);
        dPv(:,:,j) = dPv(:,:,k);
        dPw(:,:,j) = dPw(:,:,k);
        dPm(:,:,j) = dPm(:,:,k);
        Phat(:,:,j) = Phat(:,:,k);
        P(:,:,j) = P(:,:,k);
        Sig_a(:,:,j) = Sig_a(:,:,k);
    end

end

%% Variance Sandpiles
% Generate "variance sandpiles," which are stacked area charts showing the
% the time series of each solve-for variance's contribution from _a
% priori_ error variance, measurement noise variance, and process noise
% variance.  This can be done is several ways.  The true variance and the
% formal variance are
%
% $$ P = P_a + P_v + P_w, \quad \hat{P} = \hat{P}_a + \hat{P}_v $$
%
% and the delta variances are
%
% $$ \Delta P_a = S P_a S' - \hat{P}_a, \quad  \Delta P_v = S P_v S' -
% \hat{P}_v, \quad \Delta P_w = S P_w S' $$
%
% When all the delta variances are positive, the sandpile should show the
% formal variance, and the deltas due to each component.  When all the
% deltas are negative, the sandpile should show the true variance, and
% negatives of the deltas.  Otherwise, plot the components of the true
% variance as a positive sandpile, and the components of the formal
% variance as a negative sandpile, and relabel the negative y-axes to
% indicate this.

% NOTE: For the demomode examples plotted below, the pre- and
% post-multiplication of P_a, P_v, and P_w by S and S',
% respectively, have been ignored for simplicity since the
% solve-for states are the first ns states.

if demomode,
    for i = 1:ns,
        figure(i)
        clf
        varpiles(titer,dPa(i,i,:),dPv(i,i,:),dPw(i,i,:),...
            Pa(i,i,:),Pv(i,i,:),Pw(i,i,:),...
            Phata(i,i,:),Phatv(i,i,:),Phatw(i,i,:),...
            P(i,i,:),Phat(i,i,:))
    end
end

%% Sensitivity Mosaics
% Generate "sensitivity mosaics," which are checkerboard plots of the
% sensitivity matrices.  For some reason, Matlab's pcolor function does not
% plot the final row and column, so append an extra row and column to get
% the correct plot.  Initially plot the sensitivity at the final time,
% and put up a slider that lets the user scan through sensitivities over
% the time span.

% First map the full sensitivities to the solve-for state space:
Sig_sa = NaN(ns,n,lentr);
for j = lentr:-1:1,
    Sig_sa(:,:,j) = S(:,:,j)*Sig_a(:,:,j); 
end
if demomode,
    % Now create plot and slider:
    lastfig = ns;
    figure(lastfig+1)
    hold off
    pcolor(eye(ns+1,ns)*Sig_sa(:,:,end)*eye(n,n+1))
    set(gca,'xtick',1:n,'ytick',1:ns)
    xlabel('{\it A Priori} State Index')
    sa = uicontrol(gcf,'style','slider',...
        'max',lentr,'min',1,...
        'value',lentr,...
        'sliderstep',[1/lentr,0.5],...
        'units','normalized','position',...
        get(gca,'position')*[1 0 0 0;0 1 0 -.1;0 0 1 0;0 0 0 .1]');
    ca = @(h,e) pcolor(eye(ns+1,ns)*...
        Sig_sa(:,:,round(get(sa,'value')))*eye(n,n+1));
    set(sa,'callback', ca)
    ylabel('Solve-For State Index')
    set(gca,'ydir','rev','xaxisloc','top')
    axis equal
    colorbar
    title('Sensitivity Mosaic')
    hold on
end

varargout{1}=P;
varargout{2}=Pa;
varargout{3}=Pv;
varargout{4}=Pw;
varargout{5}=Pm;
varargout{6}=Phata;
varargout{7}=Phatv;
varargout{8}=Phatw;
varargout{9}=Phatm;
varargout{10}=Sig_a;
varargout{11}=Pdyt;
