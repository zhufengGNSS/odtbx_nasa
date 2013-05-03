function varargout = estpar(varargin)
% ESTPAR Particle Filter Estimator.
%
% ESTPAR is a Bayesian Estimator that is used to calculate the state 
% estimate based on "point-mass" or particle representations of probability
% densities. Particle Filters are sequential Monte Carlo methods that rely
% on random sampling to compute their results.
%
% [T,X,P] = ESTPAR(DYNFUN,DATFUN,TSPAN,X0,P0) with TSPAN = [T0 TFINAL] 
% integrates the system of differential equations x' = f(t,x) from T0 to TFINAL 
% for N particles seeded from the initial conditions X0 & P0. At TFINAL
% each particle is weighted based on the measurement (y) generated 
% by y = h(t,x) + v, and the measurement covariance (R). The weighted 
% average of the particles then provides X & P at TFINAL. Resampling 
% is utilized to eliminate diverging particles. To obtain updates at 
% multiple times T1, T2, ..., TFINAL, use TSPAN = [T0 T1 T2 ... TFINAL].
%
%   Function [F,A,Q]=DYNFUN(T,X,options) must return a column vector 
%   corresponding to f(t,x). If DYNFUN is "vectorized," then f(t,x) must be
%   a 2-D array with each column corresponding to f(t(i),x(t(i)). Function 
%   DYNFUN(T,X) must return an additional output if called with two output
%   arguments, which may either be a matrix corresponding to A(t), or else 
%   an empty matrix, in which case INTEG will numerically compute A(t) with
%   NUMJAC. If A(t) is supplied, it must be a 3-D array for the vectorized
%   case, with each "slice" corresponding to A(t(i)).  To include process
%   noise, DYNFUN must return the process noise spectral density matrix, Q
%   = E[ww'], where x' = f(t,x) + w, as an additional output.  If DYNFUN is
%   vectorized, then Q must be a 3-D array, with each "slice" corresponding
%   to Q(t(i)).
%
%   Function [h,H,R]=DATFUN(T,X,options) must return a column vector 
%   corresponding to h(t,x), and two additional outputs corresponding to 
%   the measurement partials, H(t) = dh(t,x)/dx, and the measurement noise 
%   covariance, R = E[vv'], where y = h(t,x) + v.  As an alternate to 
%   supplying H(t), DATFUN may return an empty matrix as its second output, 
%   in which case ESTINV will numerically compute H(t) using NUMJAC. If 
%   DATFUN is "vectorized," then h(t,x) must return as its first output a 
%   2-D array with each column corresponding to h(t(i),x(t(i)); its next 
%   two outputs must be 3-D arrays with each "slice" corresponding to 
%   H(t(i)) and R(t(i)), respectively.
%
%   The rows in the solution array X correspond to times returned in the 
%   column vector T, which are chosen by the integrator.  The rows in the 
%   solution covariance array P correspond to the unique lower triangular 
%   elements of P at the times T, appended by row from column 1 to the main 
%   diagonal.  The Ith row may be reformed into a matrix using 
%   UNSCRUNCH(P(I,:)).
%
%   [T,X,P] = ESTPAR(DYNFUN,DATFUN,TSPAN,X0,P0,OPTIONS) performs as above
%   with default properties replaced by values in OPTIONS, an argument
%   created with the SETODTBXOPTIONS function.  See ODTBXOPTIONS for
%   details. Commonly used options allow one to specify parameters or
%   features of the estimator, force model, and measurment model.  
%   Note that as of ODTBX R2013a, OPTIONS can be either a standard
%   ODTBXOPTIONS structure, or a struct of two ODTBXOPTIONS structures that
%   are used separately for truth and estimated computations. The latter
%   method is achieved by settings OPTIONS as a struct:
%      >> options.tru = setOdtbxOptions(...);
%      >> options.est = setOdtbxOptions(...);
%      >> [...] = estseq(..., options, ...);
%   Using this method, all options common to truth and estimated
%   computations are taken from the options.est structure.
%
%   [T,X,P] = ESTPAR(DYNFUN,DATFUN,TSPAN,X0,P0,OPTIONS,DYNARG,DATARG)
%   passes DYNARG to DYNFUN and DATARG to DATFUN as DYNFUN(T,X,DYNARG) and
%   DATFUN(T,X,DATARG), respectively.  Use OPTIONS = [] as a place holder
%   if no options are set.
%
%   [T,X,P] = ESTPAR(DYNFUN,DATFUN,TSPAN,X0,P0,OPTIONS,DYNARG,DATARG,S,C)
%   passes in solve-for and consider mapping matrices, S and C,
%   respectively.  These matrices partition the state into a solve-for
%   partition, S*x, and a consider partition, C*x.  Only parameters in the
%   former partition will be updated from the meaurements.  Use [] as a 
%   place holder for OPTIONS, DYNARG, and/or DATARG as necessary if these 
%   inputs are not required.  S and C can be 2-D or 3-D arrays.  If 3-D, 
%   the 3rd dimension corresponds to the time vector. However, time-
%   varying S and C are currently not implemented; therefore, constant S 
%   and C corresponding to the first time vector will be used.
%
%   To handle the different dimensions of the full state vs. the solve-
%   for and consider partitions, the user can either design DYNFUN and 
%   DATFUN to check for this, or specify DYNFUN and/or DATFUN as 
%   structures, whose fields are *.tru and *.est. The function specified 
%   in *.tru will be used to evaluate the full state, and the one in 
%   *.est will be used for the solve-for partition.  This is also a way 
%   to specify differences between the true and the estimator models of the
%   dynamics and the measurement data.  Similar conventions may be used 
%   for X0, P0, DYNARG, and DATARG, i.e. X0.Xo and X0.Xbaro, P0.Po and 
%   P0.Pbaro, DYNARG.tru and DYNARG.est, DATARG.tru and DATARG.est. 
%   INACTIVE ~ INCLUDED FOR FUTURE IMPLEMENTATION
%
%   [T,X,P,E] = ESTPAR(DYNFUN,DATFUN,TSPAN,X0,P0,...) also returns the
%   estimation errors, E.
%
%   [T,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW] = ESTPAR(...) returns the 
%   innovations, DY, and several addditional covariance matrices: PA, PV, 
%   and PW are the true covariances that arise only from the true  
%   _a priori_ covariance, the true measurement noise covariance, and the 
%   true process noise covariance, respectively;  PHATA, PHATV, PHATW are 
%   the estimator's covariances that arise only from the design values of 
%   the  _a priori_ covariance, the measurement noise covariance, and the 
%   process noise covariance.
%
%   [T,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW,EFLAG] = ESTPAR(...) also
%   returns EFLAG, the array of edit flag values for all cases for all 
%   measurement types for all times. The edit flags may have 
%   the following values:
%
%       0   = Measurement was rejected (based on the edit ratio settings)
%       1   = Measurement was checked and passed the edit ratio test
%       2   = Measurement was forced to be accepted (based on edit flag
%             settings)
%   INACTIVE ~ INCLUDED FOR FUTURE IMPLEMENTATION
%
%   [T,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW,EFLAG,PDY] = ESTPAR(...) 
%   also returns the formal covariance of the measurement 
%   innovations DY.
%   INACTIVE ~ INCLUDED FOR FUTURE IMPLEMENTATION
%
%   keyword: Estimation,
%
%   See also
%      options handling:      ODTBXOPTIONS, SETODTBXOPTIONS,
%                             GETODTBXOPTIONS
%      evaluating solutions:  ESTVAL
%      other ODEAS filters:   ESTBAT, ESTSEQ, ESTSRIF
%      other ODEAS utilties:  INTEG, OBSERV
%      ODE solvers:           ODE113, ODE23, ODE45
%      covariance storage:    SCRUNCH, UNSCRUNCH
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

% estpar.m developed for ODTBX by
% John Gaebler
% NASA Goddard Space Flight Center
%
% Adapted from particlefilter.m by
% Alinda Kenyana Mashiku
% NASA Goddard Space Flight Center

% NOTE: The interface design of this function is based on ODE45.M, which
% is copyrighted by The MathWorks, Inc.

% Modification History
% ---------------------
% 09/26/2011 Adapted existing code to ODTBX standards
% 04/02/2012 Cleaned up code and added HELP section
 %                      Commited to ODTBX Source Forge

if nargin >= 5, % ESTPAR(PROPFUN,MEASFUN,TSPAN,X0,P0...)
    if all(isfield(varargin{1}, {'tru','est'})),
        dynfun = varargin{1};
    else
        dynfun.tru = varargin{1};
        dynfun.est = varargin{1};
    end
    if all(isfield(varargin{2}, {'tru','est'})),
        datfun = varargin{2};
    else
        datfun.tru = varargin{2};
        datfun.est = varargin{2};
    end
    tspan = varargin{3};
    if isstruct(varargin{4}),
        Xo = varargin{4}.Xo;
        Xbaro = varargin{4}.Xbaro;
    else
        Xo = varargin{4};
        Xbaro = varargin{4};
    end
    if isstruct(varargin{5}),
        Po = varargin{5}.Po;
        Pbaro = varargin{5}.Pbaro;
    else
        Po = varargin{5};
        Pbaro = varargin{5};
    end
    if isempty(Po) || isempty(Pbaro)
        error('Initial covariance must be set!')
    end
elseif nargin ~= 0 && nargin ~= 1,
    error('There must be at least 5 inputs! (dynfun,datfun,tspan,X0,PO)')
end
if nargin >= 6,
    if all(isfield(varargin{6}, {'tru','est'})),
        options = varargin{6};
    else
        options.tru = varargin{6};
        options.est = options.tru;
    end
else
    options.tru = setOdtbxOptions('OdeSolvOpts',odeset);
    options.est = options.tru;
end
ncases = getOdtbxOptions(options.est,'MonteCarloCases',1); % Set to at least 1
N = getOdtbxOptions(options.est,'Particles',20);
% refint = getOdtbxOptions(options.tru,'refint',0);

if nargin >= 7,
    if all(isfield(varargin{7}, {'tru','est'}))
        dynarg = varargin{7};
    else
        dynarg.tru = varargin{7};
        dynarg.est = varargin{7};
    end
elseif nargin >= 5,
    dynarg.tru = [];
    dynarg.est = [];
end
if nargin >= 8,
    if all(isfield(varargin{8}, {'tru','est'}))
        datarg = varargin{8};
    else
        datarg.tru = varargin{8};
        datarg.est = varargin{8};
    end
elseif nargin >= 5,
    datarg.tru = [];
    datarg.est = [];
end
if nargin >= 9,
    if isa(varargin{9},'function_handle'),
        mapfun = varargin{9}; %#ok<NASGU> %TODO
    elseif isa(varargin{9},'numeric') % constant solve-for map
        S = varargin{9};
        C = []; %zeros(0,0,length(tspan)); % in case C is not input, solve for all states
    end
elseif nargin >= 5, % If S & C not input, solve for all states
    S = eye(size(Po));
    C = []; %zeros(0,0,length(tspan));
end
if nargin >= 10, % constant consider map
    C = varargin{10}; %repmat(varargin{10},[1,1,length(tspan)]);
end

%% PREP
% *Solve-For and Consider Mapping*
%
% The mapping of the state-space into solve-for and consider subspaces is
% defined according to
%
% $$ s(t) = S(t) x(t), \quad c(t) = C(t) x(t) $$
%
% $$ M(t) = \Bigl[ S(t);\, C(t) \Bigr], \quad
% M^{-1}(t) = \Bigl[ \tilde{S}(t)\, \tilde{C}(t) \Bigr]$$
%
% $$ x(t) = \tilde{S}(t) s(t) + \tilde{C}(t) c(t) $$% Determine the number of measurements

ytemp = feval(datfun.tru,tspan(1),Xo,datarg.tru);
nmeas = length(ytemp);

% Get montecarloseed and edit options
eopts = chkestopts(options.est,ncases,nmeas);

if(~isnan(eopts.monteseed(1)))
%     RandStream.setGlobalStream(eopts.monteseed(1));
    randn('state', eopts.monteseed(1));
end

Minv = inv([S;C]);
ns = size(S,1);               % Number of solve-for states
nx = size(Minv,1);            % Total number of states
% Stilde = Minv(:,1:ns);

%% Pre-allocation

tpp = reshape(repmat(tspan(:)',2,1),1,[]);  % double each time step for pre- and post- measurements
lent = length(tpp);
[~,ipp] = ismember(tpp,tspan);

X_est = [Xbaro,NaN(ns,lent-1)];              % Allocate solve-for states, same i.c. each case

% FOR NOW ONLY TRACK Pa AND P_est_a ; ALL OTHERS WILL REMAIN ZEROS
P_est_a = NaN([size(Pbaro),lent]);              % Allocate formal covariance due to a priori
[P_est_v,P_est_w] = deal(zeros(size(P_est_a))); % Allocate formal variance due to noise

Pa = NaN([size(Po),lent]);                      % Allocate true covariance of solve-for due to a priori
[Pv,Pw] = deal(zeros(size(Pa)));                % Allocate true variance of s due to noise

dy = NaN(nmeas,lent); 

%% Reference trajectory

[~,x_tru,~] = integ(dynfun.tru,tspan,Xo,options.tru,dynarg.tru);

[Y_tru,~,R_tru] = feval(datfun.tru,tspan,x_tru,datarg.tru);

%% Initialization

Xk = zeros(nx,N);
Wk_x = ones(1,N) * 1/N;

for i = 1:N
    Xk(:,i) = Xbaro + covsmpl(Pbaro); % Generation of n initial particles
end

X_est(:,1) = sum(Xk.*repmat(Wk_x,nx,1),2);
P_est_a(:,:,1) = Wk_x(ones(1,nx),:).*(Xk - X_est(:,ones(1,N)))*(Xk - X_est(:,ones(1,N)))';
        
n_thr = 0.25*N; % Threshold of choice for resampling
ReCount = 0;

%% Filtering
for i = 1:2:lent
    % i is pre-update
    % i+1 is post-update
    % i+2 is after propagation
        
    irnan = isnan(Y_tru(:,ipp(i)));
    Y_est = Y_tru(~irnan,ipp(i)) + covsmpl(R_tru(~irnan,~irnan,ipp(i)));

    Wk_y = zeros(1,N);
    
    %% Update
    if all(irnan) % if no tru measurement
        P_est_a(:,:,i+1) = P_est_a(:,:,i);
        X_est(:,i+1) = X_est(:,i);
    else
        for j = 1:N
            % Loop through N particles, find measurement esimate
            % for each particle and the measurement weight.
            [Yk,~,Rk] = feval(datfun.est,tpp(i),Xk(:,j),datarg.est);
            inan = isnan(Yk);
            Yk(inan)=0;
            P_y = Rk(~inan,~inan);
            GauDY = ((2*pi)^(size(Y_est,1))*det(P_y)).^(-0.5);
            Wk_y(j) = GauDY.*exp(-0.5*((Yk(~irnan)-Y_est)'/(P_y)*(Yk(~irnan)-Y_est)));
        end
        
        % Normalize weights
        Wk_y = Wk_y/sum(Wk_y);
        
        % Update the particle weights
        Wk_x = Wk_y .* Wk_x/(sum(Wk_y.*Wk_x));
        
        % Calculate the current X and P.
        X_est(:,i+1) = sum(Xk.*Wk_x(ones(1,nx),:),2);
        P_est_a(:,:,i+1) = (Xk - X_est(:,(i+1)*ones(1,N)))*(Xk - X_est(:,(i+1)*ones(1,N)))'/N;
       
        %% Resampling
        n_eff = 1/(sum(Wk_x.^2));
        if n_eff < n_thr
            [~,id_x] = sort(Wk_x); % sort in ascending
            Nrep = N - ceil(n_eff); % number of particles to replace
            id_x = id_x(1:Nrep); % index of lowest weights
            Xk(:,id_x) = X_est(:,(i+1)*ones(1,length(id_x))) + covsmpl(P_est_a(:,:,i+1),length(id_x)); % this assumes a gaussian resampling... 

            % easier based on equal weight distribution 
            Wk_x(id_x) = ones(1,length(id_x))/N;
            
            Wk_x = Wk_x / sum(Wk_x); % normalization
            ReCount = ReCount + 1;
        end
    end
    
    if i+2 > lent
        continue % end of timespan reached
    end
    
    %% Propagation
    for j = 1:N
        
        % Particle propagation
        [~,Xtemp,~,Qd_est] = integ(dynfun.est,[tpp(i+1),tpp(i+2)],Xk(:,j),options.est,dynarg.est);
        Xk(:,j) = Xtemp(:,end) + covsmpl(Qd_est(:,:,end)); %The samples are added with process noise covariance
        
    end
    
    X_est(:,i+2) = sum(Xk.*Wk_x(ones(1,nx),:),2);
    P_est_a(:,:,i+2) = Wk_x(ones(1,nx),:).*(Xk - X_est(:,(i+2)*ones(1,N)))*(Xk - X_est(:,(i+2)*ones(1,N)))';
end         % Go through the time span

P_est = scrunch(P_est_a) + scrunch(P_est_v) + scrunch(P_est_w );

err = x_tru(:,ipp) - X_est; % need better method

eflag = []; % for future implementation
Pdy = []; % for potential future implementation

%% Assign output variables
if nargout >= 3,
    varargout{1} = {tpp};
    varargout{2} = {X_est};
    varargout{3} = {P_est};
end
if nargout >=4,
    varargout{4} = {err};
end
if nargout >=5,
    varargout{5} = {dy};
end
if nargout >=6,
    varargout{6} = Pa;
end
if nargout >=7,
    varargout{7} = Pv;
end
if nargout >=8,
    varargout{8} = Pw;
end
if nargout >=9,
    varargout{9} = P_est_a;
end
if nargout >=10,
    varargout{10} = P_est_v;
end
if nargout >=11,
    varargout{11} = P_est_w;
end
if nargout >=12,
    varargout{12} = eflag;
end
if nargout >=13,
    varargout{13} = Pdy;
end

fprintf('Particles resampled %d times!\n',ReCount)

end % function
