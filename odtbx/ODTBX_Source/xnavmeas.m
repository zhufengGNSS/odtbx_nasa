function [Y,H,R] = xnavmeas(t,X,options)
% XNAVMEAS Calculates spacecraft to distant RSO range difference & range difference rate measurements. 
%
% Calculate range difference and range difference rate between a spacecraft
% and a distant inertially fixed Resident Space Object (RSO).
%
%   [y,H,R] = xnavmeas(t,X,options)
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each type are grouped together in rows.
% Each range_i, rangeRate_i pair corresponds to measurement i = 1:num_meas
% Thus, using the default options settings, the output y will look like:
%   y = [  range_1(t1)      range_1(t2)     ... ;
%          rangeRate_1(t1)  rangeRate_1(t2) ... ;
%          range_2(t1)      range_2(t2)     ... ;
%          rangeRate_2(t1)  rangeRate_2(t2) ... ;
%              :                  :         ... ;
%          range_I(t1)      range_I(t2)     ... ;
%          rangeRate_I(t1)  rangeRate_I(t2) ... ]
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t        (1xN)	Times from epoch (secs)
%      X        (6xN)   States (km & km/s) of user satellite at times t
%      options  (1x1)   data structure
%
% OUTPUTS
%      y        (MIxN)   measurements
%      H        (MIx6xN) measurement partials matrix
%      R        (MIxMxN) measurement covariance
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. Default values are 
% enclosed in parentheses. The options parameters that are valid for this 
% function are:
%
%   PARAMETER             DESCRIPTION                   VALID VALUES             
%  xnav.useRangeRate     flag to output range diff     {true(default), false}
%                        rate                
%  xnav.obs_time	     observation time in secs      t_obs > 0 (1e5  secs)        
%  xnav.num_meas         Number of measurements to     1 <= num_meas (1)
%                        RSO
%                        
%  xnav.sigma_r          1xnum_meas array of            sigma_r >0 (0.0136640)
%                        range measurement noise for   
%                        each k=1:num_meas RSO in Km	                         
%  xnav.sigma_rr         1xnum_meas array of measurement standard 
%                        deviations corresponding to range     sigma_rr>=0 (1.9324e-007)
%                        rate measurement noise for each 
%                        k=1..num_meas RSO in	Km/s	
%  xnav.rso_eci          3xnum_meas array of unit vectors
%                        to the RSO's.  
%                        DEFAULT values are given below:
%
%    rso_eci =  [ 0.104817533369789   0.094966457808533   0.389831073451460   0.607045484518442
%                 0.921144947109562  -0.902283088049447  -0.844226330971047   0.417623670844230
%                 0.374840327489981  -0.420555110432931   0.367850018712583   0.676081540404973];
%
%  xnav.C_1              1xnum_meas array of RSO time of 
%                        arrival (TOA) error attribute
%                        DEFAULT values are given below:
%    C_1 = [0.000000000207737   0.000050430244604   0.000002312225139   0.004047715093601];
%
%   Schedule          [(ncx3) matrix]        XNAV Tracking Schedule. 
%               Schedule restricts the xnavmeas measurement model to only
%               provide measurements for specific pulsars during specific 
%               time intervals. The format for each row is:
%               [ pulsar_index start_time stop_time ]
%               The start_time and stop_time must be in seconds from epoch.
%               The matrix can be any length (nc = number of contacts)
%       Example: 
%       Sched = {1 '28-Jan-2010 06:56:35' '28-Jan-2010 12:56:35'
%                2 '28-Jan-2010 12:56:35' '28-Jan-2010 19:56:35'
%                1 '28-Jan-2010 19:56:35' '29-Jan-2010 02:56:35'};
%       Sched = cell2mat(TrackSched(:,1)); % pulsar numbers
%       Sched(:,2) = (datenum(TrackSched(:,2))-epoch)*86400; %start(epsecs)
%       Sched(:,3) = (datenum(TrackSched(:,3))-epoch)*86400; %end(epsecs)
% 
% VALIDATION/REGRESSION TEST
%  The validation/regression testing capability has been extracted to
%  xnavmeas_test.m in the regression testing framework.
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, RRDOTLT, GSMEAS, 
%  ODTBXOPTIONS
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Ksenia Kolcio       03/25/2010      Original xnavmeas.m
%   Kevin Berry         04/06/2010      Added xnav options structure
%   Sun Hur-Diaz        07/08/2010      Added Sched option
%   Sun Hur-Diaz        09/22/2010      Added fields to xnav structure
%   Ravi Mathur         08/27/2012      Extracted regression test

% specify defaults here 
t_obs_def = 1e5; % sec

% measurement characteristic parameter for noise calc.
C_1 = [0.000000000207737   0.000050430244604   0.000002312225139];%   0.004047715093601];

% ECI LOS unit vectors to RSO.
l_rso_eci =  [ 0.104817533369789   0.094966457808533   0.389831073451460%   0.607045484518442
               0.921144947109562  -0.902283088049447  -0.844226330971047%   0.417623670844230
               0.374840327489981  -0.420555110432931   0.367850018712583];%   0.676081540404973];

% Set the default options           
num_meas = 3; 
t_obs =  t_obs_def ; % sec
useRangeRate = true ;
sigma_r_user = []; % let this be empty if none is specified, they will be calcuated
sigma_rr_user =[]; 

% Get the options data or use defaults if none are specified,
% and constants
xnav = getOdtbxOptions(options, 'xnav',[]); %ddor structure      
if isfield(xnav,'useRangeRate'),    useRangeRate  = xnav.useRangeRate;  end
if isfield(xnav,'obs_time'),        t_obs         = xnav.obs_time;      end
if isfield(xnav,'num_meas'),        num_meas      = xnav.num_meas;      end
if isfield(xnav,'sigma_r'),         sigma_r_user  = xnav.sigma_r;       end
if isfield(xnav,'sigma_rr'),        sigma_rr_user = xnav.sigma_rr;      end
if isfield(xnav,'rso_eci'),         l_rso_eci     = xnav.rso_eci;       end
if isfield(xnav,'C_1'),             C_1           = xnav.C_1;           end
if num_meas ~= size(l_rso_eci,2) || num_meas ~= length(C_1)
    error('Inconsistent number of RSOs specified.')
end
Sched = getOdtbxOptions(options,'Schedule',[]); % Contact schedule

% JAT constants
c_km = JATConstant('c')*1e-3; %speed of light in kilometers

% vector/matrix sizes
N = length(t); 
M = 1+ useRangeRate; % number of measurement types; range diff is always an output
P = size(X,1); 

% Error checking user input

% If length of user specified range sigma vector does not match num_meas
% use default values
if ((length(sigma_r_user) ~= num_meas) && (~isempty(sigma_r_user)))
    warning('XNAVMEAS:Inconsistency',...
    ' Number of user-specified range measurement noise sigmas does not match number of user-specified measurements so calculated values will be used')
    sigma_r_user = [];
end

% Size output matrices and create pulsar LOS vectors
Y = NaN(M*num_meas,N);
H = zeros(M*num_meas,P,N);
R = zeros(M*num_meas,M*num_meas,N);

sigma_r = zeros(1,num_meas);
sigma_rr = zeros(1,num_meas);

% Loop on the pulsars 1:K
i=1;
for k=1:num_meas
    
    % Check schedule for times when tracking is done
    if ~isempty(Sched)
        gSch = Sched(k==Sched(:,1),2:3);
        tind = [];
        for m=1:size(gSch,1)
            tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
        end
    else
        tind = 1:length(t);
    end
    t1  = t(tind);
    x1  = X(:,tind);
    
    len=length(t1);
    
    if ~isempty(t1)
        
        % Range
        rdiff_k = l_rso_eci(:,k)'*x1(1:3,:);
        
        % Range Measurement Jacobian
        Hr_k = [l_rso_eci(:,k)' zeros(1,3)];
        
        % Measurement Covariance
        % Rr_k
        if isempty(sigma_r_user)
            % noise variance based on observation time and measurement
            % characteristics
            var_toa_k = C_1(k)/t_obs;
            R_k = c_km^2*var_toa_k;
            sigma_r(k) = sqrt(R_k);
            
        else % user-specified sigmas
            R_k = sigma_r_user(k)^2;
            sigma_r(k) = sigma_r_user(k);
        end
        if useRangeRate
            rdiffr_k = l_rso_eci(:,k)'*x1(4:6,:);
            if isempty(sigma_rr_user)
                sigma_rr(k) = sqrt(2)*sigma_r(k)/t_obs;
                R_k = [sigma_r(k)^2 sigma_rr(k)^2];
            else % user specified sigmas
                sigma_rr(k) = sigma_rr_user(k);
                R_k = [sigma_r(k)^2 sigma_rr(k)^2];
            end
            % Range Diff rate Measurement Jacobian
            Hv_k = [zeros(1,3) l_rso_eci(:,k)' ];
        else
            rdiffr_k = [];
            Hv_k = [];
        end
        
        % stack the measurements for pulsar k into output matrix Y
        Y(i:M*k,tind) = [rdiff_k;rdiffr_k];
        H(i:M*k,1:P,tind) = reshape(repmat([Hr_k;Hv_k],1,length(tind)),[M,P,len]);
        R(i:M*k,i:M*k,tind) = reshape(repmat(diag(R_k),1,length(tind)),[M,M,len]);
        
    end
    
    i=M*k+1; % advance row index to next RSO by number of measurements

end % loop over num_meas

end