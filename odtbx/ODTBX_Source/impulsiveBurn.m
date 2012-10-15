function restartRecord = impulsiveBurn(varargin)
% IMPULSIVEBURN  Impulsive Maneuver (for RESTARTRECORD).
%
%   [RESTARTRECORD] = IMPULSIVEBURN(RESTARTRECORD,DV,REL,ABS) applies the
%   impulsive maneuver DV to the RESTARTRECORD with the relative execution 
%   error REL (magnitude only) and the absolute execution error ABS 
%   (magnitude and direction). The size of DV is [3xN+1] where N is the 
%   number of Monte Carlo cases in RESTARTRECORD.  DV[:,1] is associated 
%   with the linear covariance reference trajectory.  DV[:,2:N+1] is 
%   associated with Monte Carlo case trajectories [1:N].
%
%   keyword: Impulsive, Maneuver, Restart
%
%   See also
%      Estimation:      ESTSEQ
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

% Kenneth M. Getzandanner
% NASA Goddard Space Flight Center
%
% Modification History
% ---------------------
% Kenneth Getzandanner    10/01/2010      Original impulsiveBurn.m
%
% Kenneth Getzandanner    04/05/2011      Added relative and absolute
%                                         error as inputs
%
% Kenneth Getzandanner    05/12/2011      Added optional input flag to
%                                         reset a priori covariance 

% Read function input
restartRecord = varargin{1};
dv = varargin{2};
rel = varargin{3};
abs = varargin{4}; 

if length(varargin)>4
    flag = varargin{5};
else
    flag = [];
end

% Loop through the linear covariance (i=1) and Monte Carlo (i=2:n+1)
% reference trajectories
for i=1:size(dv,2)
    
    % Define the covariance associated with the maneuver execution error
    Pdv = zeros(6,6);
    if norm(dv(:,i))>0
        
        u = unit(dv(:,i));
        y = norm(dv(:,i));
        
        % Define the maneuver execution error relative to the DV direction
        % (Gibbs method)
        Pg = diag([(abs + rel*y) abs abs].^2);
        
        % Define the DCM to rotate from the 'DV' frame to inertial
        ref = [0;0;1];
        
        m = u;
        n = unit(cross(ref,u));
        l = unit(cross(m,n));
        
        DCM(1,:) = m;
        DCM(2,:) = n;
        DCM(3,:) = l;
        
        % Rotate the maneuver execution covariance into the inertial frame
        Pdv(4:6,4:6) = DCM'*Pg*DCM;
    end
    
    % Sample from the maneuver execution covariance to calculate the
    % maneuver execution error
    xerr = covsmpl(Pdv(4:6,4:6));
    
    % Linear Covariance portion of the restart record
    if i == 1
        % Apply DV to the estimate/true reference trajectory
        restartRecord.Xsrefo(4:6,1) = restartRecord.Xsrefo(4:6,1) + dv(:,i);
        restartRecord.Xrefo(4:6,1) = restartRecord.Xrefo(4:6,1) + dv(:,i);
        
        % Add maneuver execution contributions to the estimate/true 
        % covariance
        restartRecord.Pmo(1:6,1:6,1) = restartRecord.Pmo(1:6,1:6,1) + Pdv;
        restartRecord.Phatmo(1:6,1:6,1) = restartRecord.Phatmo(1:6,1:6,1)...
            + Pdv;
        
        % If the 'clear' flag is set, assign all covariance as a priori and
        % set process and measurement covariance partitions to zero
        if strcmp(flag,'clear')
            restartRecord.Pao = restartRecord.Pao + restartRecord.Pvo +...
                restartRecord.Pwo;
            
            restartRecord.Pvo = zeros(size(restartRecord.Pvo));
            restartRecord.Pwo = zeros(size(restartRecord.Pwo));
        end
    
    % Monte Carlo portion of the restart record
    else
        % Apply DV to the estimate/true reference trajectory
        restartRecord.Xhato(4:6,i-1) = restartRecord.Xhato(4:6,i-1)+dv(:,i);
        restartRecord.Xo(4:6,i-1) = restartRecord.Xo(4:6,i-1)...
            + dv(:,i) + xerr;
        
        % Add maneuver execution contributions to the estimate covariance
        restartRecord.Phato(1:6,1:6,i-1) = restartRecord.Phato(1:6,1:6,i-1)...
            + Pdv;
    end
end

end