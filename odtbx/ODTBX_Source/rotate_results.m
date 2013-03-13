function [xhat_out,e_out,P_out,sigsa_out] = ...
    rotate_results(xhat,e,P,sigsa,flag)
% 3D rotation for estimator output
%
% ROTATE_RESULTS
% Rotates Monte Carlo state estimates, covariances, and linear 
% sensitivities into a specified frame using the DCM function.
%
% [xhat_out,e_out,P_out,sigsa_out] = ROTATE_RESULTS(t,x,options) returns
% rotated state estimate, errors, covariance, and sensitivities in the
% frame specified by the 'flag' option.
%
%   INPUTS
%   VARIABLE        SIZE        DESCRIPTION
%      xhat         {N}(6xm)	Monte-Carlo state estimates [pos;vel]
%      e            {N}(6xm)    Monte-Carlo state errors
%      P            {N}(21xm)   Monte-Carlo covariance matrices
%      sigsa        (6x6xm)     Linear sensitivites
%      flag         'XXX'       DCM flag (see DCM for more information)
%
%   OUTPUTS
%      xhat_out     {N}(6xm)	State estimates in new frame [pos;vel]
%      e_out        {N}(6xm)    State errors in new frame
%      P_out        {N}(21xm)   Covariance matrices in new frame
%      sigsa        (6x6xm)     Linear sensitivites in new frame
%
% keyword: Coordinate Transformations, DCM, Utilities
% See also: DCM, JATDCM
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
%   Kenneth Getzandanner    03/17/2013      Original rotate_results.m

for i=length(xhat):-1:1
    x_tru{i} = xhat{i} - e{i};
    C_IR{i} = dcm(flag,x_tru{i}(1:3,:),x_tru{i}(4:6,:));
end

for j=length(xhat):-1:1
    for i=size(x_tru{1},2):-1:1
        
        A = blkdiag(C_IR{j}(:,:,i),C_IR{j}(:,:,i));
        
        xhat_out{j}(:,i) = A*xhat{j}(:,i);
        e_out{j}(:,i) = A*e{j}(:,i);
        
        temp = unscrunch(P{j}(:,i));
        temp = A*temp*A';
        
        P_out{j}(:,i) = scrunch(temp);
        
        if ~isempty(sigsa) && i==1
            sigsa_out(:,:,i) = A*sigsa(:,:,i)*A';
        end
    end
end