function [xdot A Q] = cameradyn(~,x,options)
% CAMERADYN Models dynamics for camera calibration parameters
%
% [XDOT A Q] = CAMERADYN(~,X,OPTIONS) returns the derivative, jacobian, and
% the process noise spectral density for camera model parameters based on
% the state X and OPTIONS.
%
% keyword: measurement
% See also OPNAVMEAS, CAMERA, ATTITUDE, BODY
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
%   Author      		    Date         	Comment
%   Kenneth Getzandanner    09/07/2011      Original cameradyn.m


% Get sigmas and time constants from ODTBX options structure
attsig = getOdtbxOptions(options,'AttitudeSigma',1e-9);
attT = getOdtbxOptions(options,'AttitudeTimeConstant',3600);
bsig = getOdtbxOptions(options,'OpticalBiasSigma',1);
bT = getOdtbxOptions(options,'OpticalBiasTimeConstant',600);

% Set process noise based on steady-state noise
attQ = 2*attsig^2/attT;
bQ = 2*bsig^2/bT;

% Check if CCD model is defined
isCCD = getOdtbxOptions(options,'isCCD',1);

% Calculate derivative, jacobian, and PSD
if ~isCCD
    xdot = [0;
            -x(2:4,1)/attT;
            -x(5:6,1)/bT];
    
    if nargout > 1
        A = zeros(6,6);
        A(2:4,2:4) = -1/attT*eye(3,3);
        A(5:6,5:6) = -1/bT*eye(2,2);
    end
    
    if nargout > 2
        Q = zeros(6,6);
        Q(2:4,2:4) = attQ*eye(3,3);
        Q(5:6,5:6) = bQ*eye(2,2);
    end
else
    xdot = [zeros(8,1);
            -x(9:11,1)/attT;
            -x(12:13,1)/bT];
    
        
    if nargout > 1
        A = zeros(13,13);
        A(9:11,9:11) = -1/attT*eye(3,3);
        A(12:13,12:13) = -1/bT*eye(2,2);
    end
    
    if nargout > 2
        Q = zeros(13,13);
        Q(9:11,9:11) = attQ*eye(3,3);
        Q(12:13,12:13) = bQ*eye(2,2);
    end
end

end