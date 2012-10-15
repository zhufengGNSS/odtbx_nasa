function [y H R] = opnavmeas(t,x,options)
% OPNAVMEAS  Simulates measurements from optical landmark tracking.
%
% [y,H,R] = OPNAVMEAS(tspan,x,options) generates measurements of landmark
% positions in the camera focal plane for a given celestial body, camera,
% and attitude profile defined in 'options'.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch)
%      x            (6xN)   Body-centered spacecraft state [pos;vel] (km)
%      options      (1x1)   data structure, see below
%
%   OUTPUTS
%      y            (MxN)   measurements
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%
% The measurements are output in y. Each column corresponds to a different
% time. Each landmark produces two measurements corresponding to the 
% horizontal and vertical distance from the boresite in the camera focal
% plane.  Therefore, measurements from 3 landmarks will produce:
%
%   y = [   lmk1_x(t1)   lmk1_x(t2)...;
%           lmk1_y(t1)   lmk1_y(t2)...;
%           lmk2_x(t1)   lmk2_x(t2)...;
%           lmk2_y(t1)   lmk2_y(t2)...;
%           lmk3_x(t1)   lmk3_x(t2)...;
%           lmk3_y(t1)   lmk3_y(t2)... ]
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES         NOTES
%   Camera                Camera object      See CAMERA
%   Attitude              Attitude object    See Attitude
%   OpticalSigma          Scalar>=0          {6e-10}
%
% keyword: measurement
% See also CAMERA, ATTITUDE, BODY
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
%   Kenneth Getzandanner    03/29/2011      Original opnavmeas.m
%

% Get camera and attitude objects from options structure
camera = getOdtbxOptions(options,'Camera',[]);
attitude = getOdtbxOptions(options,'Attitude',[]);

% Get measurement error from options structure
if isa(camera,'CCDCamera')
    operr = getOdtbxOptions(options,'OpticalSigma',1);
else
    operr = getOdtbxOptions(options,'OpticalSigma',6e-10);
end

% Check camera calibration flag
cflag = getOdtbxOptions(options,'CalibrationFlag',0);

% Check if camera and attitude are valid objects
if ~isa(camera,'Camera') || ~isa(attitude,'Attitude')
    error('Error: user must specify valid camera and attitude object')
end

% Number of landmarks
n = size(camera.body.lmk,1);

% Initialize measurement vector and partials
y = NaN(n*2,length(t));
H = NaN(n*2,size(x,1),length(t));

% Loop through each timestep
for i=1:length(t)
    
    % Set camera attitude using attitude object
    camera.C_CI = attitude.C_CI(:,:,t(i) == attitude.tspan);
    
    % Set camera position using state vector
    camera.R = x(1:3,i);
    
    if cflag
        camera = setStates(x(:,i),camera);
    end
    
    % Get camera frame to generate measurements
    [pixel line] = getCameraFrame(camera,t(i));
    y(1:2:(end-1),i) = pixel;
    y(2:2:end,i) = line;
    
    % Get measurement partials
    % Check calibration flag
    if cflag
        % If CCDCamera, return partials for focal length, pixel mapping,
        % zero pixel offset, and error coefficients
        if isa(camera,'CCDCamera')
            [Hc Hf Hk Ho He Ha Hb]  = getCameraPartials(camera,t(i));
            H(:,:,i) = [Hc Hf Hk Ho He Ha Hb];
        % If Camera, return partials for focal length and attitude
        else
            [Hc Hf Ha Hb] = getCameraPartials(camera,t(i));
            H(:,:,i) = [Hc Hf Ha Hb];
        end
    else
        H(:,:,i) = getCameraPartials(camera,t(i));
    end
end

% Set measurement covariance using operr
R = repmat(eye(2*n,2*n)*operr^2,[1 1 length(t)]);

end

function camera = setStates(x,camera)

if isa(camera,'CCDCamera')
    camera.f = x(7,1);
    camera.Kx = x(8,1);
    camera.Ky = x(9,1);
    camera.s0 = x(10,1);
    camera.l0 = x(11,1);
    camera.Er = x(12,1);
    camera.Em = x(13,1);
    camera.En = x(14,1);
    camera.theta = x(15:17,1);
    camera.bias = x(16:18,1);
else
    camera.f = x(7,1);
    camera.theta = x(8:10,1);
    camera.bias = x(11:12,1);
end

end
