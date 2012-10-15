function [xDot,A,Q] = nbodypm(t,x,options)
% NBODYPM N-Body Point Mass dynamics equations.
% 
%   xDot = r2bp(t,x,options) returns the derivatives of the nbody point
%   masses force models specified in the given options structure. 
% 
%   [xDot,A] = r2bp(t,x,options) also returns the Jacobian (the state
%   derivatives of the equations of motion) if possible; otherwise an []
%   empty matrix is returned.
% 
%   [xDot,A,Q] = r2bp(t,x,options) also returns the process noise matrix 
%   that is hardcoded within this file.
% 
%   
% 
%   INPUTS
%   VARIABLE       SIZE    		DESCRIPTION (Optional/Default)
%      t           (1xn)       Time since start of simulation (secs)
%      x           (6xn)       Input state (x,y,z position in kilometers 
%                              followed by x,y,z velocity in km/s)
%      options     structure   Options structure created by
%                              initialize_nbody
% 
% 
%   OUTPUTS 
%      xDot        (6xn)       Derivatives of state (kilometers, seconds)
%      A           (6x6xn)     State transition matrix (kilometers, seconds)
%      Q           (6x6xn)     Process noise (kilometers, seconds)
% 
% 
%    keyword: nbody, point mass, force model 
%    See also ODTBXOPTIONS, INITIALIZE_NBODY
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

% Kevin Berry
% NASA Goddard Space Flight Center
%
%  REVISION HISTORY
%   Author     			Date         	Comment
%               		(MM/DD/YYYY)
%   Kevin Berry          05/10/2011     Original
%   Benjamin Asher       05/10/2011     Added initialize_nbody and changed
%                                       the inputs accordingly

epoch = getOdtbxOptions(options,'epoch',[]);
PointMasses = getOdtbxOptions(options,'PointMasses',[]);
CentralBody = getOdtbxOptions(options,'CentralBody',[]);
GM_PM = getOdtbxOptions(options,'GM_PM',[]);
GM_CB = getOdtbxOptions(options,'GM_CB',[]);

%% Get the point mass locations
et = cspice_str2et(datestr(epoch, ...
    'dd mmm yyyy HH:MM:SS.FFF')); %Converts the epoch to the number of TDB seconds past the J2000 epoch
et = et+t(:)'; %Adds the simulation times to the epoch
for n=length(PointMasses):-1:1
    X_pm{n} = cspice_spkezr(PointMasses{n}, et, 'J2000', 'NONE', CentralBody);
end

%% Combine the gravity contributions from each body
I = eye(3);
R = x(1:3,:); %User position vectors
V = x(4:6,:); %User velocity vectors
[Ur,r] = unit(R); %Unit vectors and magnitudes of user positions
G = ([1;1;1]*(-GM_CB./r.^2)).*Ur; %Gravity contribution from central body
A = zeros(6,6,length(t));
if nargout > 1
    for k = length(t):-1:1
        A(4:6,1:3,k) = -GM_CB./r(k).^3*(I - 3*Ur(:,k)*Ur(:,k)');
        A(1:3,4:6,k) = I;
    end
end
for n=1:length(PointMasses)
    %Gravity contribution from each point mass
    R_bn = X_pm{n}(1:3,:); %Body n position vectors wrt cb
    [Ur_bn,r_bn] = unit(R_bn); %Unit vectors and magnitudes of body n wrt cb
    [Ur_ubn,r_ubn] = unit(R_bn-R); %Unit vectors and magnitudes of user wrt body n
    G = G + GM_PM{n} * ( ([1;1;1]*r_ubn.^-2).*Ur_ubn - ([1;1;1]*r_bn.^-2).*Ur_bn );
    if nargout > 1
        for k = length(t):-1:1
            A(4:6,1:3,k) = A(4:6,1:3,k) - GM_PM{n}./r_ubn(k).^3*(I - 3*Ur_ubn(:,k)*Ur_ubn(:,k)');
        end
    end
end
xDot = [V; G];

if nargout == 3,
    Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2),[1 1 length(t)]);
end

