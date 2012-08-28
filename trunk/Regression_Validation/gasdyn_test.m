function f = gasdyn_test
% GASDYN_TEST Regression test for gasdyn.
%
% F = GASDYN_TEST runs the regression test
%
% See also GASDYN
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Benjamin Asher          08/09/2012      Original
%   Ravi Mathur             08/27/2012      Rename to conform to new
%                                           regression test format


% Initial conditions

tic

% Load relavent SPICE files to obtain planetary/asteroid states
cspice_furnsh(which('wld4525.15'))
cspice_furnsh(which('de421.bsp'))
cspice_furnsh(which('naif0009.tls'))

% Define the epoch
start = '25 Feb 2024 04:00:19.000';

% Load test data
load gasdyn_data.mat

% Define scenario epoch
epoch  = datenum(start);
tspan = linspace(0,3600,60);
rho = 5.1977e11;

RA = 325.88*pi/180;
DEC = 72.71*pi/180;
w = 1136.8421/86400*pi/180;

PRA = 4.6511;

SPIN = [RA DEC PRA w]';

SPICE_ID1 = '1000109';

mass1 = 1000;
area = 20;
Cr = 1.0;

dynfun = @testdyn;

dynOpts = odtbxOptions('force');
dynOpts = setOdtbxOptions(dynOpts, 'epoch', epoch);
dynOpts = setOdtbxOptions(dynOpts, 'cR', Cr);
dynOpts = setOdtbxOptions(dynOpts, 'mass', mass1);
dynOpts = setOdtbxOptions(dynOpts, 'srpArea', area);
dynOpts.SPICE = SPICE_ID1;
dynOpts.tr = tr;
dynOpts.rho = rho*2;
dynOpts.RA = SPIN(1);
dynOpts.DEC = SPIN(2);
dynOpts.PRA = SPIN(3);
dynOpts.w = SPIN(4);
dynOpts.T_sig = 0;

x0 = [0;1.5;0;0;0;0.00011];

[~,X] = integ(dynfun,tspan,x0,[],dynOpts);

% Uncomment to generate regression test data
% X_test = X;
% save gasdyn_data.mat X_test tr 


% Compare generated and test data
tol = 1e-9;
dX = abs(X_test-X);

if dX < tol
    f = 0;
else
    f = 1;
    fprintf('OPT_SCHED regression test failed! dSens = %g\n\n',max(max(max(dX))))
end
end

function [xdot A Q] = testdyn(t,x,options)
lent = length(t);
if nargout > 1
    [xdot1(1:6,:,:) A1(1:6,1:6,:)] = polydyn(t,x,options);
else
    xdot1 = polydyn(t,x,options);
end


if nargout > 1
    [xdot2(1:6,:,:) A2(1:6,1:6,:)] = gasdyn(t,x,options);
    A = A1+A2;
else
    xdot2 = gasdyn(t,x,options);
%     xdot2(4:6,:,:) = xdot2(4:6,:,:)/1000;
end



Po = 4.51e-6/1000;
au = 1.49548e8;
cR = options.cR;
srpArea = options.srpArea;
mass = options.mass;
SPICE_ID = options.SPICE;
tstr = datestr(options.epoch + ...
    t/86400,'dd mmm yyyy HH:MM:SS.FFF');
et = cspice_str2et(tstr);
xdot(1:6,:,:) = xdot1+xdot2;
Rs = cspice_spkezr(SPICE_ID, et, 'J2000','NONE', 'SUN') + x(1:6,:);

    
for i = lent:-1:1
    [uRs rs] = unit(Rs(1:3,i));
    asc = cR*srpArea/mass*Po*au^2*Rs(1:3,i)/rs^3;
    if nargout > 1
        dvdr = cR*srpArea/mass*Po*au^2/rs^3*(eye(3,3)-3*(uRs*uRs'));
        A(:,:,i) = A(:,:,i) + [zeros(3,6); dvdr zeros(3,3)];
        Q(:,:,i) = diag([0 0 0 1e-9 1e-9 1e-9]).^2;
    end
    xdot(:,:,i) = xdot(:,:,i)+[zeros(3,1); asc];
end

end