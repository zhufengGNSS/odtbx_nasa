function f = opt_sched_test
% OPT_SCHED_TEST Regression test for opt_sched.
%
% F = OPT_SCHED_TEST runs the regression test
%
% See also OPT_SCHED
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
tspan = 0:60:3600;
Xo = [3677.25392172411;3128.98832260195;5068.2576564521;-6.42096333845847;2.08243295166916;3.37307323108473];
mu = 3.986004415e5;
sig_rad = 1e-3;
sig_spd = 1e-6;
Po = initialize_cov(Xo(1:3),Xo(4:6),sig_rad,sig_spd,mu);
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'MonteCarloSeed',1);
S = eye(6);
C = [];
Xnot.Xo = Xo;
Xnot.Xbaro = S*Xo;
Pnot.Po = Po;
Pnot.Pbaro = S*Po*S';
dynfun.tru = @r2bp;
dynfun.est = dynfun.tru;
dynarg.tru = mu;
dynarg.est = dynarg.tru;
datfun.tru = @simplemeas;
datfun.est = datfun.tru;
datarg.tru = [];
datarg.est = datarg.tru;

% Generate test data
[~,~,P1] =estbat(dynfun,datfun,tspan,Xnot,Pnot,eOpts,dynarg,datarg,S,C);
Pa = unscrunch(P1{1}(:,:));
Sens = opt_sched(dynfun.tru,datfun.tru,tspan,Xo,eOpts,dynarg.tru,datarg.tru,Pa(:,:,end));

% Uncomment to generate regression test data
% Sens_test = Sens;
% save opt_sched_data.mat Sens_test

% Load test data
load opt_sched_data.mat

% Compare generated and test data
tol = 1e-9;
dSens = abs(Sens_test-Sens);

if dSens < tol
    f = 0;
else
    f = 1;
    fprintf('OPT_SCHED regression test failed! dSens = %g\n\n',max(max(max(dSens))))
end
end

function [y H R] = simplemeas(~,x,~)
y = x;
lentx = size(x,2);
for ii = lentx:-1:1   
    H(:,:,ii) = eye(6);
    R(:,:,ii) = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2);
end
end