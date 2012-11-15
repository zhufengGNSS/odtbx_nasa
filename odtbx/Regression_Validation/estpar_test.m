function failed = estpar_test

% Particle Filter regression test
% keyword: estpar regression
% See also estpar
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

% REVISION HISTORY
%  Author:                  Date:          Comments:
%  John Gaebler         June 2012
%  Ravi Mathur             08/27/2012      Rename to conform to new
%                                          regression test format

failed = 0;
tol = 1e-9;

disp('Running ESTPAR regression test case 1...')

% Now check the nonlinear example of the planar 2-body problem

load estpar_reg_data.mat

% It uses the restricted 2-body problem 'r2bp' for the dynamics model, and 
% the range and range rate from two ground stations located in 
% Antarctica, Syoma and Sweden, Kiruna using 'gsmeas' for the 
% measurement model.
%
% Specify the dynamics models for the truth and the estimator.  Here, we 
% are going to use the same models for both.
dynfun.tru = @r2bp; 
dynfun.est = @r2bp;
%
% Specify the inputs to the dynamics models
mu = 3.986004415e5; % km^3'sec^2
dynarg.tru = mu;
dynarg.est = mu;

% Define the initial reference state and a priori covariance matrix:
%xo = [-29793.471825;21600.776094;-39473.116501;0.497358;-0.676146;-1.684046]; %km  & km/s Polar elliptical orbit      
xo = [-3.061505424000000e+005;2.651129163110612e-011;2.651129163110612e-011;-5.869440191664465e-016;-1.920887429475406e-001;-1.920887429475406e-001]; %km  & km/s Circular orbit
%
Po = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6].^2); % km^2 & km^2/s^2

% Set up the solve-for and consider mapping
S = eye(6);                       % Solve-for map - solve for all 6 states
C=[];                         % Consider map - no consider states
 
% The next four lines say to use the same a priori for truth and estimator:
Xnot.Xo = xo;                 % True initial state
Xnot.Xbaro = S*xo;            % Estimator initial state
Pnot.Po = Po;                 % True initial covariance
if isempty(Po)
    Pnot.Pbaro = [];
else
    Pnot.Pbaro = S*Po*S';         % Estimator initial covariance
end
%
% Specify the measurement models for the truth and the estimator.  Here, we 
% are going to use the same models for both.
datfun.tru = @gsmeas;%
datfun.est = @gsmeas;%

% Specify the inputs to the measurement models, which in this case is the 
% epoch date and time, ground station IDs and the measurement standard
% deviations

Sched = [...
    1   60      300
    2   360    660
    3   720    1020
    ];

EpochString = 'July 9, 2010 16:00:00';
epoch = datenum(EpochString);
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
datarg = setOdtbxOptions(datarg,'gsID',{'SYOQ','KI2S','MCMS'}); %,'MC1S','BLKQ'
datarg = setOdtbxOptions(datarg,'rSigma',[1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]); % 1e-2 1e-5 1e-2 1e-5
datarg = setOdtbxOptions(datarg,'gsList',createGroundStationList());
datarg = setOdtbxOptions(datarg,'Schedule',Sched);

tspan = 0:60:1140;

% Specify the number of particles/samples for the filter
N =20;

eopts = odtbxOptions('estimator');
eopts = setOdtbxOptions(eopts,'MonteCarloCases',1);
eopts = setOdtbxOptions(eopts,'UpdateIterations',1);
eopts = setOdtbxOptions(eopts,'MonteCarloSeed',1);
eopts = setOdtbxOptions(eopts,'EditRatio',ones(6,1)*9);
eopts = setOdtbxOptions(eopts,'EditFlag',ones(6,1));
eopts = setOdtbxOptions(eopts,'Particles',N);

% Now run (and time) the particle filter:

[~,X,P,E,DY,PA,PV,PW,PHATA,PHATV,PHATW]=estpar...
    (dynfun,datfun,tspan,Xnot,Pnot,eopts,dynarg,datarg,S,C);

% Uncomment to generate test data per the request of the ODTBX dictator
% Xs = X; Ps = P; Es = E; DYs = DY; PAs = PA; PVs = PV; PWs = PW; PHATAs = PHATA; PHATVs = PHATV; PHATWs = PHATW;
% save('estpar_reg_data.mat','Xs','Ps','Es','DYs','PAs','PVs','PWs','PHATAs','PHATVs','PHATWs')

if any(any(any(abs(cell2mat(P) - cell2mat(Ps)) > tol)))
    warning('estpar_regression: Phat mismatch.');
    failed = 1;
end
if any(any(any(abs(PA - PAs) > tol)))
    warning('estpar_regression: PA mismatch.');
    failed = 1;
end
if  any(any(any(abs(PV - PVs) > tol)))
    warning('estpar_regression: PV mismatch.');
    failed = 1;
end
if any(any(any(abs(PW - PWs) > tol)))
    warning('estpar_regression: PW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATW - PHATWs) > tol)))
    warning('estpar_regression: PHATW mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATV - PHATVs) > tol)))
    warning('estpar_regression: PHATV mismatch.');
    failed = 1;
end
if any(any(any(abs(PHATA - PHATAs) > tol)))
    warning('estpar_regression: PHATA mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(X) - cell2mat(Xs)) > tol)))
    warning('estpar_regression: X mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(E) - cell2mat(Es)) > tol)))
    warning('estpar_regression: E mismatch.');
    failed = 1;
end
if any(any(any(abs(cell2mat(DY) - cell2mat(DYs)) > tol)))
    warning('estpar_regression: DY mismatch.');
    failed = 1;
end

end
% End of function.
