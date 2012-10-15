function failed = jatWorldPropagatorRKF78_case1_test()
% Regression Test Case
% Function(s) jatWorldPropagatorRKF78, jatForces, createJatWorld
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
%
%   Ravi Mathur         08/29/2012      This test fails, and needs to be
%                                       fixed or removed pending CCB.

failed      = 0;
tol         = 1e-6;
integrator  = @jatWorldPropagatorRKF78;

% Set time period
tspan   = 0:300:86400;

% Set initial state vector [x y z xdot ydot zdot] in km & km/s
x0 = [-4453783.586; -5038203.756; -426384.456; 3831.888; -2887.221; -6018.232]/1000;

% Select Integrator
odeopts     = odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10);
odeOptions  = setOdtbxOptions('OdeSolver',integrator,'OdeSolvOpts',odeopts);

% Initialize options structure
jOptions = jatForcesRegressionOptions();

% jatWorld is an ODToolboxJATModel object that stores relevant force and
% spacecraft information
jatWorld = createJATWorld(jOptions);

% Run simulation
x0 = x0*1000;  % using JAT integrators, need to convert input state to m

[tPlot,xPlot] = integ(@jatForces,tspan,x0,odeOptions,jatWorld);

xPlotTest = xPlot/1000;  % convert back to km for comparisons

load jatWorldPropagatorRKF78RegressionBaseline

xDiff = xPlotTest - xPlot;

if( any( any( abs(xDiff)>tol ) ) )
    failed = 1;
end