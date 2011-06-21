function elapsedTime = jatForcesRegression(useJAT, integrator, filename)

% This is a demonstration function of how to access the JAT force models
% and integrators for Earth-centric orbits
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

if( nargin < 1 )
    useJAT   = false;
end;

if( (nargin < 2) || ~isa(integrator,'function_handle') )
    integrator = @ode45;
end

if( nargin < 3 )
    filename = 'forcesDemo';
end

jatIntegrators = {
    'jatWorldPropagatorRK8'
    'jatWorldPropagatorRK4'
    'jatWorldPropagatorRKF78' };

if( isempty( getIndex(func2str(integrator),jatIntegrators,false) ) )
    useJatIntegrator = false;
else
    useJatIntegrator = true;
end

tic

% Set time period
tspan   = 0:300:86400;
% tspan = [0 86400];

% Set initial state vector [x y z xdot ydot zdot] in km & km/s
x0 = [-4453783.586; -5038203.756; -426384.456; 3831.888; -2887.221; -6018.232]/1000;
% x0 = [6878;0.00;0.00;0.00;0.00;8.339];

% Select Integrator
odeopts = odeset('reltol',1e-9,'abstol',1e-9,'initialstep',10); 
odeOptions = setOdtbxOptions('OdeSolver',integrator,'OdeSolvOpts',odeopts); % if a JAT integrator is selected, the matlab force model is ignored

% Simulation Demo
if( useJAT )
    % Initialize options structure
    jOptions = jatForcesRegressionOptions();
    
    % jatWorld is an ODToolboxJATModel object that stores relevant force and
    % spacecraft information
    jatWorld = createJATWorld(jOptions);

    % Run simulation
    if( useJatIntegrator ) % using JAT integrators, need to convert input state to m
        x0 = x0*1000;
    end

    [tPlot,xPlot] = integ(@jatForces,tspan,x0,odeOptions,jatWorld); 

    if( useJatIntegrator ) % convert back to km for comparisons
        xPlot = xPlot/1000;
    end   
         
    save(filename, 'tPlot', 'xPlot', 'useJAT', 'jOptions')
else
    mu            = JATConstant('muEarth','JGM')/1000^3;  % km^3'sec^2
    [tPlot,xPlot] = integ(@r2bp,tspan,x0,odeOptions,mu);
    save(filename, 'tPlot', 'xPlot', 'useJAT')
end;

elapsedTime = toc;

end
