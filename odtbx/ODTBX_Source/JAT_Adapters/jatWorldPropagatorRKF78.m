function [tOut,yOut] = jatWorldPropagatorRKF78(ode,tspan,y0,odeSolverOptions,varargin)

% JATWORLDPROPAGATORRKF78 JAT Runge-Kutta-Fehlberg variable step integrator (Earth orbits)
%
% [tOut,yOut] = jatWorldPropagatorRKF78(ode,tspan,y0,odeSolverOptions,varargin)
%   integrates the ode over the tspan.  The variable step integrator uses the
%   'InitialStep' parameter within odeSolverOptions as the initial step size.
%   The 'RelTol' parameter is used as the accuracy tolerance. The minimum
%   step size is hardcoded within this file to be 1e-9.
%
%   Output is provided at the points input in tspan. If the output
%   timepoints don't match the stepsize, Lagrangian interpolation is used 
%   to compute the values at the output times.
%
% jatWorldPropagatorRKF78 only works with the JAT force models and is 
% intended for Earth-centered orbits with a single spacecraft. As such, the
% inputs ode and dynfun are ignored.  The forces specified in jatWorld are
% used. 
%
%   INPUTS 
%   VARIABLE            SIZE    DESCRIPTION (Optional/Default)
%     ode               (1x1)   IGNORED
%     tspan             (1xN)   Vector of input times (seconds)
%     y0                (6x1)   Initial state vector [x y z xdot ydot zdot]
%                               in meters and m/s
%     odeSolverOptions  (1x1)   ode options structure (see odeset)
%     varargin                  cell array containing following:
%       options         (1x1)   estimator options structure (see ODTBXOPTIONS)
%       dynfun          (1x1)   IGNORED - handle to dynamics function
%       neqns           (1x1)   number of states/equations
%       jatWorld        (1X1)   java object created by createJATWorld
%
%   OUTPUTS 
%     tOut              (Nx1)   Vector of output times (seconds)
%     yOut              (Nx6)   Output states (meters and m/s)
%
%   keyword: integrator
%   See also jatWorldPropagatorRK4, jatWorldPropagatorRK8, JATFORCES, 
%   CREATEJATWORLD
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Derek Surka         08/13/2007   	Original
%   Allen Brown         02/10/2009      Renamed to jatWorldPropagatorRKF78

if nargin < 8
    error('MATLAB:jatWorldPropagatorRKF78:NotEnoughInputs',...
        'An ODToolboxJATModel object from createJATWorld must be input.');
end

jatWorld = varargin{4};

% need to verify that this is a JAT integrator with a JAT force model
if( ~isjava(jatWorld) || ...
        ~strcmp(class(jatWorld),'jat.matlabInterface.ODToolboxJATModel') )
    error('MATLAB:jatWorldPropagatorRKF78:NoJATObject An ODToolboxJATModel object must be the eighth input.');
end

stepSize    = odeget(odeSolverOptions,'InitialStep',10,'fast');
accuracy    = odeget(odeSolverOptions,'RelTol',1e-9,'fast');
minStepSize = 1e-9;

% Now do integrator stuff...
jatStates = jat.matlabInterface.ODToolboxIntegrators.RKF78(tspan,y0,stepSize,minStepSize,accuracy,jatWorld);

tOut = jatStates(:,1);
yOut = jatStates(:,2:end);

return

