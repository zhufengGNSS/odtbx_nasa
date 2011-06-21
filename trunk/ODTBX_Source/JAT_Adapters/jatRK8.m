function [t,y] = jatRK8(eoms, time, x0, varargin)

% JATRK8 The wrapper for the JAT Runge-Kutta 8th order integrator.
%
%    [t,y] = jatRK8(eoms, time, x0, varargin);
%
% This is a wrapper for the Runge-Kutta 8th order fixed-step integrator 
% implemented in JAT. It interprets the Java return and returns it the way 
% Matlab would return its output from the ODE integrators. It also calls 
% getJatRK8Options to check the optional input structure and make sure it 
% has all the necessary values.  
%
% This wrapper can be used as a 'generic' RK8 integrator with the client 
% specifying the equations of motion and initial state OR it can be used 
% as an Earth-centric single spacecraft propagator via jatWorld force 
% models.
%
%   INPUTS
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%   eoms       string       Function name of the EOM function to use, or 
%                           special JatUniverse names, see below
%   time       (1xM)        Start/finish times (M=2) or the specified output 
%                           times (M>2) 
%   x0         (1xN)        Initial state defined by the user, see below
%                           (Units aligned with the initial state and eoms). 
%   varagin    structure    Options structure that is setup using
%                           setJarRK8Options
%
%   OUTPUTS 
%   t          (1XM)        Output from function, output time, stepsize 
%                           dictated by user or default set within to
%                           correspond to 0.5 deg true anomaly
%   y          (NxM)        States corresponding to the output time
%
%   INPUT SELECTION AND OPERATION:
%   First you must choose how you are going to use the JAT 
%   integrator and if you will use any of the JAT force models.
%   
%   METHOD 1: Custom models and state:
%   If you are going to use your own force models then you will 
%   enter in whatever state is needed for your force model, and 
%   the first argument to the integrator will be the name of your
%   Matlab EOM file without the .m extension. 
%
%      jatRK8('EOMName', ...
%
%   Make sure the state size, state units, and time units correspond with 
%   the equations in the 'EOMName' function.  The format of the function 
%   has to be of the form: 
%
%      f=EOMName(t,x).  
%
%   **NOTE: Inside this EOM function, you must add the following statement 
%   before any calculation is done with the states:
%
%      if iscell(x)
%         x=cell2mat(x);
%      end
% 
%   See TESTEOM.M for an example.
%
%   METHOD 2: JatUniverse models and state:
%   If you decide to use the predefined Jat Force Model 
%   cases with the JAT integrator you will need to specify which
%   force model case you would like to use. Currently there are 
%   two options: 'JatUniverseJGM2' and 'JatUniverseJGM3'.
%
%   JatUniverseJGM2 is made up of the following:
%   Sun & Moon Perturbations
%   NRL Atmospheric Drag
%   Solar Radiation Pressure Forces
%   A 20X20 JGM2 Harmonic Gravity Model
%
%   JatUniverseJGM3 is made up of the following:
%   Sun & Moon Perturbations
%   NRL Atmospheric Drag
%   Solar Radiation Pressure Forces
%   A 20X20 JGM3 Harmonic Gravity Model
%
%   To use these force model cases respectively, you need to type:
%
%      jatRK8('JatUniverseJGM2',...
%      jatRK8('JatUniverseJGM3',...
%
%   You must also define the state in an array depending on the 
%   force models being used. If using the JAT predefined cases
%   you must use the following format:
%  
%      State = [x0 y0 z0 x0dot y0dot z0dot] (meters and seconds)
%
%   Next define the time output, you can choose to define only 
%   the starting and ending values (see time_1 below) in which case the 
%   output frequency will be determined by the stepSize, either specified 
%   by the user, or computed internally to correspond to 0.5 deg true 
%   anomaly. You can also choose to dictate a different output frequency 
%   by creating an array of times you would like to output (see time_2 
%   below) or by using the colon operator to create an array of time 
%   (see time_3 below). 
%
%   time_1 = [ti, tf]
%   time_2 = [ti, t2, t5, t7, t8, ..., tf]
%   time_3 = [ti:dt:tf]  Where dt is the stepsize of output
%
%   Finally, you will need to set the adaptor options using the 
%   option setting function.
%
%   setJatRK8Options('NAME',Value, 'NAME2', Value2,...)
% 
%   OPTIONAL INPUTS: 
%    
%     stepSize = This is the stepSize for the integrator in seconds.
%                **NOTE: If stepSize is not specified, the default is 
%                computed to correspond to 0.5 deg true anomalay with the 
%                default gravitational parameter corresponding to Earth.
%
%     coefD    = This is the coefficient of drag for the spacecraft in orbit.
%
%     cr       = Coefficient of reflectivity for the spacecraft in orbit.
%
%     mass     = The mass of the spacecraft in orbit in kilograms.
% 
%     cArea    = The cross-sectional area of the spacecraft in orbit in m^2.
%
%     mjd_utc  = The Modified Julian Date in UTC time.
%
%   If you are using the JAT predefined cases and you want to have
%   polar motion turned on and not set to zero then you will need to add
%   your own EOP.dat file to match the dates you are running your
%   simulation over. You will need to put your EOP.dat file into
%   ODTBX1\Jat\jat\spacetime\EOP.dat. Your file must be names EOP.dat you
%   may want to rename the EOP.dat file already contained in that directory
%   for future use so it is not written over by your new file. 
%
%   In the case of the JAT predefined cases all of the optional
%   inputs are required and as such default values will be used if they 
%   are not defined by the user. 
% 
%   Refer to matlabjatinterface.m for an example program using the JAT 
%   propagator through the adaptor. It can be found in the examples directory
%   within the installation directory Orbit_Determination_Toolbox.
%
%   You can also choose the degree and order for your Spherical
%   Harmonic Gravity Model. The default setting is a 20x20 field for
%   both JGM2 and JGM3. Up to a 69x69 size field is allowed. If you
%   would like a larger field you will need to supply your own .grv
%   file. For an example of the file format refer to
%   ODTBX1/Jat/jat/forces/earthGravity/JGM2.grv or you can refer to
%   JGM3.grv located in the same directory. If you would like to put
%   your own file in, rename the respective .grv file in the folder
%   and name you file either JGM2.grv (for JGM2 simulations) or
%   JGM3.grv (for JGM3 simulations). \\
%
%   keyword: JAT Adaptor, RK8
%   See also getJatRK8Options, setJatRK8Options, TESTEOM
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
%               		(MM/DD/YYYY)
%   Kathryn Bradley     02/23/2006      Original
%   Kathryn Bradley	    04/17/2006	    Added documentation and call to 
%                                       adaptorGetOptions
%   Allen Brown         02/10/2007      Renamed to jatRK8, updated
%                                       comments.
%   Sun Hur-Diaz        04/10/2009      Modifed check for options and how
%                                       the stepsize is determined
%   Sun Hur-Diaz        04/14/2010      Replaced call to
%                                       calculatetrueanomalystep by a call
%                                       to trueanomalystep that is defined
%                                       here only

if(nargin>3)
    if(isstruct(varargin{1}))
        Options = varargin{1}; 
    else
        error('The 4th argument must be a structure.')
    end
else
    Options=setJatRK8Options; % Create an empty structure
end

if isempty(getJatRK8Options(Options, 'stepSize'))
    disp('Calculating stepsize corresponding to 0.5 deg true anomaly.')
    if length(x0)>=6
        stepSize = trueanomalystep(x0);
    else
        error('Default stepsize cannot be determined - specify with Options')
    end
else
    stepSize = getJatRK8Options(Options, 'stepSize');
end

disp(['Integrator step size is ',num2str(stepSize),' sec'])

if ((strcmpi(eoms, 'JatUniverseJGM2') || strcmpi(eoms, 'JatUniverseJGM3')))
    
    coefD = getJatRK8Options(Options, 'coefD', 2.2);
    cr = getJatRK8Options(Options, 'cr', 1.2);
    mass = getJatRK8Options(Options, 'mass',1000);
    cArea = getJatRK8Options(Options, 'cArea', 20);
    mjd_utc = getJatRK8Options(Options, 'mjd_utc', 53157.5);
    order = getJatRK8Options(Options, 'JGMOrder',20);
    degree = getJatRK8Options(Options,'JGMDegree',20);
    if(order>degree)
        disp('Error you cannot have an order greater than your degree. Order will be set to equal your degree.')
        order = degree;
    end
    jatPath = whereIsJat();
    oldState = jat.matlabInterface.PropagatorAdaptor.RK8(eoms,time, x0, stepSize, coefD, cr, mass, cArea, mjd_utc, jatPath, order, degree);

else % Using user-specified function

    oldState = jat.matlabInterface.PropagatorAdaptor.RK8(eoms,time, x0, stepSize);

end

t = oldState(:,1);
y = oldState(:,2:end);

end

function newStep=trueanomalystep(x)
GM = 3.9860044*10^14; % m^3/s^2 Gravitational Constant for Earth
delta_ta = .5; % deg
OE = kepel(x,[],GM);
n = sqrt(GM/OE.sma^3);
OE1 = OE;
OE1.tran = OE1.tran + delta_ta*pi/180;
newStep = mod((kepanom(OE1,'M') - kepanom(OE,'M')),2*pi)/n;
end


