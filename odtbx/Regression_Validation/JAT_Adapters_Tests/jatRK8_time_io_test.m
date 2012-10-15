function fail = jatRK8_time_io_test(doplot)
% Regression test and script to replicate and then test for regressions
% with the jatRK8 JAT Adaptor's handling of time inputs.
% (Modified from a script to replicate Mantis Issue 147.)
%
% The intent of this function is to show how the code failed (with the
% first jatRK8 invocation) before the fix, AND now to show that jatRK8
% operates as designed.  
%
% This function tests the call to the jatRK8 propagator with several 
% different time specifications.  The results should match even though the
% reported output times are different.  This shows that 
% propagator accuracy is maintained and separable from the time 
% specification format.
%
% Allen Brown, May and June 2009, ODTBX Mantis issue 147.
%
% Inputs:
% doplot (Optional) if set and not empty then comparison plots are shown
%                   (don't use with regression suite)
%
% Outputs:
% fail: 0 for success, 1 for failure (to meet the regression test
% interface)
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

fail = 0;


% This first section replicated the problem with 
% ODTBX /trunk r633 and JAT /branches/ODTBX r94.
x = 8881.95135138262; % km
y = 9430.77519139708; % km
z = 4162.93561538469; % km
xdot = -4.71684695624137; % km/s
ydot = -2.37093913899615; % km/s
zdot = 3.34006991068724; % km/s

R0 = 1000.0*[x y z xdot ydot zdot];

tspan = [0.0 86400.0];
[t_two, y_two] = jatRK8('JatUniverseJGM2', tspan, R0);
%
% End first section, JAT exception occurred in the jatRK8 call, above.
% After the fix, this call operates as expected.  Output points are
% determined by JAT.  (In this case, 3359 total output points.)
%

% Testing three specified output points:
tspan = [0.0 86400.0/2 86400.0];
[t_three, y_three] = jatRK8('JatUniverseJGM2', tspan, R0);

% Testing 101 specified output points:
tspan = 0.0:(86400.0/100):86400.0;
[t_101, y_101] = jatRK8('JatUniverseJGM2', tspan, R0);

% Testing 101 specified output points, the same time span,
% but changing the default step size (roughly doubling it, from 25.7 to 50)
% (do this just like in jatRK8.m).  The propagator accuracy shouldn't
% really differ too much with these EOMs and these time step sizes.
Options = setJatRK8Options;
Options = setJatRK8Options(Options, 'coefD', 2.2);
Options = setJatRK8Options(Options, 'cr', 1.2);
Options = setJatRK8Options(Options, 'mass',1000);
Options = setJatRK8Options(Options, 'cArea', 20);
Options = setJatRK8Options(Options, 'mjd_utc', 53157.5);
Options = setJatRK8Options(Options, 'JGMOrder',20);
Options = setJatRK8Options(Options,'JGMDegree',20);

Options = setJatRK8Options(Options, 'stepSize',50);

[t_101ss, y_101ss] = jatRK8('JatUniverseJGM2', tspan, R0, Options);

% comparison plots
if exist('doplot','var') && ~isempty(doplot)
    % All the data points should lie right on top of each other.  This means
    % that the propagator is "doing the right thing" no matter how you specify
    % time.  Also, it shows that manually setting the stepSize can be done
    % independent of the time specification.
    figure;
    hold on;
    plot(t_two, sqrt(sum(y_two(:,1:3).*y_two(:,1:3),2)), 'kd');
    plot(t_three, sqrt(sum(y_three(:,1:3).*y_three(:,1:3),2)), 'gx-');
    plot(t_101, sqrt(sum(y_101(:,1:3).*y_101(:,1:3),2)), 'ro');
    plot(t_101ss, sqrt(sum(y_101ss(:,1:3).*y_101ss(:,1:3),2)), 'b.');
    legend('two times', 'three times', 'def stepSize','double stepSize');
    hold off;
    title('Position Mag Co-Plot');
    xlabel('Sim time (sec)');
    ylabel('Position (m)');

    figure;
    hold on;
    plot(t_two, sqrt(sum(y_two(:,4:6).*y_two(:,4:6),2)), 'kd');
    plot(t_three, sqrt(sum(y_three(:,4:6).*y_three(:,4:6),2)), 'gx-');
    plot(t_101, sqrt(sum(y_101(:,4:6).*y_101(:,4:6),2)), 'ro');
    plot(t_101ss, sqrt(sum(y_101ss(:,4:6).*y_101ss(:,4:6),2)), 'b.');
    legend('two times', 'three times', 'def stepSize','double stepSize');
    hold off;
    title('Velocity Mag Co-Plot');
    xlabel('Sim time (sec)');
    ylabel('Velocity (m/s)');
end

% Perform interpolation and differencing
% Interpolate the finest-grained data (t_two, y_two) with the times
% of the other data

% tolerances:
postol = 20; % just above the max difference for all data sets
veltol = .1; % just above the max difference for all data sets
y3idiff = zeros(size(y_three));
for i = 1:size(y_two,2)
    y3idiff(:,i) = y_three(:,i) - interp1(t_two, y_two(:,i), t_three,'cubic');
    if i < 4
        tol = postol;
    else
        tol = veltol;
    end
    if any(abs(y3idiff(:,i)) > tol)
        fail = 1;
        disp(['jatRK8_time_io_test: Failed tolerance check for y_three dataset, index' i]);
    end
end

y101idiff = zeros(size(y_101));
for i = 1:size(y_two,2)
    y101idiff(:,i) = y_101(:,i) - interp1(t_two, y_two(:,i), t_101,'cubic');
    if i < 4
        tol = postol;
    else
        tol = veltol;
    end
    if any(abs(y101idiff(:,i)) > tol)
        fail = 1;
        disp(['jatRK8_time_io_test: Failed tolerance check for y_101 dataset, index' i]);
    end
end

y101ssidiff = zeros(size(y_101ss));
for i = 1:size(y_two,2)
    y101ssidiff(:,i) = y_101ss(:,i) - interp1(t_two, y_two(:,i), t_101ss,'cubic');
    if i < 4
        tol = postol;
    else
        tol = veltol;
    end
    if any(abs(y101ssidiff(:,i)) > tol)
        fail = 1;
        disp(['jatRK8_time_io_test: Failed tolerance check for y_101ss dataset, index' i]);
    end
end

%
% This section tests the user-defined EOM input method, instead of using
% the JatUniverseJGM2 method.
%
optj=setJatRK8Options('stepSize',.01);
[ta,ya] = jatRK8('testrk8',[0 10],10,optj); % let JAT select and return lots of time points
[tb,yb] = jatRK8('testrk8',[0 1 5 10],10,optj); % specify some time points

ybadiff = yb - interp1(ta, ya, tb, 'cubic');
if any(abs(ybadiff) > 1e-9)
        fail = 1;
        disp('jatRK8_time_io_test: Failed tolerance check for ybadiff dataset');
end