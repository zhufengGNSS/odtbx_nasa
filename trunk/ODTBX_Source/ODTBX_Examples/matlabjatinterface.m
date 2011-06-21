% Example using the JAT RK8 Integrator and Matlab.
%
% This is a test file for the MATLAB JAT Interface using jatRK8. Refer to
% it for an example of how to access the JAT propagator.
%
%    keyword: Example Programs, JAT Adaptor,
%    See also setJatRK8Options, getJatRK8Options, jatRK8, TESTEOM
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
%   Kathryn Bradley     02/09/2006   	Original
%   Kathryn Bradley	04/17/2006		Added documentation and adaptorSetOptions 
%   Allen Brown         02/11/2009      Brought back to work with jatRK8.
%   Allen Brown         03/24/2009      Updated documentation.
      
close all
clear all
%clc

%For Standard two state example
%x0 = [1, 1];

% For Jat Universe Case (m/s) x0 = [x, y, z, xdot, ydot, zdot];
x0 = [-4453783.586, -5038203.756, -426384.456, 3831.888, -2887.221, -6018.232];

% Can choose this time option if the user would like to output times not 
% in sequential order.
%time = [0 .3 .6 .11 .5 .72 .9 2.5 4.7 8.2 8.7 9.2 10]; 
% Can choose this time option if the user would like to change the frequency 
% of output from the default.
%time = [0:1:10];

% When this time Option is chosen the frequency will be determined by the
% stepSize.
time= [0,10];

% Here is where the Options need to be set by the user, if any of these
% options are not set and are required by a force model the default will be
% used.
Options = setJatRK8Options('stepSize', 1, 'coefd', 2.2, 'cr', 1.2, 'mass', 1000, 'cArea', 20, 'mjd_utc', 53157.5, 'JGMOrder', 20, 'JGMDegree', 20);
%Options = setJatRK8Options('forces','jat.forces.JGM2,jat.forces.Sun,Moon.m,jat.forces.NRL,SRP.m');

%For User defined example, where testEOM is the name of the user defined
%EOM function.
%[t,y] = jatRK8('testEOM',time, x0, Options);

%For Jat Universe Example, see jatRK8
[t,y] = jatRK8('JatUniverseJGM2',time, x0, Options);

% Jat Universe Output, You can choose to output however you would like.
% fid = fopen('JGM_MATLAB_1day.txt', 'w');
% fprintf(fid,'Time (s)\t     X\t            Y\t             Z\t           Vx\t         Vy\t          Vz\n--------    ------------     ------------   ------------   ------------  -----------   -----------\n');
% for i=1:length(t);
% fprintf(fid,'%3.4f %17.13f %16.13f %13.13f %12.13f %13.13f %13.13f\n', t(i), y(i,1), y(i,2), y(i,3), y(i,4), y(i,5), y(i,6));
% end
% fclose(fid);

fprintf('Time (s)\t     X\t            Y\t             Z\t           Vx\t         Vy\t          Vz\n--------    ------------     ------------   ------------   ------------  -----------   -----------\n')
for i=1:length(t);
fprintf('%3.4f %17.6f %16.6f %13.6f %12.6f %13.6f %13.6f\n', t(i), y(i,1), y(i,2), y(i,3), y(i,4), y(i,5), y(i,6))
end

% Two Body Output, You can choose to output however you would like.
% fid = fopen('MATLAB_twoBody.txt', 'w');
% fprintf(fid,'  Time (s)\t       X1\t            X2\n----------     ----------       ----------\n');
% for i=1:length(t);
%     fprintf(fid,'%3.4f %17.16f %16.16f\n', t(i), y(i,1), y(i,2));
% end
% fclose(fid);

% fprintf('  Time (s)\t       X1\t            X2\n----------     ----------       ----------\n');
% for i=1:length(t);
%     fprintf('%3.4f %17.6f %16.6f\n', t(i), y(i,1), y(i,2));
% end




