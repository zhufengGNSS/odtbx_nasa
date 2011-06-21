% Script for performing timing on gpsmeas.
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

% the user must check some things before running:
disp('HEY! Before we begin, you did set the MATLAB process to have');
disp('single-CPU affinity, right? (for more accurate results)');
disp('...Do that now before continuing.  Hit the spacebar to continue.');
pause;

% reset the profiler to recognize any CPU affinity changes
profile -timer real
profile -timer cpu

% init the workspace
close all;
clear; % just variables

% call once without profiling to ensure Java classes are loaded:
GPSMEAS_PROFILE_OVERRIDE = 1; % we only need to call the loop once, use an override
gpsmeas_sequential_fortiming
clear GPSMEAS_PROFILE_OVERRIDE; % done with the override, run normally

% now, start the profiler and time it
disp('Starting sequential timing...')
profile on
gpsmeas_sequential_fortiming
profile viewer
profile off
disp('...done with sequential timing.')

stats = profile('info');
dirname = 'results_profile_gpsmeas_seq';
disp(sprintf('Profile results saved in %s',dirname));
profsave(stats,dirname);

disp('Reminder: Don''t forget to change back any altered CPU affinity settings...');
