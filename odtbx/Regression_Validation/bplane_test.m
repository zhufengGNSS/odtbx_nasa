function failed = bplane_test()
% Test function for kepel.m.
%
% Tests one elliptical trajectory about the Earth with different argument
% presentation.  Didn't test passing in GM as an argument.
%
% From Vallad, David A., " Fundamentals of Astrodynamics and Applications",
% Third Ed., 2007, Example problem 2-5, pages 122-124.
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

failed = 0;
tol = 1e-9;

% Define test variables
Bt = nan(3,1);
Br = nan(3,1);
LTOF = nan(3,1);

% Test #1

% Initial spacecraft state (Mars-Centered Inertial)
R = [-64961.204483358 -133190.207105847 0]';
V = [1.37629832487732 2.39127882923125 0]';

% Mars orbit normal vector
N = [0 0 1]';

% Mars gravitational parameter
mu_mars = 42828;

% Calculate B-Plane parameters for current state
[Bt(1),Br(1),LTOF(1)] = bplane([R;V],N,mu_mars);

% Test #2

% Calculate B-Plane parameters for current state
[Bt(2),Br(2),LTOF(2)] = bplane([R;V],[],mu_mars);

% Test #3

% Calculate B-Plane parameters for current state
[Bt(3),Br(3),LTOF(3)] = bplane([R;V],[],[]);

% Uncomment to generate test data
% Bt_test = Bt;
% Br_test = Br;
% LTOF_test = LTOF;
% save DataFiles/bplane_test_data Bt_test Br_test LTOF_test

load DataFiles/bplane_test_data

Bt_err = Bt_test-Bt;
Br_err = Br_test-Br;
LTOF_err = LTOF_test-LTOF;

if any(Bt_err>tol) || any(Br_err>tol) || any(LTOF_err>tol)
    failed = 1;
    fprintf('B-Plane calculation failed!\n')
end

end