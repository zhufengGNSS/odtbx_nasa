% Sequential Estimator Demo
% 
% This demo script calls various examples that show the use of the
% sequential estimator function estseq.  
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

% Sun Hur-Diaz
% Emergent Space Technologies, Inc.
% Greenbelt, Md
% July 30, 2009
%
echo on
%
% ************************************************************************
% This demo script calls various examples that show the use of the
% sequential estimator function estseq.  
% ************************************************************************
%
% EXAMPLE 1 shows the use of estseq with user-supplied dynamics model
% and measurement model functions.  In addition, some of the diagnostic
% tools like estval, plot_ominusc, and varpiles are demonstrated through
% the function plot_results.  This example solves the same problem as
% Example 1 of the batch_estimator_demo.
%
% In this example, the dynamics model 'r2bp' is the restricted 2-body 
% problem:
%
pause % Hit RET to continue
%
type r2bp
%
% There are 6 states which are the 3D inertial position and velocity.
% Note that the Jacobian and the process noise are specified in this model.
% If the Jacobian output by the dynamics function had been an empty 
% matrix, then estseq would have computed a numerical Jacobian.  
%
pause % Hit RET to continue
%
% The measurement model 'rrdot3D_1way' is the range and range-rate from a 
% single fixed ground station at the Earth's equator.  
%
pause % Hit RET to continue
%
type rrdot3D_1way
%
% Again, the Jacobian is specified.  If it were output as an empty matrix,
% estseq would have computed the numerical Jacobian.
%
pause % Hit RET to continue
%
% Before we can execute the estseq command we need to specify certain
% inputs, parameters, and options.  Note that in this example, the true and
% the estimator models are identical except that the estimator measurement
% noise is 3 times that of the truth.  Furthermore, the initial covariance 
% is specified.  The orbit is a 500km altitude LEO polar orbit and the time 
% span is only 5 minutes at 10 sec measurement intervals.
%
pause % Hit RET to continue
%
sequential_example_1
%
echo on
%
% It is interesting to compare the results of this example with the
% similar one in Example 1 of the batch_estimator_demo.  Note that the
% estseq plots of the Monte Carlo simulations have 3 extra time points 
% between measurement times.  This is because of the parameter 'refint'
% which specifies the number of  points between measurements to propagate 
% and is currently hard-coded at 3.  Becasue estseq sequentially processes
% the measurements, a valid comparison with the batch filter can really 
% only be made at the end of the time span.  We find that the covariances
% are fairly comparable in both cases at the end time.
%
% ************************************************************************
%
% EXAMPLE 2 solves the same problem as Example 1 but using ODTBX-supplied 
% dynamics model and measurement model functions.  The dynamics model is 
% based on jatForces, and the measurement model is based on gsmeas.  It 
% also shows how the integrator options can be set instead of using the 
% default settings as in Example 1.
%
pause % Hit RET to continue
%
sequential_example_2
%
echo off

