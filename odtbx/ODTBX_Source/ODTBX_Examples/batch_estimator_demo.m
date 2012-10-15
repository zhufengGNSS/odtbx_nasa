% Batch Estimator Demo
% 
% This demo script calls various examples that show the use of the batch
% estimator function estbat.  
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
% July 24, 2009
%
echo on
%
% ************************************************************************
% This demo script calls various examples that show the use of the batch
% estimator function estbat.  
% ************************************************************************
%
% EXAMPLE 1 shows the use of estbat with user-supplied dynamics model
% and measurement model functions.  In addition, some of the diagnostic
% tools like estval, plot_results, and varpiles are demonstrated.
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
% matrix, then estbat would have computed a numerical Jacobian.  Note that 
% the process noise is used in the computation of the true covariance only.
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
% estbat would have computed the numerical Jacobian.
%
pause % Hit RET to continue
%
% Before we can execute the estbat command we need to specify certain
% inputs, parameters, and options.  Note that in this example, the true and
% the estimator models are identical except that the estimator measurement
% noise is 3 times that of the truth.  Furthermore, the initial covariance 
% is specified.  The orbit is a 500km altitude LEO polar orbit and the time 
% span is only 5 minutes at 10 sec measurement intervals.  For such a short 
% data span compared to the orbit period (95 min.), there is not enough 
% observability from the measurements.  The initial covariance, however, 
% helps to keep the normal matrix of the batch filter nonsingular.
%
pause % Hit RET to continue
%
batch_example_1
%
echo on
%
% ************************************************************************
%
% EXAMPLE 2 solves the same problem as Example 1 but using ODTBX-supplied 
% dynamics model and measurement model functions.  The dynamics model is 
% based on jatForces, and the measurement model is based on gsmeas.  It 
% also shows how the integrator options can be set instead of using the 
% default settings as in Example 1.  (Note this examples takes about 80
% times as long as Example 1 to complete.)
%
pause % Hit RET to continue
%
batch_example_2
%
echo on
%
% ************************************************************************
%
% EXAMPLE 3 is similar to Examples 1 and 2 but the time span is extended to
% 60,000 seconds to process more data.  The time interval is still 10 
% seconds.  To avoid out-of-memory problems arising from such a large time
% vector, the process noise is set to zero.  The dynamics model is given
% by the function jatForces_km_noQ which is a wrapper function around the 
% ODTBX-supplied function jatForces_km to force the process noise to be 
% zero.  (The memory issue may also be avoided, even with nonzero process 
% noise, by specifying the time vector to be only those times when 
% measurements are available.) 
%      The measurement model is given by rrdot3D_1way as in Example 1.  
% With the additional measurement data, the performance of batch filter is
% improved.
%
pause % Hit RET to continue
%
batch_example_3
%
echo off

