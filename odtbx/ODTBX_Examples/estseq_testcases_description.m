%% estseq.m Self-Test Cases Descriptions
%
%% Case 1
% This is a scalar constant state case.  Both the true and the assumed
% dynamics are the same (i.e., S=1, C=0) and are given below:  
%
% $$ \begin{array}{ll} \dot{x} = \dot{\hat{x}} = 0, & df/dx = d\hat{f}/d\hat{x} = 0 \end{array}$$
%
% The process noise PSD of the truth and the assumed are
%
% $$ \begin{array}{ll} Q = 1, & \hat{Q} = 0 \end{array} $$
%
% Both the true and the assumed measurement models are the same and are given below:
%
% $$ \begin{array}{ll} Y = x, & H = 1 \\
%                      \hat{Y} = \hat{x}, & \hat{H} = 1\end{array}$$
% 
% The measurement noise covariance of the truth and the assumed are the
% same
%
% $$ \begin{array}{ll} R = 1, & \hat{R} = 1 \end{array} $$
%
% The initial state and covariance for the truth and assumed are
% 
% $$ \begin{array}{ll} x_o = 0, & P_o = 10\\
%                      \bar{x}_o = 0, & \bar{P}_o = 10 \end{array} $$
% 
% The measurement updates occur at t = [1:1:5]. There are 125 Monte Carlo 
% simulations.  In the Monte Carlo simulations, there are 31 interior time 
% steps for the integrator between measurement updates, and each 
% measurement update is iterated 3 times.
%
% Since the only difference between the true and the assumed models is in
% the process noise, we expect the variances due to _a priori_ error and
% the measurement noise to be the same between the truth and the assumed
% while the variance due to process noise of the truth to be higher than
% the assumed.  With no assumed process noise, we expect the total assumed
% error covariance to decrease with each measurement update.  On the other
% hand, the true error covarince increases over time since the gains based
% on the assumed decrease over time and not enough error reduction is
% performed to cancel the growth during the time update.
%% Case 2
% This is an 8-state case of simple 3-D kinematic equations of motion with
% the measurement bias (b) and its rate (r) included:
%
% $$ \vec{x} = [\begin{array}{cccccccc}x & y & z & u & v & w & b & r\end{array}]^T $$
%
% where the true dynamics and its true process noise are given by
%
% $$
% \frac{d}{dt}\left[\begin{array}{c}x\\y\\z\\u\\v\\w\\b\\r\end{array}\right
% ]
%                = \left[\begin{array}{cccccccc} 
%                    0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0    
%                    \end{array}\right]\vec{x} $$
% 
% $$
% Q = \left[\begin{array}{cccccccc} 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 1e-6 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 1e-6 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 1e-6 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 & 0 & 1e-6    
%                    \end{array}\right] $$
% 
% The true measurement noise model is given by
% 
% $$ Y = \left[\begin{array}{cccccccc}
%        1 & 0 & 0 & 0 & 0 & 0 & -1 & 0 \\
%        0 & 1 & 0 & 0 & 0 & 0 & -1 & 0 \\
%        0 & 0 & 1 & 0 & 0 & 0 & -1 & 0
%        \end{array}\right]\vec{x} $$
%
% $$ R = \left[\begin{array}{ccc}
%        1 & 0 & 0 \\
%        0 & 1 & 0 \\
%        0 & 0 & 1
%        \end{array}\right] $$
%
% In the estimator, the dynamics model and its process noise consist only
% of the solve-for states corresponding to the 3-D position and velocity.
%
% $$ \vec{x} = [\begin{array}{cccccccc}x & y & z & u & v & w\end{array}]^T $$
%
% $$
% \frac{d}{dt}\left[\begin{array}{c}x\\y\\z\\u\\v\\w\end{array}\right
% ]
%                = \left[\begin{array}{cccccc} 
%                    0 & 0 & 0 & 1 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 1 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 1 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0     
%                    \end{array}\right]\vec{x} $$
% 
% $$
% Q = \left[\begin{array}{cccccc} 
%                    0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 0 \\ 
%                    0 & 0 & 0 & 1e-2 & 0 & 0 \\ 
%                    0 & 0 & 0 & 0 & 1e-2 & 0 \\ 
%                    0 & 0 & 0 & 0 & 0 & 1e-2 
%                    \end{array}\right] $$
% 
% The estimator measurement model also consists of only the solve-for
% states:
% 
% $$ Y = \left[\begin{array}{cccccc}
%        1 & 0 & 0 & 0 & 0 & 0 \\
%        0 & 1 & 0 & 0 & 0 & 0 \\
%        0 & 0 & 1 & 0 & 0 & 0 
%        \end{array}\right]\vec{x} $$
%
% $$ R = \left[\begin{array}{ccc}
%        1e4 & 0 & 0 \\
%        0 & 1e4 & 0 \\
%        0 & 0 & 1e4
%        \end{array}\right] $$
%
% The partition matrices are
%
% $$ S = \left[\begin{array}{cccccccc}
%              1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
%              0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
%              0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
%              0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
%              0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
%              0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\end{array}\right] $$
%
% $$ C = \left[\begin{array}{cccccccc}
%              0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
%              0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\end{array}\right] $$
%
% The initial states and covariances are
% 
% $$ x_o = [\begin{array}{cccccccc}0 & 0 & 0 & 0 & 0 & 0 & 0 &
% 0\end{array}]^T$$
%
% $$ P_o = \left[\begin{array}{cc}
%                I_{6x6} & O_{6x2}\\
%                O_{2x6} & 2I_{2x2}\end{array}\right] $$
%
% $$ \bar{x}_o = [\begin{array}{cccccc}0 & 0 & 0 & 0 & 0 & 0\end{array}]^T$$
%
% $$ \bar{P}_o = 2I_{6x6} $$
%
% The measurement updates occur at t = [1:1:30]. There are 12 Monte Carlo 
% simulations.  In the Monte Carlo simulations, there are 3 interior time 
% steps for the integrator between measurement updates, and each 
% measurement update is iterated 3 times.
%
%% Case 3
% Case 3 is a Schmidt-Kalman filter implementation of Case 2.  Both the true
% and the estimator models are identical corresponding to the true model of
% Case 2 above.  
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

