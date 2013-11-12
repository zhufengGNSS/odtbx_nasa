function [y,H,R,yabs,Habs,Rabs] = gpsmeas_test1()
% This demonstrates the use of the gpsmeas measurement model with several
% inputs times and several gps satellites.  It is useful for producing a
% repeatable data set via a batch-style call.
% The example holds the position of the spacecraft stationary and rotates
% the Earth and GPS satellites over time.
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

%           state_j2k        the J2000 position and velocity
%                    .t      the time of day in seconds
%                    .mjd    the modified julian day (integer)
%                    .p      the mean of J2000 position in km
%                    .v      the mean of J2000 velocity in km/s

x1 = [6.661728208278700e+003    6.661278931448849e+003    6.660828815452385e+003    6.660377860346100e+003
    4.974733893323415e+003    4.971291888137325e+003    4.967849256683786e+003    4.964405999396554e+003
    1.344395137958589e+004    1.346747804712036e+004    1.349100301806504e+004    1.351452628945648e+004
    -4.488572276416573e-001   -4.496964226153426e-001   -4.505355608480570e-001   -4.513746422337851e-001
    -3.441691907798407e-001   -3.442318392101571e-001   -3.442944442685567e-001   -3.443570059471385e-001
    2.352751484151581e+000    2.352581973349160e+000    2.352412166194112e+000    2.352242062707958e+000];

x1(4:6,:) = zeros(3,size(x1,2));

epoch  = datenum('Jan 1 2006');

v1.state_j2k.mjd    = epoch;
v1.state_j2k.tSim   = [0 10 20 30];
v1.state_j2k.p      = x1(1:3,:);
v1.state_j2k.v      = x1(4:6,:);
v1.clk.b            = [0 1e-5 2e-5 3e-5];
v1.clock.b          = v1.clk.b;
v1.clock.tSim       = v1.state_j2k.tSim;

param.measrate      = 10;
param.noise.sgs     = 0;
param.noise.sign    = 1e-3;
param.noise.bias    = 1e-3;

measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useiono',false);
measOptions = setOdtbxOptions(measOptions,'usetropo',false);

tPos  = v1.state_j2k.tSim;

% This demonstrates batch measurement generation
[y,H,R] = gpsmeas(tPos, x1, measOptions)
[yabs,Habs,Rabs] = gpsmeasAbs(tPos, x1, measOptions,[])


