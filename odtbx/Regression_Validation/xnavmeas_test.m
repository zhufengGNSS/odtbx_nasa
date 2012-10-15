function failed = xnavmeas_test
% xnavmeas_test Performs regression testing on xnavmeas function
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
%
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur             08/27/2012      Extract regression test from
%                                           original xnavmeas.m

disp(' ')
disp('Performing Test....')
disp(' ')

% set up state vector and time
Ra = 6378137.; % m
mu = 3.986005e14;
alt_0 = 35786*1e3; % GEO alt in meters
Pis = Ra + alt_0;
Vis = sqrt(mu/Pis);
X(1:3,1) = [Pis 0 0]';
X(4:6,1) = [0 Vis 0]';
t=0;

tol = 1e-7;
xnav.useRangeRate = true;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006')); %UTC
options = setOdtbxOptions(options, 'xnav', xnav );

ExMeas = [4419540.83700585
          2832.20891721189
          4004178.73744371
         -2774.21508508712
          16436890.7878644
         -2595.71020850086
%            25595548.974467
%           1284.05142786188
];
y = xnavmeas(t,X,options);

fprintf('%s\n',char(ones(1,37)*'-'));
disp('Calculated XNAV Measurements by time:')
fprintf('%s\n',char(ones(1,37)*'-'));
fprintf('%16s %18s\n','Expected','Observed');
fprintf('%s\n',char(ones(1,37)*'-'));
for n=1:length(y)
    fprintf('%18.6f %18.6f\n',ExMeas(n),y(n))
end
fprintf('%s\n',char(ones(1,37)*'-'));

failed1 = any(tol < max( abs( ExMeas - y ) ));
if failed1
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end

% Regression Test of the Sched capability
Sched = [...
    1 10 20
    3 30 40];
options = setOdtbxOptions(options,'Schedule',Sched);
tt = 0:10:40;
XX = repmat(X,1,length(tt));
y = xnavmeas(tt,XX,options);
ytest =  [NaN 4419540.83700585 4419540.83700585                 NaN                NaN
          NaN 2832.20891721189 2832.20891721189                 NaN                NaN
          NaN              NaN              NaN                 NaN                NaN
          NaN              NaN              NaN                 NaN                NaN
          NaN              NaN              NaN    16436890.7878644   16436890.7878644
          NaN              NaN              NaN   -2595.71020850086  -2595.71020850086
%           NaN              NaN              NaN                 NaN                NaN
%           NaN              NaN              NaN                 NaN                NaN
          ];
failed2 = any(tol < max( abs( ytest - y ) ));
if failed2
    disp(' ')
    disp('Sched Test Failed!')
else
    disp(' ')
    disp('Sched Test Passed.')
end
          
failed = any([failed1 failed2]);

end
