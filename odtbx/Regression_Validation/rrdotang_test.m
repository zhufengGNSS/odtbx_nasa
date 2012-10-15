function failed = rrdotang_test()
%
% rrdotang_test Regression test for rrdotang
% See also: rrdotang.m
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%
%   Ravi Mathur         08/29/2012      Extracted from rrdotang.m

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%12s%17s%19s%17s%22s%23s\n','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (km km/s)','Calculated (km km/s)')
fprintf('%s\n\n',char(ones(1,110)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006')); %UTC
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);
options = setOdtbxOptions(options, 'useRange', true );
options = setOdtbxOptions(options, 'useRangeRate', true );
options = setOdtbxOptions(options, 'useDoppler', false );
options = setOdtbxOptions(options, 'useUnit', false );
options = setOdtbxOptions(options, 'useAngles', true );
t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [repmat(e1,1,3) repmat(e2,1,3) repmat(e3,1,3)]*9;
r2 = repmat([e1 e2 e3],1,3)*5;    
v1 = repmat(e1,1,9);
v2 = repmat([e1 e2 e3],1,3);
x1 = [r1;v1]; x2 = [r2;v2];

ExMeas = [               4                         0                         0                         0
           10.295630140987          1.35979800452198        -0.507098504392337                         0
           10.295630140987          1.35979800452198                         0        -0.507098504392337
           10.295630140987                         0          2.07789483118723                         0
                         4         -1.00000333565208           1.5707963267949                         0
           10.295630140987         0.485642144472135           1.5707963267949        -0.507098504392337
           10.295630140987                         0          3.14159265358979          1.06369782240256
           10.295630140987         0.485642144472135          -1.5707963267949          1.06369782240256
                         4         -1.00000333565208                         0                         0]';
y = rrdotang(t,x1,x2,options);

fprintf('%4.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %11.6f %9.6f %11.6f %9.6f\n',...
    [r1; r2; v1; v2; ExMeas; y]);

failed = any(tol < max( abs( ExMeas - y ) ));
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
