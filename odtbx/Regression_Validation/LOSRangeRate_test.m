function failed = LOSRangeRate_test()
%
% LOSRangeRate_test Regression test for LOSRangeRate
% See also: LOSRangeRate.m
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
%   Ravi Mathur         08/28/2012      Extracted from LOSRangeRate.m

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%12s%17s%19s%17s%20s%20s\n','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (km/s)','Calculated (km/s)')
fprintf('%s\n\n',char(ones(1,105)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006'));
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);

t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [e1 e1 e1 e1 e1 e1 e1 e1 e1];
r2 = [e2 e2 e2 e2 e2 e2 e2 e2 e2];    
v1 = [e1 e1 e1 e2 e2 e2 e3 e3 e3];
v2 = [e1 e2 e3 e1 e2 e3 e1 e2 e3];

ExRangeRate = [0 1.41421022674001 0.707106781186547 -1.41421689802191 ...
    0 -0.707106781186547 -0.707108449010957 0.707105113370005 0];
RangeRate = LOSRangeRate(t, r1, v1, r2, v2, options);

fprintf('%4.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %16.6f %16.6f\n',...
    [r1; r2; v1; v2; ExRangeRate; RangeRate]);

passed = tol > max( abs( ExRangeRate - RangeRate ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end