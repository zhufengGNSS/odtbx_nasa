function failed = LOSRange_test()
%
% LOSRange_test Regression test for LOSRange
% See also: LOSRange.m
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
%   Ravi Mathur         08/28/2012      Extracted from LOSRange.m

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%20s%28s%25s%20s\n','Pos 1 (km)','Pos 2 (km)','R-Expected (km)','R-Calculated (km)')
fprintf('%s\n\n',char(ones(1,92)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006'));
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);

t=(1:9)*60*60;
e1 = [10000; 0; 0]; e2 = [0; 10000; 0]; e3 = [0; 0; 10000];
r1 = [e1 e1 e1 e2 e2 e2 e3 e3 e3];
r2 = [e1 e2 e3 e1 e2 e3 e1 e2 e3];    

ExRange = [0 sqrt(2e8) sqrt(2e8) sqrt(2e8) 0 sqrt(2e8) sqrt(2e8) sqrt(2e8) 0];
Range = LOSRange(t, r1, r2, options);

fprintf('%8.2f %8.2f %8.2f %10.2f %8.2f %8.2f %16.6f %16.6f\n',...
    [r1; r2; ExRange; Range]);

passed = tol > max( abs( ExRange - Range ) );
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
