function failed = lightTimeCorrection_test()
%
% lightTimeCorrection_test Regression test for lightTimeCorrection
% See also: lightTimeCorrection.m
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
%   Ravi Mathur         08/28/2012      Extracted from lightTimeCorrection.m

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%8s%12s%17s%19s%17s%19s%18s\n','t (sec)','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (sec)','Calculated (sec)')
fprintf('%s\n\n',char(ones(1,110)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006'));
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);
options = setOdtbxOptions(options, 'useRange', true );
options = setOdtbxOptions(options, 'useRangeRate', true );
options = setOdtbxOptions(options, 'useDoppler', false );
t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [repmat(e1,1,3) repmat(e2,1,3) repmat(e3,1,3)]*9;
r2 = repmat([e1 e2 e3],1,3)*5;    
v1 = repmat(e1,1,9);
v2 = repmat([e1 e2 e3],1,3);
x1 = [r1;v1]; x2 = [r2;v2];

Ex_t2 = [3599.99998665744 7199.99996565747 10799.9999656575 14399.9999656575 ...
    17999.9999866574 21599.9999656575 25199.9999656575 28799.9999656575 32399.9999866574];
t2_lt = lightTimeCorrection(t, x1, t, x2, options,-1,1e-10,10);

fprintf('%6i %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %16.6f %16.6f\n',...
    [t; r1; r2; v1; v2; Ex_t2; t2_lt]);

failed = tol < max( abs( Ex_t2 - t2_lt ) );
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
