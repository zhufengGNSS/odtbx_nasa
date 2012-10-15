function failed = rrdotlt_test()
%
% rrdotlt_test Regression test for rrdotlt
% See also: rrdotlt.m
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
%   Ravi Mathur         08/28/2012      Extracted from rrdotlt.m

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
t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [repmat(e1,1,3) repmat(e2,1,3) repmat(e3,1,3)]*9;
r2 = repmat([e1 e2 e3],1,3)*5;    
v1 = repmat(e1,1,9);
v2 = repmat([e1 e2 e3],1,3);
x1 = [r1;v1]; x2 = [r2;v2];

ExMeas = [3.99999998004317 -3.99137830726775e-009; 10.295630124683 1.35979800279089; ...
    10.295630181257 1.35979801189009; 10.2956301806334 7.92926180965375e-009; ...
    4 -1.00000334028493; 10.2956301013406 0.485642143045818; ...
    10.295630100717 -8.05399233005487e-009; 10.295630157291 0.485642138321872; ...
    4.00000001995683 -1.00000332203861]';
y = rrdotlt(t,x1,t,x2,options);

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
