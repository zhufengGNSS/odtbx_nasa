function [f,failed]=jatWorldPropagatorRKF78_regression
% Regression test for interpolation in jatWorldPropagatorRKF78
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
% Emergent Space Technologies
% Modified 

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Sun Hur-Diaz           			Original
%   Stephen D Metcalfe  06/09/2010   	Modified to document which 
%							test case calls which interpolator
%							and how to update the test data.
%   Stephen D Metcalfe  09/29/2010      Modified to explicitly reference 
%                                       test data in the DataFiles folder.
%                                       Added timestamp and name to test 
%                                       data to identify who generated it 
%                                       and when. Added print of when test 
%                                       data was generated.
%                                       Added pass/fail printout.

failed = zeros(2,1);
load DataFiles/jatWorldPropagatorRKF78_regression_data
fprintf('Test data generated on %s\n', GENTIME);
odeopts = odeset('RelTol',1e-3);

fprintf('\nCase 1 - using HermiteInterpolator\n');
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', datenum('01 Jan 2010 00:00:00'));
jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM2');
jatWorld = createJATWorld(jOptions);
[tt,xprop] = jatWorldPropagatorRKF78([],0:10:60,[15000 0 0 0 4 1]'.*1000,odeopts,[],[],[],jatWorld);
if max(abs(tt-tt1) > 1e-10) || max(max(abs(xprop - xprop1))) > 1e-10;
   failed(1) = 1;
   % uncomment the following 4 lines to update the regression data
   % WHOGEN =                                  %Insert your name here
   % GENTIME = datestr(now);
   % xprop1 = xprop;
   % save DataFiles/jatWorldPropagatorRKF78_regression_data WHOGEN GENTIME tt1 tt2 xprop1 xprop2
   fprintf('Failed\n');
else
   fprintf('Passed\n');
end

fprintf('\nCase 2 - using LagrangianInterpolator\n');
jOptions = odtbxOptions('force');
jOptions = setOdtbxOptions(jOptions, 'epoch', datenum('01 Jan 2010 00:00:00'));
jOptions = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM2');
jatWorld = createJATWorld(jOptions);
[tt,xprop] = jatWorldPropagatorRKF78([],0:10:60,[7000 0 0 4 0 0]'.*1000,odeopts,[],[],[],jatWorld);
if max(abs(tt-tt2) > 1e-10) || max(max(abs(xprop - xprop2))) > 1e-10;
   failed(2) = 1;
   % uncomment the following 3 lines to update the regression data
   % WHOGEN =                                  %Insert your name here
   % GENTIME = datestr(now);
   % xprop2 = xprop;
   % save DataFiles/jatWorldPropagatorRKF78_regression_data WHOGEN GENTIME tt1 tt2 xprop1 xprop2
   fprintf('Failed\n');
else
   fprintf('Passed\n');
end

f=any(failed);


