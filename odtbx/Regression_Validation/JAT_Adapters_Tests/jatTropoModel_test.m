function failed = jatTropoModel_test(type)
% jatTropoModel_test Regression test for jatTropoModel
% See also: jatTropoModel.m
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%
%   Ravi Mathur         08/29/2012      Extracted from jatTropoModel.m

if strcmpi(type,'ValidationTest')

	failed = jatTropoModel_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = jatTropoModel_regression_test();

else
    disp('jatTropoModel_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = jatTropoModel_validation_test()
% jatTropoModel_validation_test - jatTropoModel Validation Tests
% 
% This function validates the jatTropoModel function against the test case
% contained in TropModel.java.
%


%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/04/2008	 	Original
%   Ravi Mathur            08/29/2012      Extracted from jatTropoModel.m

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-10;

c = 299792458.0;
L1_freq = 1575.42E+06;
lambda = c / L1_freq;
els = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0]*pi/180;
p = 938.0;
T = 286.0;
fh = 0.73;
rho = 1000000.0;

[dRho, dE] = jatTropoModel(rho, els, lambda, p, T, fh );

disp(['Elevation (deg) ',...
      ' Range Err (m)  ',...
      'El Err (arcsec) '])

disp(['----------------',...
      '----------------',...
      '----------------'])
disp([els'*180/pi,dRho',dE'])


if exist('jatTropoModelJatOutputs6_08.txt') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	% Load the JAT Tropo Model outputs
	load jatTropoModelJatOutputs6_08.txt

	format short e

	Diff = jatTropoModelJatOutputs6_08(:,2:3) - [dRho',dE']

	format short g

	if any(any( abs(Diff) > tol )) |  any(any(isnan(Diff)))
		failed = 1;
	end

else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function failed = jatTropoModel_regression_test()
% jatTropoModel_regression_test - jatTropoModel Validation Tests
% 

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/05/2008	 	Original
%   Ravi Mathur            08/29/2012      Extracted from jatTropoModel.m

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-10;

c = 299792458.0;
L1_freq = 1575.42E+06;
lambda = c / L1_freq;
els = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0]*pi/180;
p = 938.0;
T = 286.0;
fh = 0.73;
rho = 1000000.0;

[dRho, dE] = jatTropoModel(rho, els, lambda, p, T, fh );

disp(['Elevation (deg) ',...
      ' Range Err (m)  ',...
      'El Err (arcsec) '])

disp(['----------------',...
      '----------------',...
      '----------------'])
disp([els'*180/pi,dRho',dE'])


if exist('jatTropoModelRegressionData6_08.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load the JAT Tropo Model outputs
	load jatTropoModelRegressionData6_08

	format short e

	Diff = PreviousValues - [dRho',dE']

	format short g

	if any(any( abs(Diff) > tol )) |  any(any(isnan(Diff)))
		failed = 1;
	end

else
	failed = 1;
end


end