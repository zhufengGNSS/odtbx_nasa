function failed = jatChargedParticleModel_test(type)

% jatChargedParticleModel_test Regression test for jatChargedParticleModel
% See also: jatChargedParticleModel.m
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
%   Ravi Mathur         08/29/2012      Extracted from jatChargedParticleModel.m

if strcmpi(type,'ValidationTest')

	failed = jatChargedParticleModel_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = jatChargedParticleModel_regression_test();

else
    disp('jatChargedParticleModel_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jalChargedParticleModel_validation_test - Validation Tests
%
% tests against output from directly from java as well as testing multiple satellites

function failed = jatChargedParticleModel_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/09/2008	 	Original
%   Ravi Mathur            08/29/2012      Extracted from jatChargedParticleModel.m

disp(' ')
disp('Performing Test....')
disp(' ')

% Case 1 - test against Jat
Ephem.satPos = [2.056020937918350   1.411064433075980   0.160945891394200]'*1e7;
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

% from JAT's ChargedParticleModel.main() output: 0.0016458877075534145 % (m)
Delay = jatChargedParticleModel(Ephem);

tol = 1e-10;
failed = 0;

Diff = Delay - 0.0016458877075534145 %#ok<NOPRT>

if any( abs(Diff) > tol ) || any(isnan(Diff))
    failed = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jalChargedParticleModel_regression_test - jatChargedParticlModel Regression Tests
%
% tests output against previous values

function failed = jatChargedParticleModel_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/09/2008	 	Original
%   Ravi Mathur            08/29/2012      Extracted from jatChargedParticleModel.m

disp(' ')
disp('Performing Test....')
disp(' ')

% Case 1 - test against Jat
Ephem.satPos = [2.056020937918350   1.411064433075980   0.160945891394200]'*1e7;
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

% from JAT's ChargedParticleModel.main() output: 0.0016458877075534145 (m)
Delay(1) = jatChargedParticleModel(Ephem); 

% Case 2 - test multiple satellites, multiple epochs and frequencies
Ephem.satPos = [[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7,...
               2*[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7];
Ephem.Epoch = [datenum('June 9, 2008'),datenum('June 10, 2008')];
Ephem.SignalFreq = [1.6e9,1.8e9];

Delay(2:3) = jatChargedParticleModel(Ephem);

% Case 2 - test multiple satellites, single epoch and frequency
Ephem.satPos = [[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7,...
               2*[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7];
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

Delay(4:5) = jatChargedParticleModel(Ephem) %#ok<NOPRT>

% uncomment to correct regression data:
%PreviousDelay = Delay;
%save jatCPM_RegressionData6_08 PreviousDelay

tol = 1e-10;
failed = 0;
if exist('jatCPM_RegressionData11_10.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load regression values
	load jatCPM_RegressionData11_10

	Diff = Delay - PreviousDelay %#ok<NOPRT>

	if any( abs(Diff) > tol ) || any(isnan(Diff))
		failed = 1;
	end
else
	failed = 1;
end


end

