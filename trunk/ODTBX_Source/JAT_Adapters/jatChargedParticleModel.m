function Delay = jatChargedParticleModel(Ephem)

% JATCHARGEDPARTICLEMODEL Returns the light time delay error due to charged particles
%
%    Delay = jatChargedParticleModel(Ephem)
%
%  This function calculates the light time delay due to charged particles.  A
%  series of satellite position vectors can be input in one function call with
%  corresponding epochs and signal frequencies.  If the epoch or signal
%  frequency is a scalar rather than a vector, the single value will be used
%  for all spacecraft positions.
%
% INPUTS 
%   VARIABLE            SIZE      DESCRIPTION
%   Ephem.satPos        (3xN)     ECI satellite coordinates (m)
%   Ephem.Epoch      (1xN or 1x1) UTC epoch in Matlab datenum time format
%   Ephem.SignalFreq (1xN or 1x1) Signal Frequency (Hz)
%
% OUTPUT 
%   Delay              (1xN)     Signal delay error due to charged particles (m)
%
% REFERENCE
%
%  Burkhart, Paul D., "Adaptive Orbit Determination for Interplanetary
%  Spacecraft," Ph.D. dissertation, The University of Texas at Austin, May 1995,
%  pp. 27-30.
%
% VALIDATION TEST
%
%  To perform a validation test, replace the satPos variable with 
%  'ValidationTest' as the input argument.  If the data file is not in the path
%  this will perform as an example.
%
% REGRESSION TEST
%
%  To perform a regression test, replace the Ephem structure with 
%  'RegressionTest' as the input argument.  If the data file is not in the path
%  this will perform as an example
%
%  keywords: JAT Adapter, delay
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

% REVISION HISTORY
%  Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%  Keith Speckman       06/09/2008   	Original
%  Kevin Berry          06/22/2009      Corrected the help to specify the
%                                       input time scale of UTC  


% Determine whether this is an actual call to the program or a test

if strcmpi(Ephem,'ValidationTest')

	Delay = jatChargedParticleModel_validation_test();

elseif strcmpi(Ephem,'RegressionTest')

	Delay = jatChargedParticleModel_regression_test();

else

	Delay = getDelay(Ephem);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function Delay = getDelay(Ephem)

% Determine size of output
N = size(Ephem.satPos,2);

% Initialize Delay for quicker execution time
Delay = zeros(N,1);

% Create charged particle java object
ModelObject = jat.ground_tracking.ChargedParticleModel(Ephem.SignalFreq(1));

% Create counter for multiple epochs
m = 1;

for n = 1:N

	% Update charged particle model if multiple frequencies are provided
	if length(Ephem.SignalFreq) > 1
		ModelObject = jat.ground_tracking.ChargedParticleModel(Ephem.SignalFreq(n));
	end

	% Update epoch counter if multiple epochs are provided
	if length(Ephem.Epoch) > 1
		m = n;
	end

	% Create java objects for time and satellite position
	TimeObject = jat.spacetime.Time(matlabTime2MJD(Ephem.Epoch(m)));
	satObject = jat.matvec.data.VectorN(Ephem.satPos(:,n));

	% Compute the delay distance
	Delay(n) = ModelObject.returnDelay(satObject,TimeObject);

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

disp(' ')
disp('Performing Test....')
disp(' ')

% Case 1 - test against Jat
Ephem.satPos = [2.056020937918350   1.411064433075980   0.160945891394200]'*1e7;
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

% from JAT's ChargedParticleModel.main() output: 0.0016458877075534145 % (m)
Delay = getDelay(Ephem);

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

disp(' ')
disp('Performing Test....')
disp(' ')

% Case 1 - test against Jat
Ephem.satPos = [2.056020937918350   1.411064433075980   0.160945891394200]'*1e7;
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

% from JAT's ChargedParticleModel.main() output: 0.0016458877075534145 (m)
Delay(1) = getDelay(Ephem); 

% Case 2 - test multiple satellites, multiple epochs and frequencies
Ephem.satPos = [[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7,...
               2*[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7];
Ephem.Epoch = [datenum('June 9, 2008'),datenum('June 10, 2008')];
Ephem.SignalFreq = [1.6e9,1.8e9];

Delay(2:3) = getDelay(Ephem);

% Case 2 - test multiple satellites, single epoch and frequency
Ephem.satPos = [[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7,...
               2*[2.056020937918350   1.411064433075980   0.160945891394200]'*1e7];
Ephem.Epoch = datenum('June 9, 2008');
Ephem.SignalFreq = 1.6e9;

Delay(4:5) = getDelay(Ephem) %#ok<NOPRT>

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

