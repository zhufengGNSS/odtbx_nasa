function [dRho, dE] = jatTropoModel(rho, elev, lambda, pressure, temp, relHumidity )

% JATTROPOMODEL  Returns tropospheric refraction errors for range and elevation measurements.
%
%   [dRho, dE] = jatTropoModel(rho, elev, lambda, pressure, temp, relHumidity )
% returns the refraction errors from JAT for range and elevation measurements.
%   
% The tropo model is the Hopfield-Goad model described in Section 6.3 of 
% Montenbruck's book. It says the model is applicable for radar and 
% optical measurements. This model may not be relevant for space-to-space since 
% the troposphere is the lower layer, up to ~42 km in altitude.
%
%   INPUTS 
%   VARIABLE        SIZE   	DESCRIPTION (Optional/Default)
%      rho          (1xN)   range in m
%      elev         (1xM)   elevation angle in radians (default = pi/2)
%      lambda       (1x1)   wavelength in meters (default = L1)
%      pressure     (1x1)   pressure in hPa (default = 938)
%      temp         (1x1)   temperature in deg K (default = 286)
%      relHumidity  (1x1)   relative humidity (0 <= relHumidity <= 1)
%                             (default = 0.73) 
%
%   OUTPUTS 
%      dRho         (NxM)   tropospheric refraction error for range (m)
%      dElev        (NxM)   elevation error (arcsec)
%
%   REFERENCE
%
%   Montenbruck, Oliver and Gill, Eberhard, Satellite Orbits: Methods, Models
%   and Applications, Springer-Verlag, 2001, pp. 221-224.
%
%  VALIDATION TEST
%
%   To perform a validation test, enter 'ValidationTest' as the only input 
%   argument.  If the data file is not in the path, this will perform as an 
%   example.
%
%  REGRESSION TEST
%
%   To perform a regression test, enter 'RegressionTest' as the only input 
%   argument.  If the data file is not in the path, this will perform as an 
%   example.
%
%   keyword: JAT Adapter, measurement model, atmosphere model
%   See also JATGPSIONOMODEL
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%                      (MM/DD/YYYY)
%   Derek Surka         08/20/2007      Original
%   Keith Speckman      06/04/2008      Added Validation and Regression Tests
%   Kevin Berry         12/18/2008      Moved nargin checks out of the
%                                       getdRhodE function so that the user
%                                       can still pass in fewer inputs.

if strcmpi(rho,'ValidationTest')

	dRho = jatTropoModel_validation_test();

elseif strcmpi(rho,'RegressionTest')

	dRho = jatTropoModel_regression_test();
else

    % Set defaults
    if( nargin < 2 )
        elev = pi/2;
    end
    if( nargin < 3 )
        lambda = JATConstant('L1wavelength');
    end
    if( nargin < 4 )
        pressure = 938; % hPa
    end
    if( nargin < 5 )
        temp = 286; % deg K
    end
    if( nargin < 6 )
        relHumidity = 0.73;
    end

    [dRho, dE] = getdRhodE(rho, elev, lambda, pressure, temp, relHumidity );
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function[dRho, dE] = getdRhodE(rho, elev, lambda, pressure, temp, relHumidity )

n = length(rho);
m = length(elev);
dRho = zeros(n,m);
dE   = zeros(n,m);

% Get errors from JAT
tropoModel = jat.ground_tracking.TropoModel(lambda);

for i=1:n 
    for j=1:m 
        y = tropoModel.corrections(pressure, temp, relHumidity, elev(j), rho(i));
        dRho(i,j) = y(1);
        dE(i,j)   = y(2);
    end
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

[dRho, dE] = getdRhodE(rho, els, lambda, p, T, fh );

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

[dRho, dE] = getdRhodE(rho, els, lambda, p, T, fh );

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



