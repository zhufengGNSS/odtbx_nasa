function [lat,lon,alt] = ecef2LLA(x,ErrTol)

% ECEF2LLA  Returns the latitude, longitude and altitude of an ECEF position.
%
%   [lat,lon,alt] = ecef2LLA(x) returns the latitude, longitude and
%   altitude of input ECEF positions. Latitude and longitude both range
%   from -180 to 180 (deg).
%
%   INPUTS
%   VARIABLE    SIZE        DESCRIPTION (Optional/Default)
%     x         (3xN)       ECEF Positions in km
%     ErrTol    (1x1)       (optional) error (km) above which a warning is
%                           produced; default = 1e-9 km
%
%   OUTPUTS
%     lat       (1xN)       Latitudes in degs
%     lon       (1xN)       Longitudes in degs
%     alt       (1xN)       Altitudes in km
%
% VALIDATION TEST
%
%   To perform a validation test, replace the x input with 'ValidationTest'
%   as the input argument.  If the data file is not in the path this will
%   perform as an example.
%
%   The proper invocation for this function for validation test is:
%
%   [fail] = ecef2LLA('ValidationTest')
%
% REGRESSION TEST
%
%   To perform a regression test, replace the x input with 'RegressionTest'
%   as the input argument.  If the data file is not in the path this will
%   perform as an example
%
%   The proper invocation for this function regression test is:
%
%   [fail] = ecef2LLA('RegressionTest')
%
%   keyword: Coordinate Transformations, JAT Adapter
%   See also JATDCM
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
%                           (MM/DD/YYYY)
%   Derek Surka              09/18/2007         Original
%   Keith Speckman           05/30/2008       Corrected error in JAT call
%                                             Modified JAT algorithm
%                                             Added regression and validation
%   Allen Brown              3/30/2009        Fixed regression check logic,
%                                             fixed regression data error
%                                             (now using
%                                             ecef2LLA_RegressionData6_08b)

if strcmp(x,'ValidationTest')

	lat = ecef2LLA_validation_test();

elseif strcmp(x,'RegressionTest')

	lat = ecef2LLA_regression_test();
else


	if nargin == 1
		ErrTol = 1e-9; %km
	end

	[lat,lon,alt] = getLLA(x,ErrTol);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function

function [lat,lon,alt] = getLLA(x,ErrTol);

x = x*1000; % m

n   = size(x,2);
lat = zeros(1,n);
lon = zeros(1,n);
alt = zeros(1,n);

for k=1:n
    LLA = jat.groundstations.GroundStation.getLLA(x(:,k));
    lat(1,k) = LLA(1)*180/pi;
    lon(1,k) = LLA(2)*180/pi;
    alt(1,k) = LLA(3)/1000; %km
    ecef_error = LLA(4:6); %m

    % If the error excedes the tolerance, produce a warning
    % This may be due to the JAT iteration scheme not converging resulting in a small error
    if norm(ecef_error) > ErrTol*1e3
	    warning(['Error exceeds tolerance.  Error = ',num2str(norm(ecef_error/1e3)),' km'])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ecef2LLA_validation_test - ecef2LLA validation Tests

function failed = ecef2LLA_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/03/2008	 	Original

disp(' ')
disp('Performing Validation Test....')

failed = 0;
tol = 1e-8;

Re = JATConstant('rEarth')/1e3;

LatLonAltValues = [	0	0	0
		0	90	0
		0	180	0
		0	-90	0
		0	45	0
		0	135	0
		0	-135	0
		0	-45	0
		90	0	0
		-90	0	0
		45	90	0
		-45	-90	0
		0	0	0.1*Re
		0	-90	0.1*Re
		90	0	0.1*Re	];

for k = 1:size(LatLonAltValues,1)
	x(k,:) = [LLA2ecef(LatLonAltValues(k,1),LatLonAltValues(k,2),LatLonAltValues(k,3))]';
end

[lat,lon,height] = ecef2LLA(x');

% intput a negative tolerance in order to test tolerance warning

fprintf(1, '\n');
fprintf(1, 'Note that a negative tolerance has been deliberately input\n');
fprintf(1, 'in order to ensure that it triggers a warning message.\n');
fprintf(1, '\n');

x(k+1,:) = LLA2ecef(20,35,25)';
[lat(k+1),lon(k+1),height(k+1)] = ecef2LLA(x(k+1,:)',-1);

format short g
lstring = '--------------';
disp(['     x(1)     ', ...
      '     x(2)     ', ...
      '     x(3)     ', ...
      '   lat (deg)  ', ...
      '   lon (deg)  ', ...
      '  height (km) '])
disp([lstring,lstring,lstring,lstring,lstring,lstring])
disp([x,lat',lon',height'])

if exist('ecef2LLA_ValidationData6_08.m') == 2
	disp(' ')

	clear EstimatedValues

	% Load validation values for AzEx and ElEx
	ecef2LLA_ValidationData6_08

	Diff = [LatLonAltValues - [lat',lon',height']]

	if any( any( abs(Diff) > tol )) | any(isnan(Diff))
		failed = 1;
	end
else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LLA2ecef_regression_test - LLA2ecef Regression Tests

function failed = ecef2LLA_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/03/2008	 	Original

disp(' ')
disp('Performing Regression Test...')
disp(' ')

tol = 1e-7;

Re = JATConstant('rEarth')/1e3;

x = Re * [	1		0		0
		0		1		0
		-1		0		0
		0		-1		0
		1/sqrt(2)	1/sqrt(2)	0
		-1/sqrt(2)	1/sqrt(2)	0
		-1/sqrt(2)	-1/sqrt(2)	0
		1/sqrt(2)	-1/sqrt(2)	0
		0		0		1
		0		0		-1
		0		1/sqrt(2)	1/sqrt(2)
		0		-1/sqrt(2)	-1/sqrt(2)
		1.1		0		0
		0		-1.1		0
		0		0		1.1	]';

[lat,lon,height] = ecef2LLA(x);


format short g
lstring = '--------------';
disp(['     x(1)     ', ...
      '     x(2)     ', ...
      '     x(3)     ', ...
      '   lat (deg)  ', ...
      '   lon (deg)  ', ...
      '  height (km) '])
disp([lstring,lstring,lstring,lstring,lstring,lstring])
disp([x',lat',lon',height'])

%PreviousValues = [lat',lon',height'];
%save ecef2LLA_RegressionData6_08 PreviousValues

failed = 0;
if exist('ecef2LLA_RegressionData6_08b.mat') == 2
	disp(' ')

	% Load regression values for AzEx and ElEx
	load ecef2LLA_RegressionData6_08b

    disp('Differences from regression data:');
	Diff = PreviousValues - [lat', lon', height'] 

	if any(any(abs(Diff) > tol )) || any(any(isnan(Diff)))
		failed = 1;
	end
else
	failed = 1;
end

end
