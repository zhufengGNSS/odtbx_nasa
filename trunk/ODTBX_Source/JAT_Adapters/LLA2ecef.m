function x = LLA2ecef(lat, lon, alt)

% LLA2ECEF  Returns the ECEF position from latitude, longitude and altitude
%
%  x = LLA2ecef(lat, lon, alt) returns the ECEF positions of the input latitude,
%  longitude and altitude positions.
%
% INPUTS
%   VARIABLE    SIZE        DESCRIPTION (Optional/Default)
%     lat       (1xN)       Latitudes in degs
%     lon       (1xN)       Longitudes in degs
%     alt       (1xN)       Altitudes in km
%
% OUTPUTS
%     x         (3xN)       ECEF Positions in km
%
% VALIDATION TEST
%
%   To perform a validation test, replace the lat, lon, and alt inputs with 
%   'ValidationTest' as the input argument.  If the data file is not in the path
%   this will perform as an example.
%
% REGRESSION TEST
%
%   To perform a regression test, replace the lat, lon, and alt inputs with 
%   'RegressionTest' as the input argument.  If the data file is not in the path
%   this will perform as an example
%
%   keywords: Coordinate Transformations, JAT Adapter
%   See also ECEF2LLA
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
%   Author                    Date        Comment
%                         (MM/DD/YYYY)
%   Keith Speckman         05/19/2008     Original

% Determine whether this is an actual call to the program or a test

if strcmp(lat,'ValidationTest')

	x = LLA2ecef_validation_test();

elseif strcmp(lat,'RegressionTest')

	x = LLA2ecef_regression_test();
else

	x = getECEF(lat,lon,alt);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function

function x = getECEF(lat,lon,alt)

lat = lat/180*pi;   % rads
lon = lon/180*pi;   % rads
alt = alt*1000;     % m

n = size(lat,2);
x = zeros(3,n);

for k=1:n
    station = jat.groundstations.GroundStation('tmp',lat(k),lon(k),alt(k));
    x(:,k) = station.getECEFPosition();
end


x = x/1000;         % back to km

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LLA2ecef_validation_test - LLA2ecef validation Tests

function failed = LLA2ecef_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original

disp(' ')
disp('Performing Test....')

tol = 30;  % (km) validation test data is based on a spherical earth
           % and is therefore only an approximation for validation
           % purposes

Re = JATConstant('rEarth')/1e3;

LatLonAlt = [	0	0	0
		0	90	0
		0	180	0
		0	270	0
		0	45	0
		0	135	0
		0	225	0
		0	315	0
		90	0	0
		-90	0	0
		45	90	0
		-45	270	0
		0	0	0.1*Re
		0	270	0.1*Re
		90	0	0.1*Re	];


EstimatedValues = Re * [	1		0		0
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
				0		0		1.1	];

x = zeros(size(EstimatedValues));

disp(' ')
disp('Expected values are based on a spherical earth approximation')
disp(' ')
disp('1.  Expected Values')
disp('2.  Test Values')
disp(' ')
fprintf(1,'%17s%17s%17s\n','Latitude','Longitude','Altitude')
fprintf(1,'----------------------------------------------------------\n')

for k = 1:size(LatLonAlt,1)
	x(k,:) = [getECEF(LatLonAlt(k,1),LatLonAlt(k,2),LatLonAlt(k,3))]';
	fprintf(1,'1. %17.8f%17.8f%17.8f\n',EstimatedValues(k,1),EstimatedValues(k,2),...
	        EstimatedValues(k,3))
	fprintf(1,'2. %17.8f%17.8f%17.8f\n\n',x(k,1),x(k,2),x(k,3))
end

disp(x)

failed = 0;

if exist('LLA2ecef_ValidationData5_08.m') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	clear EstimatedValues

	% Load validation values for AzEx and ElEx
	LLA2ecef_ValidationData5_08

	Diff = [EstimatedValues - x]

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end
else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LLA2ecef_regression_test - LLA2ecef Regression Tests

function failed = LLA2ecef_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original

disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;

Re = JATConstant('rEarth')/1e3;

LatLonAlt = [	0	0	0
		0	90	0
		0	180	0
		0	270	0
		0	45	0
		0	135	0
		0	225	0
		315	0	0
		90	0	0
		-90	0	0
		45	90	0
		-45	270	0
		0	0	0.1*Re
		0	270	0.1*Re
		90	0	0.1*Re	];

x = [getECEF(LatLonAlt(:,1)',LatLonAlt(:,2)',LatLonAlt(:,3)')]'

failed = 0;

if exist('LLA2ecef_RegressionData5_08.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	% Load regression values for AzEx and ElEx
	load LLA2ecef_RegressionData5_08

	Diff = PreviousValues - x

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end
else
	failed = 1;
end


end
