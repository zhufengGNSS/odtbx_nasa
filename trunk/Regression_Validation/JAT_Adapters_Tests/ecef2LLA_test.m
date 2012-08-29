function failed = ecef2LLA_test(type)
% ecef2LLA_test Regression test for ecef2LLA
% See also: ecef2LLA.m
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
%   Ravi Mathur         08/29/2012      Extracted from ecef2LLA.m

% THIS IS THE OLD REGRESSION TEST. It fails and has been replaced by
% the previously built-in tests from ecef2LLA.
% THIS SECTION SHOULD BE DELETED AFTER CCB DISCUSSION.
% failed  = 0;
% tol     = 1e-10;
% 
% x = [   0     0   6500    0
%         0    6500   0     0
%         6500  0     0   -6500 ];
%     
% latTarget = [   0                       90       0         0 ];
% lonTarget = [   90                       0       0       -90 ];
% altTarget = [  1.432476857548207e+002 121.863 121.863 1.432476857548207e+002 ]; 
%     
% [lat,lon,alt] = ecef2LLA(x);
% 
% if( (any(abs(latTarget-lat))>tol) || (any(abs(lonTarget-lon))>tol) || (any(abs(altTarget-alt))>tol) )
%     failed = 1;
% end
% END OLD REGRESSION TEST

if strcmpi(type,'ValidationTest')

	failed = ecef2LLA_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = ecef2LLA_regression_test();

else
    disp('ecef2LLA_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ecef2LLA_validation_test - ecef2LLA validation Tests

function failed = ecef2LLA_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/03/2008	 	Original
%   Ravi Mathur            08/29/2012      Extracted from ecef2LLA.m

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
%   Ravi Mathur            08/29/2012      Extracted from ecef2LLA.m

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
