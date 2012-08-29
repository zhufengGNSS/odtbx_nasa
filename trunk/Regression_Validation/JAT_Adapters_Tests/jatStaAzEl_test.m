function failed = jatStaAzEl_test(type)
% jatStaAzEl_test Regression test for jatStaAzEl
% See also: jatStaAzEl.m
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
%   Ravi Mathur         08/29/2012      Extracted from jatStaAzEl.m

if strcmpi(type,'ValidationTest')

	failed = jatStaAzEl_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = jatStaAzEl_regression_test();

else
    disp('jatStaAzEl_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jatStaAzEl_validation_test - jatStaAzEl Validation Tests

function failed = jatStaAzEl_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original
%   Keith Speckman         06/27/2008       Changed examples from SEZ to NED
%   Allen Brown            04/09/2009       Updated data for coincident
%                                           satellite and station
%   Ravi Mathur            08/29/2012      Extracted from jatStaAzEl.m

disp(' ')
disp('Performing Test....')
disp(' ')

fprintf('%17s%17s%17s%17s\n','Station Pos','Satellite Pos','Expected Az','Az Results')
fprintf('%17s%17s%17s%17s\n','/ Rearth   ','/ Rearth   ',' El (deg)','  El (deg)')
fprintf('%s\n\n',char(ones(1,17*4)*'-'));

tol = 1e-7;

Re = JATConstant('rEarth');

Ephem.Epoch = datenum('May 19, 2008');

Case = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 90;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);


Case = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 270;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = 90;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);
disp('Azimuth undefined')
disp(' ')

Case = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[0 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

% A coincident satellite and station should have az & el of zero.
AzEx = 0;
ElEx = 0;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);
disp('Azimuth undefined')
disp('Elevation undefined')
disp(' ')

Case = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[0 2 1]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = atan(2/1)*180/pi;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

AzEx = 270;
ElEx = -45;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.lat = 0;
Ephem.lon = 0;
Ephem.height = -0.0007;
Ephem.satPos = Ephem.staPos + Re*[1 0 1]';
Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

AzEx = 0;
ElEx = 45;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);

Case = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = jatDCM('ecef2eci',Ephem.Epoch) * [Ephem.staPos + Re*[1 0 0]'];
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

AzEx = 0;
ElEx = 90;

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1)/Re,Ephem.satPos(1)/Re,AzEx,...
   azimuth(Case,1)*180/pi);
fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2)/Re,Ephem.satPos(2)/Re,ElEx,...
   elevation(Case,1)*180/pi);
fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3)/Re,Ephem.satPos(3)/Re);


Case = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos(:,1) = Re*[1 0 0]';
Ephem.satPos(:,1) = Ephem.staPos(:,1) + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

Ephem.staPos(:,2) = Re*[0 1 0]';
Ephem.satPos(:,2) = Ephem.staPos(:,2) + Re*[1 1 0]';

Ephem.satPos(:,3) = Ephem.staPos(:,2) + Re*[0 1 2]';

AzEx(1,1) = 90;
ElEx(1,1) = 45;
AzEx(1,2) = 270;
ElEx(1,2) = 0;
AzEx(2,1) = 90;
ElEx(2,1) = 0;
AzEx(2,2) = 270;
ElEx(2,2) = 45;
AzEx(3,1) = 45;
ElEx(3,1) = -atan(1/sqrt(8))*180/pi;
AzEx(3,2) = 0;
ElEx(3,2) = atan(1/2)*180/pi;

[azimuthtemp,elevationtemp] = jatStaAzEl(Ephem);
fprintf('********* Case %d, Staion = %s, Sat = %s ***********\n\n',Case,Ephem.StationInfo,...
   Ephem.SatCoords)
disp('Matrix output')
disp(' ')

casenum = Case-1;
for n = 1:3
	for m = 1:2
		casenum = casenum + 1;
		azimuth(casenum,1) = azimuthtemp(n,m);
		elevation(casenum,1) = elevationtemp(n,m);
		fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(1,m)/Re,Ephem.satPos(1,n)/Re,...
		        AzEx(n,m),azimuth(casenum)*180/pi);
		fprintf('%17.5f%17.5f%17.5f%17.5f\n',Ephem.staPos(2,m)/Re,Ephem.satPos(2,n)/Re,...
		        ElEx(n,m),elevation(casenum)*180/pi);
		fprintf('%17.5f%17.5f\n\n',Ephem.staPos(3,m)/Re,Ephem.satPos(3,n)/Re);
	end
end

failed = 0;

if exist('jatStaAzEl_ValidationData4_09.m') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')

	clear AzEx ElEx

	% Load validation values for AzEx and ElEx
	jatStaAzEl_ValidationData4_09

	Diff = [AzEx - azimuth;ElEx - elevation]

	if any( abs(Diff) > tol ) | any(isnan(Diff([1:end])))
		failed = 1;
	end
else
	failed = 1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jatStaAzEl_regression_test - jatStaAzEl Regression Tests

function failed = jatStaAzEl_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         05/19/2008	 	Original
%   Keith Speckman         06/27/2008       Changed examples from SEZ to NED
%   Allen Brown            04/09/2009       Updated data for coincident
%                                           satellite and station
%   Ravi Mathur            08/29/2012      Extracted from jatStaAzEl.m

disp(' ')
disp('Performing Demonstration....')

failed = 0;
tol = 1e-7;

Re = JATConstant('rEarth');

Ephem.Epoch = datenum('May 19, 2008');

Case = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[1 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = Ephem.staPos + Re*[0 0 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[0 2 1]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[0 1 0]';
Ephem.satPos = Ephem.staPos + Re*[1 -1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.lat = 0;
Ephem.lon = 0;
Ephem.height = -0.0007;
Ephem.satPos = Ephem.staPos + Re*[1 0 1]';
Ephem.StationInfo = 'LatLonHeight';
Ephem.SatCoords = 'ECEF';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos = Re*[1 0 0]';
Ephem.satPos = jatDCM('ecef2eci',Ephem.Epoch) * [Ephem.staPos + Re*[1 0 0]'];
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECI';

[azimuth(Case,1),elevation(Case,1)] = jatStaAzEl(Ephem);

Case = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ephem.staPos(:,1) = Re*[1 0 0]';
Ephem.satPos(:,1) = Ephem.staPos(:,1) + Re*[1 1 0]';
Ephem.StationInfo = 'ECEF';
Ephem.SatCoords = 'ECEF';

Ephem.staPos(:,2) = Re*[0 1 0]';
Ephem.satPos(:,2) = Ephem.staPos(:,2) + Re*[1 1 0]';

Ephem.satPos(:,3) = Ephem.staPos(:,2) + Re*[0 1 2]';

[azimuthtemp,elevationtemp] = jatStaAzEl(Ephem);

casenum = Case-1;
for n = 1:3
	for m = 1:2
		casenum = casenum + 1;
		azimuth(casenum,1) = azimuthtemp(n,m);
		elevation(casenum,1) = elevationtemp(n,m);
	end
end
azimuth
elevation

format long g

current_values = [azimuth,elevation];

%save jatStaAzEl_RegressionData6_08 azimuth elevation

failed = 0;

if exist('jatStaAzEl_RegressionData4_09.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')

	clear azimuth elevation

	% Load regression values for AzEx and ElEx
	load jatStaAzEl_RegressionData4_09

	Diff = [azimuth,elevation] - current_values

	if any(any(abs(Diff)>tol)) | any(isnan([Diff(:,1);Diff([1:end],2)]))
		failed = 1;
	end
else
	failed = 1;
end


end