function failed = srpAccel_test(type)
%
% srpAccel_test Regression test for srpAccel
% See also: srpAccel.m
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
%   Ravi Mathur         08/28/2012      Extracted from srpAccel.m

if strcmpi(type,'ValidationTest')

	failed = srpAccel_validation_test();

elseif strcmpi(type,'RegressionTest')

	failed = srpAccel_regression_test();

else
    disp('srpAccel_test: Unsupported test type ')
    failed = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% srpAccel_validation_test - Validation for srpAccel.m

function failed = srpAccel_validation_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/5/2008	 	Original
%   Ravi Mathur         08/28/2012      Extracted from srpAccel.m

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-7;
format compact

%%%%%%%%%
Case = 1; % Earth orbit - no shadow - plate normal in sun direction
disp('Case = 1;  Earth orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 2; % Earth orbit - no shadow - 2 plates normal in sun direction
disp('Case = 2;  Earth orbit - no shadow - 2 plates normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = [1.0,1.0];
SC(Case).c_diffuse = [0.5,0.5];
SC(Case).NormVec_B = [1 0 0;1 0 0]';
SC(Case).area = [5,5];
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 3; % Saturn orbit - no shadow - plate normal in sun direction
disp('Case = 3;  Saturn orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = -ephemDE405('earth',Epoch(Case),'UTC') + ephemDE405('saturn',Epoch(Case),'UTC')
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

RsunUnit = unit(-SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC'));
CrossProd = cross(SC(Case).NormVec_B,RsunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,RsunUnit)) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 4; % Earth orbit - no shadow - plate normal edge-on to sun direction
disp('Case = 4;  Earth orbit - no shadow - plate normal edge-on to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/2,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 5; % Earth orbit - no shadow - plate normal 45 deg to sun direction
disp('Case = 5;  Earth orbit - no shadow - plate normal 45 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/4,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 6; % Earth orbit - shadow - plate normal to sun direction
disp('Case = 6;  Earth orbit - shadow - plate normal to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
CrossProd = cross([0 1 0]',rSunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
DCecef2sun = inv(dcm('axs',RotAngle+pi,RotVector));

SC(Case).satPos = DCecef2sun*[0 2 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 7; % Earth orbit - no shadow - plate normal 135 deg to sun direction
disp('Case = 7;  Earth orbit - no shadow - plate normal 135 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+3*pi/4,RotVector+3*pi/4);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%% Validation Tests %%%%

if exist('srpAccel_ValidationData6_08.mat') == 2
	disp(' ')
	disp('Performing Validation...')
	disp(' ')


	%%%%%%%%%
	Case = 1;

	JatSRPObject = jat.forces.SolarRadiationPressure(SC(Case).mass,SC(Case).area,...
                                                 1+SC(Case).c_specular+2/3*SC(Case).c_diffuse);
	Rsun = ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC');
	JatsatPos = jat.matvec.data.VectorN(SC(Case).satPos*1e3);
	JatRsun = jat.matvec.data.VectorN(Rsun*1e3);
	PredictedAccel(:,Case) = JatSRPObject.accelSRP(JatsatPos,JatRsun).x;

	%%%%%%%%%
	Case = 2;
	PredictedAccel(:,Case) = 2*PredictedAccel(:,1);

	%%%%%%%%%
	Case = 3;

	JatSRPObject = jat.forces.SolarRadiationPressure(SC(Case).mass,SC(Case).area,...
                                               	1+SC(Case).c_specular+2/3*SC(Case).c_diffuse);

	Rsun = -SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC');
	JatsatPos = jat.matvec.data.VectorN(SC(Case).satPos*1e3);
	JatRsun = jat.matvec.data.VectorN(Rsun*1e3);
	PredictedAccel(:,Case) = JatSRPObject.accelSRP(JatsatPos,JatRsun).x;

	%%%%%%%%%
	Case = 4;
	PredictedAccel(:,Case) = [0 0 0]';

	%%%%%%%%%
	Case = 5;
	PredictedAccel(:,Case) = inv(DCi2b{Case}) * (-SC(Case).NormVec_B) * ...
	                       norm(PredictedAccel(:,1))*cos(pi/4);

	%%%%%%%%%
	Case = 6;
	PredictedAccel(:,Case) = [0 0 0]';

	%%%%%%%%%
	Case = 7;
	PredictedAccel(:,Case) = [0 0 0]';


	%save srpAccel_ValidationData6_08 PredictedAccel
	clear PredictedAccel % get data from data file rather than this instance

	% Load validation values PredictedAccel
	load srpAccel_ValidationData6_08

	disp(' ')
	disp('Test Results (m/s^2)')
	disp(' ')
	disp('        Test             Predicted')
	disp('----------------------------------------')

	for k = 1:7
		fprintf('\nCase %d\n\n',k)
		disp([Accel(:,k),PredictedAccel(:,k)])
	end

	disp(' ')


	Diff = [PredictedAccel - Accel]

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end


	%%%%%%%%%
	Case = 8; % Test Earth shadowing
	disp(' ')
	disp('Case 8; Testing Earth shadowing ...')
	disp(' ')

	Re = JATConstant('rEarth')/1e3;
	Xgrid = [-3:0.1:3] * Re;
	Ygrid = [5:-10:-300] * Re;
	ShadowBodies = {'EARTH'};
	Epoch(Case) = datenum('May 16, 2008'); %UTC
	aunit = JATConstant('au')/1e3; % Astronomical Unit (km)

	rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
	CrossProd = cross([0 1 0]',rSunUnit);
	RotVector = unit(CrossProd);
	RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
	DCecef2sun = inv(dcm('axs',RotAngle,RotVector));

	c_shad_totalx = zeros(length(Xgrid),length(Ygrid));
	c_shad_totalz = c_shad_totalx;

	for y = 1:length(Ygrid)
		for x = 1:length(Xgrid)
			satPosECEF = [Xgrid(x);Ygrid(y);0];
			satPos = DCecef2sun * satPosECEF;
			c_shad_totalx(x,y) = getBodyShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
		for z = 1:length(Xgrid)
			satPosECEF = [0;Ygrid(y);Xgrid(z)];
			satPos = DCecef2sun * satPosECEF;
			c_shad_totalz(z,y) = getBodyShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
	end

	figure(1)
		subplot(211)
			pcolor(-Ygrid/2/Re,Xgrid/2/Re,c_shad_totalx)
			title(['Case 8 - Earth Shadow Plot Varied in 2 directions.  ',...
	       		'Red = Sun, Green = Penumbra, Blue = Umbra'])
			ylabel('X-axis (Earth Diameters)')
		subplot(212)
			pcolor(-Ygrid/2/Re,Xgrid/2/Re,c_shad_totalz)
			ylabel('Z-axis (Earth Diameters)')
			xlabel('Y-axis (Earth Diameters)')

	%%%%%%%%%
	Case = 9; % Test Neptune shadowing
	% this should produce a warning since no neptune shadowing is currently available
	disp(' ')
	disp('Case 9; Testing Neptune shadowing ...')
	disp(' ')
    
    fprintf(1, '\n');
    fprintf(1, 'Note that this should produce a warning since Neptune\n');
    fprintf(1, 'is not currently available.\n');
    fprintf(1, '\n');

	Epoch(Case) = datenum('May 16, 2008'); %UTC
	Rn = JATConstant('meanRadius','neptune')/1e3;
	NeptunePos = ephemDE405('neptune',Epoch(Case),'UTC') - ephemDE405('EARTH',Epoch(Case),'UTC');


%	Xgrid = [-2:0.1:2] * Rn;
%	Ygrid = [5:-500:-8000] * Rn;
	Xgrid = [0.1] * Rn;
	Ygrid = [-500] * Rn;
	ShadowBodies = {'neptune'};
	aunit = JATConstant('au')/1e3; % Astronomical Unit (km)

	rSunUnit = unit(ephemDE405('SUN',Epoch(Case),'UTC') - ephemDE405('neptune',Epoch(Case),'UTC'));
	CrossProd = cross([0 1 0]',rSunUnit);
	RotVector = unit(CrossProd);
	RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
	DCecef2sun = inv(dcm('axs',RotAngle,RotVector));

	c_shad_totalx = zeros(length(Xgrid),length(Ygrid));
	c_shad_totalz = c_shad_totalx;

	for y = 1:length(Ygrid)
		for x = 1:length(Xgrid)
			satPosECEF = [Xgrid(x);Ygrid(y);0];
			satPos = DCecef2sun * satPosECEF + NeptunePos;
			c_shad_totalx(x,y) = getBodyShadow(Epoch(Case),satPos,ShadowBodies,aunit);
		end
%		for z = 1:length(Xgrid)
%			satPosECEF = [0;Ygrid(y);Xgrid(z)];
%			satPos = DCecef2sun * satPosECEF + NeptunePos;
%			c_shad_totalz(z,y) = getBodyShadow(Epoch(Case),satPos,ShadowBodies,aunit);
%		end
	end

%	figure(2)
%		subplot(211)
%			pcolor(-Ygrid/2/Rn,Xgrid/2/Rn,c_shad_totalx)
%			title(['Case 9 - Neptune Shadow Plot Varied in 2 directions.  ',...
%			       'Red = Sun, Green = Penumbra, Blue = Umbra'])
%			ylabel('X-axis (Neptune Diameters)')
%		subplot(212)
%			pcolor(-Ygrid/2/Rn,Xgrid/2/Rn,c_shad_totalz)
%			ylabel('Z-axis (Neptune Diameters)')
%			xlabel('Y-axis (Neptune Diameters)')


else
	failed = 1;
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% srpAccel_regression_test - Regression for srpAccel.m

function failed = srpAccel_regression_test()

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman         06/8/2008	 	Original
%   Ravi Mathur         08/28/2012      Extracted from srpAccel.m

disp(' ')
disp('Performing Test....')
disp(' ')

failed = 0;
tol = 1e-12;
format compact

%%%%%%%%%
Case = 1; % Earth orbit - no shadow - plate normal in sun direction
disp('Case = 1;  Earth orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

%%%%%%%%%
Case = 2; % Earth orbit - no shadow - 2 plates normal in sun direction
disp('Case = 2;  Earth orbit - no shadow - 2 plates normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = [1.0,1.0];
SC(Case).c_diffuse = [0.5,0.5];
SC(Case).NormVec_B = [1 0 0;1 0 0]';
SC(Case).area = [5,5];
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B(:,1),unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 3; % Saturn orbit - no shadow - plate normal in sun direction
disp('Case = 3;  Saturn orbit - no shadow - plate normal in sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = -ephemDE405('earth',Epoch(Case),'UTC') + ephemDE405('saturn',Epoch(Case),'UTC')
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

RsunUnit = unit(-SC(Case).satPos + ephemDE405('geocentric_sun',Epoch(Case),'UTC'));
CrossProd = cross(SC(Case).NormVec_B,RsunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,RsunUnit)) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 4; % Earth orbit - no shadow - plate normal edge-on to sun direction
disp('Case = 4;  Earth orbit - no shadow - plate normal edge-on to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/2,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));


%%%%%%%%%
Case = 5; % Earth orbit - no shadow - plate normal 45 deg to sun direction
disp('Case = 5;  Earth orbit - no shadow - plate normal 45 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+pi/4,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 6; % Earth orbit - shadow - plate normal to sun direction
disp('Case = 6;  Earth orbit - shadow - plate normal to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

rSunUnit = unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC'));
CrossProd = cross([0 1 0]',rSunUnit);
RotVector = unit(CrossProd);
RotAngle = acos(dot([0 1 0]',rSunUnit)) * sign(asin(norm(CrossProd)));
DCecef2sun = inv(dcm('axs',RotAngle+pi,RotVector));

SC(Case).satPos = DCecef2sun*[0 2 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle,RotVector);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

%%%%%%%%%
Case = 7; % Earth orbit - no shadow - plate normal 135 deg to sun direction
disp('Case = 7;  Earth orbit - no shadow - plate normal 135 deg to sun direction')

fprintf(1, '\n');
fprintf(1, 'Note that specular and diffuse coefficients are deliberately\n');
fprintf(1, 'set to sum to > 1.0 in order to ensure that this condition is\n');
fprintf(1, 'caught by the code.\n');
fprintf(1, '\n');

Epoch(Case) = datenum('May 16, 2008'); %UTC

SC(Case).satPos = [0 10 0]'*JATConstant('rEarth')/1e3;
SC(Case).c_specular = 1.0;
SC(Case).c_diffuse = 0.5;
SC(Case).NormVec_B = [1 0 0]';
SC(Case).area = 5;
SC(Case).mass = 2000;
SC(Case).ShadowBodies = {'EARTH','MOON'};

CrossProd = cross(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')));
RotVector = unit(CrossProd);
RotAngle = acos(dot(SC(Case).NormVec_B,unit(ephemDE405('GEOCENTRIC_SUN',Epoch(Case),'UTC')))) * ...
           sign(asin(norm(CrossProd)));
DCi2b{Case} = dcm('axs',RotAngle+3*pi/4,RotVector+3*pi/4);

Accel(:,Case) = srpAccel(Epoch(Case),DCi2b{Case},SC(Case));

displayInOut(Case,Epoch(Case),DCi2b{Case},SC(Case),Accel(:,Case));

disp(' ')

%PreviousAccel = Accel;
%save srpAccel_RegressionData6_08 PreviousAccel

%%%% Validation Tests %%%%

if exist('srpAccel_RegressionData6_08.mat') == 2
	disp(' ')
	disp('Performing Regression Test...')
	disp(' ')


	% Load validation values PreviousAccel
	load srpAccel_RegressionData6_08

	disp(' ')
	disp('Test Results (m/s^2)')
	disp(' ')
	disp('        Test             Previous')
	disp('----------------------------------------')

	for k = 1:7
		fprintf('\nCase %d\n\n',k)
		disp([Accel(:,k),PreviousAccel(:,k)])
	end

	disp(' ')


	Diff = [PreviousAccel - Accel]

	if any(any( abs(Diff) > tol )) | any(any(isnan(Diff)))
		failed = 1;
	end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function displayInOut(Case,Epoch,DCi2b,SC,Accel)

disp(' ')
disp(['Case = ',num2str(Case)])
disp(['Epoch = ',num2str(Epoch),' UTC'])
DCi2b
fprintf('SC.satPos = \n')
disp(SC.satPos)
disp(['SC.c_specular = ',num2str(SC.c_specular)])
disp(['SC.c_diffuse = ',num2str(SC.c_diffuse)])
fprintf('SC.NormVec_B = \n')
disp(SC.NormVec_B)
disp(['SC.area = ',num2str(SC.area)])
disp(['SC.mass = ',num2str(SC.mass)])
fprintf('SC.ShadowBodies = \n')
disp(SC.ShadowBodies)
fprintf('Accel = \n')
disp(Accel)
disp(' ')

end
