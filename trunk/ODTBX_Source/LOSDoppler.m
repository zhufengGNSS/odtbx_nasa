function doppler = LOSDoppler(t, r1, v1, r2, v2, options)

% LOSDOPPLER  Line of sight instantaneous doppler between two objects
%
%   doppler = LOSDoppler(r1, v1, r2, v2, f) returns the doppler
% measurement between two vehicles and is simply a conversion of the range 
% rate calculation.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The parameters that are
% valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   useGPSIono          {true, false(default)} only for GPS sats as x2
%   useIono             {true, false(default)} only for groundstats as x2
%   useTropo            {true, false(default)} only for groundstats as x2
%   useChargedParticle  {true, false(default)} only for groundstats as x2
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only for Doppler and
%                                              measurement errors
%   epoch                datenum               UTC time associated with 
%                                              start of simulation.
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t         (1xN)	Times corresponding to r1 (secs)
%      r1        (3xN)	User spacecraft position (km)
%      v1        (3xN)  User spacecraft velocity (km/s)
%      r2        (3xN)  Tracking spacecraft/ground station position (km)
%      v2        (3xN)  Tracking spacecraft/ground station velocity (km/s)
%      options   (1x1)  Data structure
%
% OUTPUTS
%      doppler   (1xN)  doppler (Hz)
%
% VALIDATION TEST
%
%  To perform a validation test, pass in 'ValidationTest' as the
%  only input argument and specify only one output argument.
%
% REGRESSION TEST
%
%  To perform a regression test, pass in 'RegressionTest' as the
%  only input argument and specify only one output argument.  
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, RRDOT, RRDOTLT
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
%   Derek Surka         08/27/2007      Original
%   Kevin Berry         05/08/2008      Replaced with a simpler rdot*f/c
%                                       model
%   Kevin Berry         05/22/2008      Added options structure
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         09/08/2008      Added Charged Particle Model
%   Kevin Berry         06/25/2009      Added time scale comments

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    doppler = losdoppler_validation_test();  
else
	doppler = GetDoppler(t, r1, v1, r2, v2, options);
end
end


%% Main function
function doppler = GetDoppler(t, r1, v1, r2, v2, options)

c       = JATConstant('c')/1000;
f       = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
rdot    = LOSRangeRate(t, r1, v1, r2, v2, options);
doppler = -rdot*f/c;

end
%% Validation Test

function failed = losdoppler_validation_test()

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%12s%17s%19s%17s%18s%19s\n','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (Hz)','Calculated (Hz)')
fprintf('%s\n\n',char(ones(1,102)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006')); %UTC
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);

t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [e1 e1 e1 e1 e1 e1 e1 e1 e1];
r2 = [e2 e2 e2 e2 e2 e2 e2 e2 e2];    
v1 = [e1 e1 e1 e2 e2 e2 e3 e3 e3];
v2 = [e1 e2 e3 e1 e2 e3 e1 e2 e3];

ExDoppler = [0 -7431.72490153421 -3715.87121520219 7431.75995935723 ...
    0 3715.87121520219 3715.87997967862 -3715.8624507671 0];
Doppler = GetDoppler(t, r1, v1, r2, v2, options);

fprintf('%4.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %16.6f %16.6f\n',...
    [r1; r2; v1; v2; ExDoppler; Doppler]);

failed = tol < max( abs( ExDoppler - Doppler ) );
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
