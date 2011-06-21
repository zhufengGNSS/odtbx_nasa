function [y,H,R] = rrdotang(t,x1,x2,options)

% RRDOTANG  Calculate range, range rate, doppler, angles between two objects.
%
%   [y,H,R] = rrdotang(t,x1,x2,options)
% Returns the range and range rate (or doppler) measurement between two
% objects and is based on the GEONS DATSIM Math Spec. This function can
% also include the troposphere and ionosphere effects if requested.
%
% This function can return either '1way' or '2way' measurements based on
% the 'type' specified in the options structure. Since this measurement
% model is not including light time corrections or biases, it doesn't
% matter which object is transmitting and which is receiving. Also, two-way
% measurements are the same as one-way as far as y (the measurements) is
% concerned, but H (the measurement partial derivatives)is doubled for
% two-way.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES             NOTES
%   rangeType           {'1way'(default),'2way',
%                        '1wayFWD','1wayRTN'}    1way=1wayFWD
%   useRange            {true(default), false}
%   useRangeRate        {true(default), false}   can't combine with Doppler
%   useDoppler          {true, false(default)}   can't combine w/ RangeRate
%   useUnit             {true, false(default)}
%   useAngles           {true, false(default)}
%   frequencyTransmit   {scalar>0, 1.57542e9}    Hz, Only for Doppler and
%                                                measurement errors
%   useGPSIono          {true, false(default)}   only for GPS sats as x2
%   useIono             {true, false(default)}   only for groundstats as x2
%   useTropo            {true, false(default)}   only for groundstats as x2
%   useChargedParticle  {true, false(default)}   only for groundstats as x2
%   epoch                datenum                 UTC time associated with   
%                                                start of simulation.
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each type are grouped together in rows.
% Thus, using the default options settings, the output y will look like:
%   y = [  range(t1)      range(t2)      range(t3)      range(t4)     ... ;
%          rangeRate(t1)  rangeRate(t2)  rangeRate(t3)  rangeRate(t4) ... ]
%
% The angles data type returns an azimuth/elevation pair relative to the
% coordinates of the state vector, so that for example if the states are
% inertial, the angles would correspond to RA/Dec.
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t        (1xN)	Times corresponding to both x1 and x2 (secs)
%      x1       (6xN)   States (km & km/s) of user satellite at times t
%      x2       (6xN)   States (km & km/s) of tracking satellite/groundstat
%                       at times t
%      options  (1x1)   data structure
%
% OUTPUTS
%      y        (MxN)   measurements
%      H        (Mx6xN) measurement partials matrix
%      R        (MxMxN) measurement covariance
%
% VALIDATION TEST
%
%  To perform a validation test, replace the satPos variable with 
%  'ValidationTest' as the input argument.
%
% REGRESSION TEST
%
%  To perform a regression test, replace the Ephem structure with 
%  'RegressionTest' as the input argument.  
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, RRDOTLT, GSMEAS, 
%  ODTBXOPTIONS
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
%   Bo Naasz            12/03/2005      Original xlmeas.m
%   Derek Surka         07/20/2007      Modified xlmeas to work with OD
%                                       Toolbox
%   Derek Surka         08/22/2007      Modified input parameters, added
%                                       functionality
%   Derek Surka         09/07/2007      Added range rate capability
%   Derek Surka         09/13/2007      Added state outputs
%   Kevin Berry         04/08/2008      Simplified and renamed for user
%                                       friendlyness
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         06/25/2009      Added time scale comments
%   Kevin Berry         07/01/2009      Added angle measurement type
%   Kevin Berry         09/03/2009      Changed angle measurement to Unit
%   Kevin Berry         11/10/2009      Moved velocity calculations for 
%                                       position only measurement types
%   Russell Carpenter   02/10/2011      Added angles data type

%% Determine whether this is an actual call to the program or a test

if strcmpi(t,'ValidationTest')||strcmpi(t,'RegressionTest')
    y = rrdotang_validation_test();  
else
	[y,H,R] = Get_rrdotang(t,x1,x2,options);
end
end

%% Main function
function [y,H,R] = Get_rrdotang(t,x1,x2,options)

useRange        = getOdtbxOptions(options, 'useRange', true );
useRangeRate    = getOdtbxOptions(options, 'useRangeRate', true );
useDoppler      = getOdtbxOptions(options, 'useDoppler', false );
useUnit         = getOdtbxOptions(options, 'useUnit', false );
useAngles       = getOdtbxOptions(options, 'useAngles', false );

if( useDoppler && useRangeRate)
    error('rrdotang is not designed to handle both RangeRate and Doppler at the same time')
end

y = [];
if useRange
    r = LOSRange(t, x1(1:3,:), x2(1:3,:), options);
    y = [y; r];
end

if useRangeRate
    rdot = LOSRangeRate(t, x1(1:3,:), x1(4:6,:), x2(1:3,:), x2(4:6,:), options);
    y    = [y; rdot];
end

if useDoppler
    fd = LOSDoppler(t, x1(1:3,:), x1(4:6,:), x2(1:3,:), x2(4:6,:), options);
    y  = [y; fd];
end

if useUnit
    u = LOSUnit(t, x1(1:3,:), x2(1:3,:), options);
    y  = [y; u];
end

if useAngles
    a = [x1(1:2,:) - x2(1:2,:); zeros(size(x2(3,:)))];
    A = sqrt(sum(a.^2));
    az = atan2(a(2,:),a(1,:)); %real(acos(a(1,:)./A));
    el = atan2(x1(3,:) - x2(3,:),A); %real(acos(dot(x1(1:3,:) - x2(1:3,:),a)./r./A));
    zA = find(A==0);
    az(zA) = 0;
    el(zA) = 0;
    y = [y;az;el];
end

% Extra Outputs
if nargout > 1,
    H         = [];
    dr        = x1(1:3,:)-x2(1:3,:); %pos diff vector
    [udr rho] = unit(dr);
    if( useRange )
        Hrr = shiftdim(udr,-1);
        Hrv = zeros(size(Hrr));
        H   = [H; [Hrr Hrv] ];
    end
    if( useRangeRate )
        dv        = x1(4:6,:)-x2(4:6,:); %vel diff vector
        Hvr = shiftdim(dv.*repmat(1./rho,3,1) ...
              - dr.*repmat(dot(dr,dv)./rho.^3,3,1),-1);
          % cross(dr,cross(dr,dv))./rho.^3;
        Hvv = shiftdim(udr,-1);
        H   = [H; [Hvr Hvv] ];
    end
    if( useDoppler )
        dv        = x1(4:6,:)-x2(4:6,:); %vel diff vector
        Hvr = shiftdim(dv.*repmat(1./rho,3,1) ...
              - dr.*repmat(dot(dr,dv)./rho.^3,3,1),-1);
        Hvv = shiftdim(udr,-1);
        c = JATConstant('c')/1000;
        f = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
        H = [H; -f/c*[Hvr Hvv] ];
    end
    if ( useUnit )
        for n=size(x1,2):-1:1
            Hur(:,:,n) = (eye(3)-udr(:,n)*udr(:,n)')/rho(n);
        end
        Huv = zeros(size(Hur));
        H   = [H; [Hur Huv] ];
    end
    if useAngles
        b = cross(dr,repmat([0;0;1],[1 size(dr,2)]));
        rho = repmat(rho',[1 3]);
        A = repmat(A',[1 3]);
        Har(2,:,:) = permute(cross(b,dr)'./rho.^2./A,[3 2 1]);
        Har(1,:,:) = permute(-b'./A.^2,[3 2 1]);
        Har = squeeze(Har);
        Har(1:2,:,zA) = zeros(2,3,length(zA));
        Hav = zeros(size(Har));
        H = [H; [Har Hav] ];
    end
end
if nargout > 2,
    %get sigma out of the options
    nMTypes      = sum([useRange useRangeRate useDoppler]); %# of types
    sigmaDefault = ones(1,nMTypes)*1e-3; %default sigma
    sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
    sigma        = diag(sigma);
    
    R = repmat(sigma.^2,[1,1,length(t)]);
end

end
%% Validation Test

function failed = rrdotang_validation_test()

disp(' ')
disp('Performing Test....')
disp(' ')
fprintf('%12s%17s%19s%17s%22s%23s\n','Pos 1 (km)','Pos 2 (km)','Vel 1 (km/s)','Vel 2 (km/s)','Expected (km km/s)','Calculated (km km/s)')
fprintf('%s\n\n',char(ones(1,110)*'-'));

tol = 1e-7;
options = odtbxOptions('measurement');
options = setOdtbxOptions(options,'epoch',datenum('Jan 1 2006')); %UTC
options = setOdtbxOptions(options,'useGPSIonosphere',false);
options = setOdtbxOptions(options,'useIonosphere',false);
options = setOdtbxOptions(options,'useTroposphere',false);
options = setOdtbxOptions(options,'useChargedParticle',false);
options = setOdtbxOptions(options, 'useRange', true );
options = setOdtbxOptions(options, 'useRangeRate', true );
options = setOdtbxOptions(options, 'useDoppler', false );
options = setOdtbxOptions(options, 'useUnit', false );
options = setOdtbxOptions(options, 'useAngles', true );
t=(1:9)*60*60;
e1 = [1; 0; 0]; e2 = [0; 1; 0]; e3 = [0; 0; 1];
r1 = [repmat(e1,1,3) repmat(e2,1,3) repmat(e3,1,3)]*9;
r2 = repmat([e1 e2 e3],1,3)*5;    
v1 = repmat(e1,1,9);
v2 = repmat([e1 e2 e3],1,3);
x1 = [r1;v1]; x2 = [r2;v2];

ExMeas = [               4                         0                         0                         0
           10.295630140987          1.35979800452198        -0.507098504392337                         0
           10.295630140987          1.35979800452198                         0        -0.507098504392337
           10.295630140987                         0          2.07789483118723                         0
                         4         -1.00000333565208           1.5707963267949                         0
           10.295630140987         0.485642144472135           1.5707963267949        -0.507098504392337
           10.295630140987                         0          3.14159265358979          1.06369782240256
           10.295630140987         0.485642144472135          -1.5707963267949          1.06369782240256
                         4         -1.00000333565208                         0                         0]';
y = Get_rrdotang(t,x1,x2,options);

fprintf('%4.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %6.2f %4.2f %4.2f %11.6f %9.6f %11.6f %9.6f\n',...
    [r1; r2; v1; v2; ExMeas; y]);

failed = any(tol < max( abs( ExMeas - y ) ));
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
end
