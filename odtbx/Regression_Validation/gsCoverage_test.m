function failed = gsCoverage_test()
% gsCoverage_test The regression test for gsCoverage
% 
% See also gsCoverage.m
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
%   Ravi Mathur         08/28/2012      Extract regression test from
%                                       original gsCoverage.m

% Define the initial conditions
x_kep.sma  = 7000; %km
x_kep.ecc  = 0.01; %unitless
x_kep.incl = 35*(pi/180); %radians
x_kep.raan = 250*(pi/180); %radians
x_kep.argp = 270*(pi/180); %radians
x_kep.tran = 100*(pi/180); %radians

muEarth   = 398600.4415; %km^3/s^2
x0 = kep2cart(x_kep,muEarth); %km & km/sec

EpochString = 'April 15, 2010 13:52:24';
epoch = datenum(EpochString);

% Define the dynamics
dynfun = @r2bp;
dynarg = muEarth;

% Define the data
gsList = createGroundStationList();
gsID   = {'HBKS','USHS','USPS','WHSX'};
gsECEF = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end

datfun = @gsmeas;
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS','WHSX'});
datarg = setOdtbxOptions(datarg,'rSigma',[1e-2 1e-5 1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsECEF',gsECEF);

% Run gsCoverage
[availTimes PossibleSched] = gsCoverage(dynfun,datfun,x0,0,10,60*60*2,[],dynarg,datarg);

% Define expected results
exp_availTimes = [0:10:310 590:10:1180 3080:10:3530 4350:10:4840 5980:10:6560 6770:10:7200];
exp_PossibleSched = ...
          [2           0         310;
           4         590        1180;
           1        3080        3530;
           3        4350        4840;
           2        5980        6560;
           4        6770        7200];
       
passed = all(all(PossibleSched==exp_PossibleSched)) & all(exp_availTimes==availTimes);
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
 
end