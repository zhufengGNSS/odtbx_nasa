function Hpd = testMeasPartials(datfun,datarg)
% TESTMEASPARTIALS compares the numerically estimated H matrix to the H matrix provided by a measurement model.
%   HPD = TESTMEASPARTIALS(DATFUN,DATARG) uses the supplied measurement
%   data function DATFUN to compute measurements y and the measurement
%   partials matrix H for a random orbit using the options provided in
%   DATARG. The function then calls OMINUSC to numerically compute H using
%   NUMJAC. The function will display the median percent difference that
%   it finds between the resulting H of each method, and return the percent
%   difference of H for all time steps as HPD. If the median percent
%   difference is larger than 0.001%, the function will produce an error
%   message. For debugging purposes, "dbstop if error" can be called before
%   running TESTMEASPARTIALS to enter into Matlab's debug mode within the
%   function's workspace.
%
%   Example 
%       epoch  = datenum('Jan 1 2010');
%       gsList = createGroundStationList();
%       gsID   = {'DS16','DS46','DS66'};
%       gsECEF = zeros(3,length(gsID));
%       for n=1:length(gsID)
%           gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
%       end
%       measOptions = odtbxOptions('measurement');
%       measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
%       measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',-90); %Ignore the horizon
%       measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
%       Hpd = TestMeasPartials(@gsmeas,measOptions);
%
%   See also
%      Numerical Jacobian:    OMINUSC, NUMJAC
%      options handling:      ODTBXOPTIONS, SETODTBXOPTIONS,
%                             GETODTBXOPTIONS
%      Measurement Models:    GSMEAS, TDRSSMEAS, DDORMEAS, LNRMEAS
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
%   Kevin Berry         08/18/2010      Original
%   Kevin Berry         10/25/2011      Modified the code to work for
%                                       non-linear measurement models
%   Ravi Mathur         7/28/2015       Clear only estjac instead of all
%                                       functions, which also removes any
%                                       existing breakpoints

clear estjac % clears the persistent variables in estjac

disp(' ')
disp('Comparing the H matrix calculated in ')
disp(datfun);
disp(' with a numerically estimated H matrix using estjac (inside ominusc)')
disp(' ')

%% Set the random initial state for the satellite
R_apsis1  = 80000*rand+500+JATConstant('meanRadius','Earth')/1000; %between LEO and twice GEO
R_apsis2  = 80000*rand+500+JATConstant('meanRadius','Earth')/1000; %between LEO and twice GEO
kep1.sma  = (R_apsis1+R_apsis2)/2;
kep1.ecc  = abs(R_apsis1-R_apsis2)/(R_apsis1+R_apsis2);
kep1.incl = pi*rand;
kep1.raan = 2*pi*rand;
kep1.argp = 2*pi*rand;
kep1.tran = 2*pi*rand;
gm = JATConstant('muEarth')/1e9;
x0 = kep2cart(kep1,gm);

%% Propagate the test orbit for 1 period
T = 2*pi*sqrt(kep1.sma^3/gm);
[tspan,Xref] = integ(@r2bp,[0 T],x0,[],gm);

%% Get the measurement and partials matrix from the data function
[Y,Href] = feval(datfun,tspan,Xref,datarg);

%% Get the partials numerically from the estjac function within ominusc
eOpts = odtbxOptions('estimator');
eOpts = setOdtbxOptions(eOpts,'DatJTolerance',median(abs(Xref')));
eOpts = setOdtbxOptions(eOpts,'DatVectorized',2);
[~,Hcheck] = ominusc(@EmptyH,tspan,Xref,Y,eOpts,[],{datfun, datarg});

%% Difference the H matrices and display the maximum difference
tol = 1e-10; %tolerance for H values too close to zero to compare
Hdiff      = abs(Href - Hcheck);
isel_Hdiff = ~isinf(Hdiff)&~isnan(Hdiff)&abs(Href.*Hcheck)>tol; %Infs and NaNs show up when dividing by zero or 0/0

Hpd              = 100*Hdiff./max(abs(Href),abs(Hcheck));
Hpd(~isel_Hdiff) = 0;

disp(['The median percent difference in H is ',num2str(median(median(median(Hpd))),15)])
if median(median(median(Hpd)))>0.001
    error('The median percent difference in H is larger than 0.001%')
end

function [y,H,R] = EmptyH(tspan,Xref,options)
% This function calls datfun but empties out the H matrix output
datfun = options{1};
datarg = options{2};
[y,~,R]  = feval(datfun,tspan,Xref,datarg);
H=[];