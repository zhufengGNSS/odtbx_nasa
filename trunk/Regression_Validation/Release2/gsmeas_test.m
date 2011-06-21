function failed = gsmeas_test()
% Regression Test Case
% Function(s) gsmeas
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

failed  = 0;
tol     = 1e-8;

x1 = [6.661728208278700e+004    6.661278931448849e+004    6.660828815452385e+004    6.660377860346100e+004
    4.974733893323415e+003    4.971291888137325e+003    4.967849256683786e+003    4.964405999396554e+003
    1.344395137958589e+004    1.346747804712036e+004    1.349100301806504e+004    1.351452628945648e+004
   -4.488572276416573e-001   -4.496964226153426e-001   -4.505355608480570e-001   -4.513746422337851e-001
   -3.441691907798407e-001   -3.442318392101571e-001   -3.442944442685567e-001   -3.443570059471385e-001
    2.352751484151581e+000    2.352581973349160e+000    2.352412166194112e+000    2.352242062707958e+000];

x1(4:6,:) = zeros(3,size(x1,2));

v1.state_j2k.mjd    = datenum('Jan 1 2006');
v1.state_j2k.tSim   = [0 10 20 30];
v1.state_j2k.p      = x1(1:3,:);
v1.state_j2k.v      = x1(4:6,:);
v1.clk.b            = [0 1e-5 2e-5 3e-5];
v1.clock.b          = v1.clk.b;
v1.clock.tSim       = v1.state_j2k.tSim;

% Set up Ground Station information using the NDOSL
%--------------------------------------------------
epoch  = datenum('Jan 1 2006');
gsList = createGroundStationList();
gsID   = {  'DS16' 
            'DS46'
            'DS66' };

measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', true);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'gsList',gsList);
measOptions = setOdtbxOptions(measOptions,'gsID',gsID);
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',10);
measOptions = setOdtbxOptions(measOptions,'useIonosphere',false);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);

% First case
yTarget = [  
          62719.6215132007
        -0.142782836755335
                       NaN
                       NaN
                       NaN
                       NaN ];
kTargetNAN = [3;4;5;6];                   
tPos  = v1.state_j2k.tSim;
y = gsmeas(tPos(2), x1(:,2), measOptions);
PercentDiff = abs((yTarget(1:2)-y(1:2))./yTarget(1:2)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end

% Second test case
yTarget = [  
          676258.835263179
         0.184566269669634
          679959.450585763
        -0.340336897208893
                       NaN
                       NaN ];
measOptions = setOdtbxOptions(measOptions,'epoch',epoch+12000/86400);
y = gsmeas(tPos(2), x1(:,2).*[10;10;10;1;1;1], measOptions);
kTargetNAN = [5;6];     
PercentDiff = abs((yTarget(1:4)-y(1:4))./yTarget(1:4)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end

% Third test case
yTarget = [
                       NaN
                       NaN
                       NaN
                       NaN
          678183.274193485
        -0.299774814631325 ];
measOptions = setOdtbxOptions(measOptions,'epoch',epoch+50000/86400);
y = gsmeas(tPos(2), x1(:,2).*[10;10;10;1;1;1], measOptions);
kTargetNAN = [1;2;3;4];                   
PercentDiff = abs((yTarget(5:6)-y(5:6))./yTarget(5:6)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end


% Use Pre-Computed Ground Station information
%--------------------------------------------
gsECEF = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end

measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);

% 4th case
yTarget = [  
          62719.6215132007
        -0.142782836755335
                       NaN
                       NaN
                       NaN
                       NaN ];
kTargetNAN = [3;4;5;6];                   
tPos  = v1.state_j2k.tSim;
y = gsmeas(tPos(2), x1(:,2), measOptions);
PercentDiff = abs((yTarget(1:2)-y(1:2))./yTarget(1:2)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end

% 5th test case
yTarget = [  
          676258.835263179
         0.184566269669634
          679959.450585763
        -0.340336897208893
                       NaN
                       NaN ];
measOptions = setOdtbxOptions(measOptions,'epoch',epoch+12000/86400);
y = gsmeas(tPos(2), x1(:,2).*[10;10;10;1;1;1], measOptions);
kTargetNAN = [5;6];     
PercentDiff = abs((yTarget(1:4)-y(1:4))./yTarget(1:4)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end

% 6th test case
yTarget = [
                       NaN
                       NaN
                       NaN
                       NaN
          678183.274193485
        -0.299774814631325 ];
measOptions = setOdtbxOptions(measOptions,'epoch',epoch+50000/86400);
y = gsmeas(tPos(2), x1(:,2).*[10;10;10;1;1;1], measOptions);
kTargetNAN = [1;2;3;4];                   
PercentDiff = abs((yTarget(5:6)-y(5:6))./yTarget(5:6)) %#ok<NOPRT>
if( any(~isnan(y(kTargetNAN))) || max(PercentDiff>tol) )
    failed = 1;
end
