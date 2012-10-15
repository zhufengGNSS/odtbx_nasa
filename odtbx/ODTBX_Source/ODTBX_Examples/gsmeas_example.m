%GSMEAS_EXAMPLE This demonstrates the use of the gsmeas measurement model.
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

close all

t = [0 60 120 180];
%x = [18000 10000 10000 0 0 0]';

%Over Goldstone
x = [5158.712338     -1514.921648      3777.724544       2.754779       9.375862      -0.001963;
5310.829717      -949.016866      3768.064835       2.313986       9.479618      -0.319662;
5436.235143      -378.334492      3739.451984       1.865419       9.535073      -0.633026;
5534.648919       194.234583      3692.270615       1.415291       9.542710      -0.937961]';

% Set up Ground Station information
%----------------------------------
epoch  = datenum('Jan 1 2006');
gsList = createGroundStationList();
gsID   = {  'DS12'
            'DS16' 
            'DS46'
            'DS66' };
nGS    = length(gsID);
gsECEF = zeros(3,nGS);
for n=1:nGS
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end

measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
measOptions = setOdtbxOptions(measOptions,'frequencyTransmit', 2106406.250);
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',0);
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);

measOptions = setOdtbxOptions(measOptions,'useLightTime',true);
measOptions = setOdtbxOptions(measOptions,'useIonosphere',false);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);
[y1,H1,R1] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useIonosphere',true);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',true);
[y2,H2,R2] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useLightTime',false);
measOptions = setOdtbxOptions(measOptions,'useIonosphere',false);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',false);
[y3,H3,R3] = gsmeas(t, x, measOptions);

measOptions = setOdtbxOptions(measOptions,'useIonosphere',true);
measOptions = setOdtbxOptions(measOptions,'useTroposphere',true);
[y4,H4,R4] = gsmeas(t, x, measOptions);

disp(' ')
disp(' ')
disp([y1, y2, y3, y4])
