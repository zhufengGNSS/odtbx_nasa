function fail = gpspseudomeas_test
% Regression test for gpspseudomeas

% Test data
epoch  = datenum('Jan 1 2006');
t = [0 30];
x = [...
    6.661728208278700e+003    6.660377860346100e+003
    4.974733893323415e+003    4.964405999396554e+003
    1.344395137958589e+004    1.351452628945648e+004
    -4.488572276416573e-001   -4.513746422337851e-001
    -3.441691907798407e-001   -3.443570059471385e-001
    2.352751484151581e+000    2.352242062707958e+000
    1e+0                      2e+0
    1e-3                      2e-3];

% Configuration
measOptions = odtbxOptions('measurement');
measOptions = setOdtbxOptions(measOptions,'epoch',epoch);
measOptions = setOdtbxOptions(measOptions,'useRange', true);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', false);
measOptions.clockStateIndices = [7;8];
warning('off','ODTBX:GPSMEAS:noBodyQuat')

% Test cases
[y1,H1,R1] = gpspseudomeas(t,x,measOptions);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', true);
[y2,H2,R2] = gpspseudomeas(t,x,measOptions);
measOptions = setOdtbxOptions(measOptions,'useRangeRate', false);
measOptions = setOdtbxOptions(measOptions,'useDoppler', true);
[y3,H3,R3] = gpspseudomeas(t,x,measOptions);
measOptions = setOdtbxOptions(measOptions,'useRange', false);
[y4,H4,R4] = gpspseudomeas(t,x,measOptions);

% Create regression archive; usually this should be commented out
% y1a = y1; H1a = H1; R1a = R1;
% y2a = y2; H2a = H2; R2a = R2;
% y3a = y3; H3a = H3; R3a = R3;
% y4a = y4; H4a = H4; R4a = R4;
% save gpspseudomeas_test y1a H1a R1a y2a H2a R2a y3a H3a R3a y4a H4a R4a

% Load regression archive and test against current values
load gpspseudomeas_data
dy1 = abs(y1 - y1a); dH1 = abs(H1 - H1a); dR1 = abs(R1 - R1a); 
dy2 = abs(y2 - y2a); dH2 = abs(H2 - H2a); dR2 = abs(R2 - R2a); 
dy3 = abs(y3 - y3a); dH3 = abs(H3 - H3a); dR3 = abs(R3 - R3a); 
dy4 = abs(y4 - y4a); dH4 = abs(H4 - H4a); dR4 = abs(R4 - R4a); 

ytol = 1e-10;
yfail = any(dy1>ytol) + any(dy2>ytol) + any(dy3>ytol) + any(dy4>ytol);
Htol = 1e-10;
Hfail = any(any(dH1>Htol)) + any(any(dH2>Htol)) + any(any(dH3>Htol)) ...
    + any(any(dH4>Htol));
Rtol = 1e-10;
Rfail = any(any(dR1>Rtol)) + any(any(dR2>Rtol)) + any(any(dR3>Rtol)) ...
    + any(any(dR4>Rtol));

fail = any(yfail + shiftdim(Hfail,1) + shiftdim(Rfail,1) > 0);