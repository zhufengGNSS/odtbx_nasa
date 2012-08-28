function fail = gpslinkbudget_test()
%
% Regression and unit test for gpslinkbudget.m
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur             08/27/2012      Rename to conform to new
%                                           regression test format 

% constants
d2r = pi/180; % degrees to radians
r2d = 180/pi; % degrees to radians
fail = 0;
GAIN_TOL = 2e-11; % 1e-12 for 1D, but 1e-11 for 2D interp

%% TEST 1: simplistic regression test - omni only
RX_link.Nf = 0;
RX_link.L = 0;
RX_link.freq = 1575.42e6;    % GPS L1 frequency, Hz
RX_link.Ts = 300; % K
RX_link.As = 0;
RX_link.Ae = 0;
TX_link.P_sv = 1; % reflects the truth.sv_power, xmit_power_offset, and block in gpsmeas

% a 1D omnidirectional antenna with no gain
omni_1d = zeros(9,2);
omni_1d(:,1) = [0 5 10 30 60 90 120 150 180]';

rx_el = (0:10:180)'*d2r;
rx_az = (0:20:360)'*d2r;
tx_el = rx_el; % totally fake, but numerically useful
tx_az = rx_az; % totally fake, but numerically useful
los_mag = ones(19,1)*20000;

[CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(los_mag, RX_link, TX_link, omni_1d, rx_el, rx_az, omni_1d, tx_el, tx_az, []) ;

C = JATConstant('c') / 1000;  % km/s Speed of light
expAd = 20.*log10((C/RX_link.freq)./(4*pi.*los_mag));
expCN0 = TX_link.P_sv + expAd(1) - (10*log10(RX_link.Ts)) + 228.6;

if any(abs(Ad-expAd) > GAIN_TOL)
    fprintf(1,'Ad calculation failure.\n');
    fail = 1;
end

if any(abs(CN0-expCN0) > GAIN_TOL)
    fprintf(1,'CN0 calculation failure.\n');
    fail = 1;
end

if any(Ar > GAIN_TOL)
    fprintf(1,'Ar calculation failure.\n');
    fail = 1;    
end

if any(At > GAIN_TOL)
    fprintf(1,'At calculation failure.\n');
    fail = 1;    
end

if any(abs(AP-RP) > GAIN_TOL)
    fprintf(1,'AP and RP mismatch.\n');
    fail = 1;    
end

% Test the TX gain calcs:
CN0_bias = 13;
[CN02, ~, At] = gpslinkbudget(los_mag, RX_link, TX_link, omni_1d, rx_el, rx_az, [], tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(At-CN0_bias) > GAIN_TOL)
    fprintf(1,'TX gain calculation failure.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure.\n');
    fail = 1;
end

% Test the RX gain calcs:
[CN02, Ar] = gpslinkbudget(los_mag, RX_link, TX_link, [], rx_el, rx_az, omni_1d, tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(Ar-CN0_bias) > GAIN_TOL)
    fprintf(1,'RX gain calculation failure.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure.\n');
    fail = 1;
end

%% TEST 2: check input handling for 1D cases

% los_mag, RX_link, TX_link, omni_1d, rx_el, rx_az, omni_1d, tx_el, tx_az, CN0
arglist = [...
    1       1           1       1       1       1       1       1       1   1;... % too many args
    1       1           1       1       1       1       1       0       1   0;...
    1       1           1       1       1       1       0       1       1   0;...
    1       1           1       1       0       1       1       1       1   0;...
    1       1           1       0       1       1       1       1       1   0;...
    1       1           0       1       1       1       1       1       1   0;...
    1       0           1       1       1       1       1       1       1   0;...
    0       1           1       1       1       1       1       1       1   0;...
    ];

for i = 1:size(arglist,1)
    if arglist(i,1); arg1=los_mag; else arg1=[];end
    if arglist(i,2); arg2=RX_link; else arg2=[];end
    if arglist(i,3); arg3=TX_link; else arg3=[];end
    if arglist(i,4); arg4=omni_1d; else arg4=[];end
    if arglist(i,5); arg5=rx_el; else arg5=[];end
    if arglist(i,6); arg6=rx_az; else arg6=[];end
    if arglist(i,7); arg7=omni_1d; else arg7=[];end
    if arglist(i,8); arg8=tx_el; else arg8=[];end
    if arglist(i,9); arg9=tx_az; else arg9=[];end
    if arglist(i,10); arg10=CN0; else arg10=[];end
    try
        [CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) ; %#ok<NASGU>
        fail = 1;
    catch ex %#ok<NASGU>
        % expected
    end
end

% cases we expect to pass with no error, since az isn't used with 1D
% antennas:
% los_mag, RX_link, TX_link, omni_1d, rx_el, rx_az, omni_1d, tx_el, tx_az, CN0
arglist = [...
    1       1           1       1       1       1       1       1       0   0;...
    1       1           1       1       1       0       1       1       1   0;...
    ];

for i = 1:size(arglist,1)
    if arglist(i,1); arg1=los_mag; else arg1=[];end
    if arglist(i,2); arg2=RX_link; else arg2=[];end
    if arglist(i,3); arg3=TX_link; else arg3=[];end
    if arglist(i,4); arg4=omni_1d; else arg4=[];end
    if arglist(i,5); arg5=rx_el; else arg5=[];end
    if arglist(i,6); arg6=rx_az; else arg6=[];end
    if arglist(i,7); arg7=omni_1d; else arg7=[];end
    if arglist(i,8); arg8=tx_el; else arg8=[];end
    if arglist(i,9); arg9=tx_az; else arg9=[];end
    if arglist(i,10); arg10=CN0; else arg10=[];end
    try
        [CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) ; %#ok<NASGU>
        % expected
    catch ex %#ok<NASGU>
        fail = 1;
    end
end


%% 2D omni antennas
omni_2d = zeros(10,10);
omni_2d(1,2:10) = 0:45:360;
omni_2d(2:10,1) = 0:(45/2):180;

%% TEST 3: check input handling (2D)
% los_mag, RX_link, TX_link, omni_1d, rx_el, rx_az, omni_1d, tx_el, tx_az, CN0
arglist = [...
    1       1           1       1       1       1       1       1       1   1;... % too many args
    1       1           1       1       1       1       1       1       0   0;...
    1       1           1       1       1       1       1       0       1   0;...
    1       1           1       1       1       1       0       1       1   0;...
    1       1           1       1       1       0       1       1       1   0;...
    1       1           1       1       0       1       1       1       1   0;...
    1       1           1       0       1       1       1       1       1   0;...
    1       1           0       1       1       1       1       1       1   0;...
    1       0           1       1       1       1       1       1       1   0;...
    0       1           1       1       1       1       1       1       1   0;...
    ];

for i = 1:size(arglist,1)
    if arglist(i,1); arg1=los_mag; else arg1=[];end
    if arglist(i,2); arg2=RX_link; else arg2=[];end
    if arglist(i,3); arg3=TX_link; else arg3=[];end
    if arglist(i,4); arg4=omni_2d; else arg4=[];end
    if arglist(i,5); arg5=rx_el; else arg5=[];end
    if arglist(i,6); arg6=rx_az; else arg6=[];end
    if arglist(i,7); arg7=omni_2d; else arg7=[];end
    if arglist(i,8); arg8=tx_el; else arg8=[];end
    if arglist(i,9); arg9=tx_az; else arg9=[];end
    if arglist(i,10); arg10=CN0; else arg10=[];end
    try
        [CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) ; %#ok<NASGU>
        fail = 1;
    catch ex %#ok<NASGU>
        % expected
    end
end


%% TEST 4: Check vs 1D values

% check vs 1D antennas
% an antenna where the gain matches the angle exactly
real_1d = zeros(9,2);
real_1d(:,1) = [0 5 10 30 60 90 120 150 180]';
real_1d(:,2) = [0 5 10 30 60 90 120 150 180]';

[CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(los_mag, RX_link, TX_link, real_1d, rx_el, rx_az, real_1d, tx_el, tx_az, []) ;

% Ar & At are the same as the angles in deg
rx_el = (0:10:180)'*d2r;
if any(abs(Ar-(rx_el*r2d)) > GAIN_TOL)
    fprintf(1,'Ar calculation failure in test4.\n');
    fail = 1;
end
if any(abs(At-(rx_el*r2d)) > GAIN_TOL)
    fprintf(1,'At calculation failure in test4.\n');
    fail = 1;
end
if any(abs(CN0-(expCN0 + 2*rx_el*r2d)) > GAIN_TOL)
    fprintf(1,'CN0 calculation failure in test4.\n');
    fail = 1;
end
if any(abs(Ad-expAd) > GAIN_TOL)
    fprintf(1,'Ad calculation failure in test4.\n');
    fail = 1;
end
if any(abs(AP(1)-RP(1)) > GAIN_TOL)
    fprintf(1,'AP and RP mismatch in test4.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-(CN0-CN0(1))) > GAIN_TOL)
    fprintf(1,'RP mismatch from CN0 in test4.\n');
    fail = 1;    
end
if any(abs(AP-AP(1)-((tx_el-tx_el(1))*r2d)) > GAIN_TOL)
    fprintf(1,'AP mismatch from At in test4.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-((tx_el-tx_el(1))*r2d)-((rx_el-rx_el(1))*r2d)) > GAIN_TOL)
    fprintf(1,'RP mismatch from At & Ar in test4.\n');
    fail = 1;    
end

%% TEST 5: Check vs 2D values

% check vs 2D antennas, full coverage
% an antenna where the gain matches the angle exactly
real_2d = zeros(10,10);
real_2d(1,2:end) = ([0 5 10 30 60 90 120 150 180]*2)'; %az
real_2d(2:end,1) = [0 5 10 30 60 90 120 150 180]; % el
for i = 2:10
    for j = 2:10
        real_2d(i,j) = real_2d(1,i) + real_2d(j,1);
    end
end

[CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(los_mag, RX_link, TX_link, real_2d, rx_el, rx_az, real_2d, tx_el, tx_az, []) ;

Ar_exp = (rx_el+rx_az)*r2d;
At_exp = (tx_el+tx_az)*r2d;

% Ar & At are the same as the angles in deg
if any(abs(Ar-Ar_exp) > GAIN_TOL)
    fprintf(1,'Ar calculation failure in test5.\n');
    fail = 1;
end
if any(abs(At-At_exp) > GAIN_TOL)
    fprintf(1,'At calculation failure in test5.\n');
    fail = 1;
end
if any(abs(CN0-(expCN0 + Ar_exp + At_exp)) > GAIN_TOL)
    fprintf(1,'CN0 calculation failure in test5.\n');
    fail = 1;
end
if any(abs(Ad-expAd) > GAIN_TOL)
    fprintf(1,'Ad calculation failure in test5.\n');
    fail = 1;
end
if any(abs(AP(1)-RP(1)) > GAIN_TOL)
    fprintf(1,'AP and RP mismatch in test5.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-(CN0-CN0(1))) > GAIN_TOL)
    fprintf(1,'RP mismatch from CN0 in test5.\n');
    fail = 1;    
end
if any(abs(AP-AP(1)-(At_exp-At_exp(1))) > GAIN_TOL)
    fprintf(1,'AP mismatch from At in test5.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-(At_exp-At_exp(1))-(Ar_exp-Ar_exp(1))) > GAIN_TOL)
    fprintf(1,'RP mismatch from At & Ar in test5.\n');
    fail = 1;    
end

% Test the TX gain calcs:
CN0_bias = 13;
[CN02, ~, At] = gpslinkbudget(los_mag, RX_link, TX_link, real_2d, rx_el, rx_az, [], tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(At-At_exp-CN0_bias) > GAIN_TOL)
    fprintf(1,'TX gain calculation failure in test5.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure in test5.\n');
    fail = 1;
end

% Test the RX gain calcs:
[CN02, Ar] = gpslinkbudget(los_mag, RX_link, TX_link, [], rx_el, rx_az, real_2d, tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(Ar-Ar_exp-CN0_bias) > GAIN_TOL)
    fprintf(1,'RX gain calculation failure in test5.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure in test5.\n');
    fail = 1;
end

%% TEST 6: check masking of angles

% an antenna pattern that only extends 0-180 az, 0-90 el
cut_2d = zeros(7,7);
cut_2d(1,2:end) = ([0 5 10 30 60 90]*2)'; %az
cut_2d(2:end,1) = [0 5 10 30 60 90]; % el
for i = 2:7
    for j = 2:7
        cut_2d(i,j) = cut_2d(1,i) + cut_2d(j,1);
    end
end

[CN0, Ar, At, Ad, AP, RP] = gpslinkbudget(los_mag, RX_link, TX_link, cut_2d, rx_el, rx_az, cut_2d, tx_el, tx_az, []) ;

% any angles outside the antenna limits should be set to -100
Ar_exp = ones(length(rx_el),1)*(-100);
ind = rx_el <= pi/2;
Ar_exp(ind) = rx_el(ind)*r2d;
ind = (rx_az <= pi) & (rx_az ~= -100);
Ar_exp(ind) = Ar_exp(ind) + rx_az(ind)*r2d;

% any angles outside the antenna limits should be set to -100
At_exp = ones(length(tx_el),1)*(-100);
ind = tx_el <= pi/2;
At_exp(ind) = tx_el(ind)*r2d;
ind = (tx_az <= pi) & (tx_az ~= -100);
At_exp(ind) = At_exp(ind) + tx_az(ind)*r2d;

% Ar & At are the same as the angles in deg
if any(abs(Ar-Ar_exp) > GAIN_TOL)
    fprintf(1,'Ar calculation failure in test6.\n');
    fail = 1;
end
if any(abs(At-At_exp) > GAIN_TOL)
    fprintf(1,'At calculation failure in test6.\n');
    fail = 1;
end
if any(abs(CN0-(expCN0 + Ar_exp + At_exp)) > GAIN_TOL)
    fprintf(1,'CN0 calculation failure in test6.\n');
    fail = 1;
end
if any(abs(Ad-expAd) > GAIN_TOL)
    fprintf(1,'Ad calculation failure in test6.\n');
    fail = 1;
end
if any(abs(AP(1)-RP(1)) > GAIN_TOL)
    fprintf(1,'AP and RP mismatch in test6.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-(CN0-CN0(1))) > GAIN_TOL)
    fprintf(1,'RP mismatch from CN0 in test6.\n');
    fail = 1;    
end
if any(abs(AP-AP(1)-(At_exp-At_exp(1))) > GAIN_TOL)
    fprintf(1,'AP mismatch from At in test6.\n');
    fail = 1;    
end
if any(abs(RP-RP(1)-(At_exp-At_exp(1))-(Ar_exp-Ar_exp(1))) > GAIN_TOL)
    fprintf(1,'RP mismatch from At & Ar in test6.\n');
    fail = 1;    
end


% Test the TX gain calcs:
CN0_bias = 13;
[CN02, ~, At] = gpslinkbudget(los_mag, RX_link, TX_link, cut_2d, rx_el, rx_az, [], tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(At-At_exp-CN0_bias) > GAIN_TOL)
    fprintf(1,'TX gain calculation failure in test6.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure in test6.\n');
    fail = 1;
end

% Test the RX gain calcs:
[CN02, Ar] = gpslinkbudget(los_mag, RX_link, TX_link, [], rx_el, rx_az, cut_2d, tx_el, tx_az, (CN0+CN0_bias)) ;
if any(abs(Ar-Ar_exp-CN0_bias) > GAIN_TOL)
    fprintf(1,'RX gain calculation failure in test6.\n');
    fail = 1;
end
if any(abs(CN02-(CN0+CN0_bias)) > GAIN_TOL)
    fprintf(1,'CNO output difference during gain calculation failure in test6.\n');
    fail = 1;
end