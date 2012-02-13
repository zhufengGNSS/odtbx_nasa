function failed = kepel_test()
% Test function for kepel.m.
%
% Tests one elliptical trajectory about the Earth with different argument
% presentation.  Didn't test passing in GM as an argument.
%
% From Vallad, David A., " Fundamentals of Astrodynamics and Applications",
% Third Ed., 2007, Example problem 2-5, pages 122-124.
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


% Inputs:
rijk = [6524.834; 6862.875; 6448.296] * 1000; % m
vijk = [4.901327; 5.533756; -1.976341] * 1000; % m/s

% Vallado example results:
a = 36127.343 * 1000; % semi-major axis, m
e = 0.832853; % eccentricity
i = 87.870 / 180 * pi; % inclination, rad
raan = 227.89 / 180 * pi; % right ascension of ascending node, rad
argp = 53.38 / 180 * pi; % argument of perigee, rad
truean = 92.335 / 180 * pi; % true anomaly ,rad

% scaled tolerance:
scaledtol = 1e-4;

% Test 1: standard 3x1 args
KOE1 = kepel(rijk,vijk);
failed = check(KOE1,a,e,i,raan,argp,truean,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 1'));
end

% Test 2: 1x3 args
KOE1 = kepel(rijk',vijk');
failed = check(KOE1,a,e,i,raan,argp,truean,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 2'));
end

% Test 3: mixed args
KOE1 = kepel(rijk',vijk);
failed = check(KOE1,a,e,i,raan,argp,truean,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 3'));
end

% Test 4: reversed mixed args
KOE1 = kepel(rijk,vijk');
failed = check(KOE1,a,e,i,raan,argp,truean,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 4'));
end

% Test 5: concatenated args
KOE1 = kepel([rijk; vijk]);
failed = check(KOE1,a,e,i,raan,argp,truean,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 5'));
end

% Test 6: output elements as separate variables instead of struct
[at,et,it,Wt,wt,ft] = kepel(rijk,vijk);
failed = check(KOE1,at,et,it,Wt,wt,ft,scaledtol);
if failed
    disp(sprintf('kepel_test.m failure occurrend in test 6'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support function to do the checks
function failed = check(KOE,a,e,i,raan,argp,truean,scaledtol)

failed = 0; % default to not failed

% test support function
t = (KOE.sma - a)/a;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED sma check, expected: %f, received: %f, scaled diff: %f',a,KOE.sma,t));
end
t = (KOE.ecc - e)/e;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED eccentricity check, expected: %f, received: %f, scaled diff: %f',e,KOE.ecc,t));
end
t = (KOE.incl - i)/i;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED inclination check, expected: %f, received: %f, scaled diff: %f',i,KOE.incl,t));
end
t = (KOE.raan - raan)/raan;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED RAAN check, expected: %f, received: %f, scaled diff: %f',raan,KOE.raan,t));
end
t = (KOE.argp - argp)/argp;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED ARGP check, expected: %f, received: %f, scaled diff: %f',argp,KOE.argp,t));
end
t = (KOE.tran - truean)/truean;
if t > scaledtol
    failed = 1;
    disp(sprintf('kepel_test.m FAILED TrueAn check, expected: %f, received: %f, scaled diff: %f',treuan,KOE.tran,t));
end