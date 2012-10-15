function failed = kepprop2b_test()
%
% kepprop2b_test Regression test for kepprop2b
% See also: kepprop2b.m
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
%   Ravi Mathur         08/28/2012      Extracted from kepprop2b.m

disp(' ')
disp(' ')
disp('Performing Test....')
disp(' ')

tol = 1e-7;

KOE.epoch = datenum([2008  9 26  0  0   .000]);
KOE.sma  = 42165.3431351313; %semi-major axis
KOE.ecc  = 0.26248354351331; %eccentricity
KOE.incl = 0.30281462522101; %inclination
KOE.raan = 4.85569272927819; %right ascension of the ascending node
KOE.argp = 4.71463172847351; %argument of periapse
KOE.tran = 2.37248926702153; %true anomaly
t        = 0:8640:86400; %secs from epoch of KOE
GM       = 398600.4415; %Gravitational constant

% Set Expected Results
KOE_exp = [

          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          2.37248926702153;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          2.79722436211144;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          3.18337407409023;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          3.57400974200765;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          4.01425565759545;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          4.57232665706546;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          5.35956850972672;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351         0.137251905665217;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          1.14521863765007;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          1.86433634881636;
          42165.3431351313          0.26248354351331          0.30281462522101          4.85569272927819          4.71463172847351          2.38486787064101]';

% Display Expected Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Expected Keplerian Elements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','sma','ecc','incl','raan','argp','tran');
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','(s)','(km)','(unitless)','(rad)','(rad)','(rad)','(rad)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:11
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),KOE_exp(:,n))
end         
disp(' ')
disp(' ')

KOEf = kepprop2b(KOE,t,GM);
KOE_calc = [KOEf.sma' KOEf.ecc' KOEf.incl' KOEf.raan' KOEf.argp' KOEf.tran']';

% Display Results
fprintf('%s\n',char(ones(1,96)*'-'));
disp('Calculated Keplerian Elements by time:')
fprintf('%s\n',char(ones(1,96)*'-'));
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','Time','sma','ecc','incl','raan','argp','tran');
fprintf('%-6s %14s %14s %14s %14s %14s %14s\n','(s)','(km)','(unitless)','(rad)','(rad)','(rad)','(rad)');
fprintf('%s\n',char(ones(1,96)*'-'));
for n=1:11
    fprintf('%-6i %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f\n',t(n),KOE_calc(:,n))
end         
disp(' ')

passed = tol > max(max(abs(KOE_calc-KOE_exp)));
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end

end
