function HermiteInterpolator_regression_data
% Regression unit test for interpolation in HermiteInterpolator
% Creates test datasets for HermiteInterpolator verification
% Author Stephen D Metcalfe, Emergent Space Technologies
% Date 6-13-2010
% NOTE no acceleration is used when only 1 time entry is passed
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Stephen D Metcalfe  06/13/2010   	Original
%   Stephen D Metcalfe  09/29/2010      Added revision history block.
%                                       Changed to test around each entry
%                                       to verify change per Kevin Berry
%                                       9/28. Added timestamp for when
%                                       generated. Added name of who
%                                       generated data. Removed negative
%                                       time steps.
%   Stephen D Metcalfe  10/04/2010      Fixed MLINs and removed test for
%                                       test.current >= test.start

% please remove your name after running to force the next person to insert theirs
WHOGEN =                                %Insert your name here

% Setup the test parameters
fprintf('Configuring HermiteInterpolator_regression_data ... ');

% set the flag determining if detailed output should be displayed
% if ~= 0 enables detailed output
DETAIL=0;

% set the flag to end on first failure 
ENDONFAILURE=1;

% set the number of dimensions to test (1,2 or 3 [x,y,z])
DIMS=3;

% set the initial position
P=[1.0,2000.0,300000.0];

% set the initial velocity
V=[4.0,5000.0,600000.0];

% set the acceletation factors used for more than 1 entry
A=[7.0,8000.0,900000.0];

% set the starting times to sequence through
START=[999999,999,0,-999,-999999];

% set the time step sizes to sequence through
STEPS=[0.9, 999.0, 999999.0];  % step sizes - NO NEGATIVES

% set the number of entries to sequence through
ENTRIES = [20,15,10,8,7,6,5,4,3,2,1];

% set the percentage error allowed
%       x      y      x
PERR=[0.0001,0.0001,0.0001,...    % position
      0.0001,0.0001,0.0001]';     % velocity

% adjust arrays for number of dimensions
P = P(:,1:DIMS);                  % initial position
V = V(:,1:DIMS);                  % initial velocity
A1=[0.0,0.0,0.0];                 % acceleration when only 1 entry
A1=A1(:,1:DIMS);
An=A(:,1:DIMS);                   % acceleration when more than 1 entry
                                  % errors for position & velocity
PERR=reshape([PERR(1:DIMS,:) PERR(4:3+DIMS,:)],DIMS*2,1);

test.neqns = DIMS*2;    % number of equations
testnum=1;
for s = 1:length(START)
	test.start = START(s);                  % starting time
    for st = 1:length(STEPS)
        test.step = STEPS(st);              % time step size
        O = [-test.step/4,0,test.step/4];   % offsets before, at & after
        for e = 1:length(ENTRIES)
            test.entries = ENTRIES(e);      % number of time entries
            
            if test.entries == 1            % if 1 entry 
                A=A1;                       %   then no acceleration
            else                            % if > 1 entry 
                A=An;                       %   then use acceleration
            end
            for c = 1:test.entries          % for each entry
                t = test.start+((c-1)*test.step); 
                for o = 1:3                 % for before, at & after offsets
                    test.current = t+O(o);  % calculate current time
                    %pre-allocate test arrays
                    test.x = zeros(test.entries,1);         % time array
                    test.f = zeros(test.entries,test.neqns);% state array (x,y,z,x',y',z')
                    test.expected = zeros(test.neqns,1);    % expected result

                    % create test dataset
                    test.acceleration = A;
                    for te = 1:test.entries                 % for each time entry
                        t = test.start+((te-1)*test.step);  % calculate time
                        test.x(te) = t;
                        for p = 1:DIMS                      % for each dimension (x,y,z)
                                                            % calculate position & velocity
                            test.f(te,p)      = P(p)+(V(p)*t)+(0.5*A(p)*(t*t));
                            test.f(te,p+DIMS) = V(p)+(A(p)*t);
                        end
                    end

                    % calculate expected result
                    ct = test.current;
                    for p = 1:DIMS                          % for each dimension (x,y,z)
                                                            % calculate position & velocity
                        test.expected(p,:)      = P(p)+(V(p)*ct)+(0.5*A(p)*(ct*ct));
                        test.expected(p+DIMS,:) = V(p)+A(p)*ct;
                    end

                    % save the test in the testset
                    % speed is not an issue as this is rarely run
                    testset(testnum)=test; %#ok<AGROW>
                    testnum=testnum+1;
                end
            end
        end
    end
end

INITIALPOSITION = P;                  % initial position
INITIALVELOCITY = V;                  % initial velocity

GENTIME = datestr(now);
% WHOGEN must be set by each user
save HermiteInterpolator_regression_data WHOGEN GENTIME DETAIL ENDONFAILURE DIMS PERR INITIALPOSITION INITIALVELOCITY testset '-v6'; %#ok<USENS>
fprintf('Done\n');