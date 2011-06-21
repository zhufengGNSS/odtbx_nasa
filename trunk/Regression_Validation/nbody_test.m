function fail = nbody_test
% NBODY_TEST  N-Body force model test.
%
%   [fail] = NBODY_TEST regression test to verify results from the N-Body
%   Force Model. The test compares data from a simple Heliocentric orbit
%   example.
%
%   keyword: N-Body, Initialize
%
%   See also
%       N-Body:          NBODYPM
%       Initialize:      INITIALIZE_NBODY
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

% Benjamin Asher
% Embry-Riddle Aeronautical University
%
% Modification History
% ---------------------
% Author                  Date         	Comment
% Benjamin Asher          05/10/2011    Original nbody_test.m

fail = 0;

CentralBody = 'SUN';
PointMasses = {'EARTH BARYCENTER','MARS BARYCENTER',...
    'JUPITER BARYCENTER','SATURN BARYCENTER'};
epoch = datenum('30 MAY 2012 12:00:00');

nbodyopt = odtbxOptions('force');
nbodyopt = setOdtbxOptions(nbodyopt, 'epoch', epoch);
nbodyopt = setOdtbxOptions(nbodyopt, 'CentralBody', CentralBody);
nbodyopt = setOdtbxOptions(nbodyopt, 'PointMasses', PointMasses);
nbodyopt = initialize_nbody(nbodyopt);

KOE.sma = 149.9e6;
KOE.ecc = 0.2;
KOE.incl = 0.0863;
KOE.raan = 0;
KOE.argp = pi/4;
KOE.tran = 0;

tspan = [0 600*86400];
Xo = kep2cart(KOE,nbodyopt.GM_CB);

[t,X] = integ(@nbodypm,tspan,Xo,[],nbodyopt);

% t_test = t;
% X_test = X;
% 
% save('nbodyData.mat','*test')

plot3(X(1,:),X(2,:),X(3,:))
axis equal

load nbodyData.mat

dmax = 1e-5;

if(any(t-t_test))
    disp('N-Body Regression Test Failed: t')
    fail = 1;
end

if(max(max(abs(X-X_test)))>dmax)
    disp('N-Body Regression Test Failed: X')
    fail = 1;
end

end