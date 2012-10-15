function failed = elevation_test()
% Regression Test Case
% Function(s) elevation
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
tol     = 1e-2;     

nTests = 5;
r1 = randn(3,nTests)*6500*1000;  % use m for JAT
r2 = randn(3,nTests)*8000*1000;

elOdtbx = elevation(r1,r2);

diffGS  = zeros(1,nTests);

for k=1:size(r1,2)
    gs          = jat.groundstations.GroundStation('elTest',r1(:,k));
    azelJat     = gs.azEl(r2(:,k)-r1(:,k));
    elJat       = azelJat(2);
    diffGS(k)   = elOdtbx(k) - elJat;
end

if( any( abs(diffGS) > tol ) )
    failed = 1;
end

