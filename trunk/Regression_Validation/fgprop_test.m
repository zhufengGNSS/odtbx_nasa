function fail = fgprop_test
% Regression test for FGPROP
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

mu = 3.986004415e+14;
dt = 40*60;
ro = [1131.340;-2282.343;6672.423]*1e3;
vo = [-5.64305; 4.30333;2.42879]*1e3;
orbtyp = {'elliptical','circular','hyperbolic'};
for k = 3:-1:1
    switch orbtyp{k}
        case 'hyperbolic'
            vo = 10*vo;
        case 'circular'
            Ro = norm(ro);
            ro = Ro*[1;0;0];
            vo = sqrt(mu/Ro)*[0;1;0];
        case 'parabolic'
            Ro = norm(ro);
            ro = Ro*[1;0;0];
            vo = sqrt(2*mu/Ro)*[0;1;0];
    end
    [r,v,Phi,Phinv] = fgprop(ro,vo,dt,mu);
    [~,X,PHI] = integ(@r2bp,[0 dt],[ro;vo],[],mu);
    dX = [r;v] - X(:,end);
    dPhi = Phi - PHI(:,:,end);
    pass(k) = norm(dPhi) < 1e-5 && norm(dX(1:3)) < 1e-3 && norm(dX(4:6)) < 1e-6 ...
        && norm(Phi*Phinv - eye(6)) < 1e-9;
end

% check corner case: dt=0
faildt0 = 0;
try
[r,v,Phi,Phinv] = fgprop(ro,vo,0,mu);
    if any(r ~= ro) || any(v ~= vo) || any(any(Phi ~= eye(6))) ...
            || any(any(Phinv ~= eye(6)))
        faildt0 = 1;
    end
catch %#ok<CTCH>
    warning('test_fgprop: unexpected error with dt=0 corner case');
    faildt0 = 1;  % unexpected error
end

fail = ~all(pass) || faildt0;