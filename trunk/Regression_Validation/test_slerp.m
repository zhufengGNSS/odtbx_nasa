function fail = test_slerp()

% Regression test for slerp.m
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

fail = 0;

% two quaternions:
q0 = [0; 0; 0; 1]; % QED
q1 = [.2;.2;.2;.2]; % simple inputs, but properly corrected below
norm(q1);
q1= q1/norm(q1);

iv = [1;0;0]; % initial vector before transformation
t = 0:.01:1; % ratio between quaternions
v=zeros(length(t),3);
for tind = 1:length(t)
    q = slerp(q0,q1,t(tind));
    m = q2dcm(q);
    v(tind,:) = m*iv;
end

c0 = q2dcm(q0);
c1 = q2dcm(q1);

% figure;
% plot(t,v,'.')
% hold on;
% plot([0 0 0],c0*iv,'ro');
% plot([1 1 1]*t(end),c1*iv,'ro');
% hold off;

if norm(c0*iv - v(1,:)') > 1e-15
    fail = 1;
end
if norm(c1*iv - v(end,:)') > 1e-15
    fail = 1;
end
if norm([2; -1; 2]/3 - v(round(length(v)/2),:)') > 1e-15
    fail = 1;
end
