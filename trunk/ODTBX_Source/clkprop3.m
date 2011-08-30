function [x,rtQd] = clkprop3(tspan,xo,opt)
% [x,rtQd] = clkprop3(tspan,xo,opt) is a 3-state clock error propagagor.
% Units are km, km/s, and km/s^2; set c = 1 for s, s/s, and s/s^2.
%
% See also: gpspseudomeas, clkdyn3

% REVISION HISTORY
%   Author      		Date         	Comment
%   Russell Carpenter   08/25/2011      Original

c = JATConstant('c')/1e3;
dt = diff(tspan);
if length(dt) > 1
    error('ODTBX:CLKPROP3:multStep',...
        'CLKPROP3 propagates one step at a time; length(tspan) should be 2')
end
if nargin == 2 || isempty(opt)
    h_0 = 2.4e-22;
    hm2 = 8.0e-28;
    q1 = c^2*h_0/2;
    q2 = c^2*2*pi^2*hm2;
    q3 = c^2*1e-35;
else
    q1 = c^2*opt.q1;
    q2 = c^2*opt.q2;
    q3 = c^2*opt.q3;
end
Phi = eye(3) + dt*diag([1 1],1) + dt^2/2*diag(1,2);
x = Phi*xo;
rtQd(3,3) = sqrt(q3*dt);
rtQd(2,3) = sqrt(q3*dt)*dt/2;
rtQd(2,2) = sqrt(q2*dt + q3*dt^3/12);
rtQd(1,3) = sqrt(q3*dt)*dt^2/6;
rtQd(1,2) = sqrt(q2*dt + q3*dt^3/12)*dt/2;
rtQd(1,1) = sqrt(q1*dt + q2*dt^3/12 + q3*dt^5/720);