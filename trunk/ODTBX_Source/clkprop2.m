function [x,rtQd] = clkprop2(tspan,xo,opt)
% [x,rtQd] = clkprop2(tspan,xo,opt) is 2-state clock error propagator.
% Units are km, and km/s; set c = 1 for s, and s/s.
%
% See also: gpspseudomeas, clkdyn2

% REVISION HISTORY
%   Author      		Date         	Comment
%   Russell Carpenter   08/25/2011      Original

c = JATConstant('c')/1e3;
dt = diff(tspan);
if length(dt) > 1
    error('ODTBX:CLKPROP2:multStep',...
        'CLKPROP2 propagates one step at a time; length(tspan) should be 2')
end
if nargin == 2 || isempty(opt)
    h_0 = 2.4e-22;
    hm2 = 8.0e-28;
    q1 = c^2*h_0/2;
    q2 = c^2*2*pi^2*hm2;
else
    q1 = c^2*opt.q1;
    q2 = c^2*opt.q2;
end
Phi = eye(2) + dt*diag(1,1);
x = Phi*xo;
rtQd(2,2) = sqrt(q2*dt);
rtQd(1,2) = rtQd(2,2)*dt/2;
rtQd(1,1) = (sqrt(q1*dt + q2*dt^3/12));