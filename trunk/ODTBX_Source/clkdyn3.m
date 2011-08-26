function [xd,A,Q] = clkdyn3(t,x,opt)
% [xd,A,Q] = clkdyn3(~,x,opt) is a 3-state clock error dynamics model.
% Units are km, km/s, and km/s^2; set c = 0 for s, s/s, and s/s^2.
%
% See also: gpspseudomeas, clkdyn2

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Russell Carpenter   08/25/2011      Original

c = JATConstant('c')/1e3;
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
A = diag([1 1],1);
xd = A*x;
Q = diag([q1,q2,q3]);
n = length(t);
A = repmat(A,[1 1 n]);
Q = repmat(Q,[1 1 n]);