function [xdot A Q] = pancake_dyn(~,x,wP)

% Get state information
R = x(1:2,1);
V = x(3:4,1);
mu = x(5,1);
Rs = x(6:7,1);

r3 = norm(R)^3;

wPx = [0 -wP; wP 0];

% Calculate state derivatives
xdot = [V; 
        -mu/norm(R)^3*R; 
        0;
        wPx*Rs];

I = eye(2,2);
    
% Calculate partials
A = zeros(7,7);
A(6:7,6:7) = wPx; 
A(3:4,1:2) = -mu/r3*I + 3*mu*R*R'/r3/norm(R)^2; 
A(3:4,5) = -R/r3; 
A(1:2,3:4) = I;

% Set process noise spectral density
Q = zeros(7,7);

end