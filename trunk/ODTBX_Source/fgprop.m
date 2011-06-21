function [r,v,Phi,Phinv] = fgprop(ro,vo,dt,mu)
% FGPROP Lagrangian two-body propagation.
%   [r,v] = fgprop(ro,vo,dt,mu) propagates the position and velocity vector
%   pair (ro,vo) over dt using Lagrange's coefficients f, g, fdot, & gdot,
%   computed using a Newton-Raphson iteration on the universal anomaly
%   based on Vallado's Algorithm 8.  If mu is ommitted, a default of
%   3.986004415e+14 m^3/sec^2 will be used.
%   [r,v,Phi] = fgprop(ro,vo,dt,mu) uses (r,v) and (ro,vo) to compute a
%   two-body STM according to Battin. The inverse of Phi is an optional
%   fourth output argument.
%   [r,v,Phi,Phinv] = fgprop(...) also returns the symplectically-computed
%   inverse of the STM.
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

%   Russell Carpenter, NASA GSFC

tol = 1e-8; % Vallado's suggestion
maxiter = 50; % Vallado's suggestion
if nargin < 4
    mu = 3.986004415e+14; % m^3/sec^2
end
if nargin == 0 % Self-test from Vallado
    ro = [1131.340;-2282.343;6672.423]*1e3;
    vo = [-5.64305; 4.30333;2.42879]*1e3;
    dt = 40*60;
end
% Precalculate frequently used variables:
rtmu = sqrt(mu);
Ro = norm(ro);
V = norm(vo);
rdotv = dot(ro,vo);
sig = rdotv/rtmu;
E = V^2/2 - mu/Ro;
alpha = -2*E/mu; % inverse sma (defined for all orbits)
% Find initial guess for Sundman's anomaly (chi):
if alpha == 0 % parabolic
    p = norm(cross(ro,vo))^2/mu;
    s = (pi/2 - atan(3*sqrt(mu/p^3)*dt))/2;
    w = atan(tan(s)^(1/3));
    chi = sqrt(p)*2*cot(2*w);
else
    a = 1/alpha;
    if alpha > 0 % elliptical (or circular)
        % Need to resolve dt for multiple revs, accounting for negative dt:
        Tp = 2*pi*sqrt(a^3/mu);
        % not needed: nrevs = sign(dt)*floor(abs(dt)/p);
        dt = rem(dt,Tp); % rem accounts for negative dt; mod does not
        chi = rtmu*dt*alpha;
        % Vallado's states that when using canonical units for case when
        % Ro = R_planet, which produces alpha = 1, the guess above is so
        % close to its true value that convergence is slow.  He suggests
        % the following:
        if abs(1-alpha) < tol
            chi = .97*chi;
        end
    elseif alpha < 0 % hyperbolic
        chi = sign(dt)*sqrt(-a)*log(-2*mu*alpha*dt/(rdotv ...
            + sign(dt)*sqrt(-mu*a)*(1 - Ro*alpha)));
    end
end
i = 0;
dchi = chi;
% Use Newton-Raphson to solve for Sundman's anomaly:
while abs(dchi) > tol
    if i > maxiter
        error('FGPROP:maxit','Max Iterations Reached in FGPROP.')
    end
    chi2 = chi^2;
    psi = chi2*alpha;
    [c2,c3] = cfuns(psi);
    R = chi2*c2 + sig*chi*(1 - psi*c3) + Ro*(1 - psi*c2);
    dchi = (rtmu*dt - chi^3*c3 - sig*chi2*c2 - Ro*chi*(1 - psi*c3))/R;
    chi = chi + dchi;
    i = i + 1;
end
% Use Lagrange's coefficients to advance position and velocity vectors:
f = 1 - chi^2/Ro*c2;
gdot = 1 - chi^2/R*c2;
g = dt - chi^3/rtmu*c3;
fdot = rtmu/R/Ro*chi*(psi*c3 - 1);
r = f*ro + g*vo;
v = fdot*ro + gdot*vo;
if nargout == 2
    return
end
% Battin uses universal functions in terms of sqrt(chi) for the STM
% computation.  These are related to his c-functions by chi^n*c_n = U_n.
U2 = chi^2*c2;
if alpha == 0
    U4 = chi^4/24;
    U5 = chi^5/120;
else
    U3 = chi^3*c3;
    % Use recursions, since we already have U2 and U3:
    U4 = ((chi^2/2) - U2)/alpha; 
    U5 = ((chi^3/6) - U3)/alpha;
end
C = (3*U5 - chi*U4 - rtmu*dt*U2)/rtmu;
% Battin's partitions of the STM:
I = eye(3);
Ro2 = Ro^2;
Ro3 = Ro2*Ro;
R2 = R^2; R3 = R2*R;
dv = v - vo;
dvdvt = dv*dv';
Rolmf = Ro*(1 - f);
Phi(4:6,4:6) = Ro/mu*dvdvt + (Rolmf*r*ro' - C*r*vo')/R3 + gdot*I;
Phi(1:3,1:3) =  R/mu*dvdvt + (Rolmf*r*ro' + C*v*ro')/Ro3 + f*I;
Phi(1:3,4:6) = Rolmf/mu*((r - ro)*vo' - dv*ro') + C/mu*v*vo' + g*I;
Phi(4:6,1:3) = -dv*ro'/Ro2 - r/R2*dv' + ...
    fdot*(I - r/R2*r' + (r*v' - v*r')*r*dv'/mu/R) - mu*C/R3/Ro3*r*ro';
Phinv = [Phi(4:6,4:6)', -Phi(1:3,4:6)'; - Phi(4:6,1:3)', Phi(1:3,1:3)'];

function [c2,c3] = cfuns(psi)
% CFUNS 3rd and 4th of Battin's c-functions.
%   Battin's "c functions" relate his universal anamoly to Sundman's
%   universal anomaly.  The 3rd and 4th functions, c2 and c3, are used in
%   solving Kepler's problem.  Battin's original work called these C and S,
%   respectively.  These are more generally known as Stumpff functions.  
tol = 1e-6; % Vallado's suggestion
abpsi = abs(psi);
if abpsi > tol
    rtpsi = sqrt(abpsi);
    rtps3 = sign(psi)*rtpsi^3;
    if psi > tol % elliptical case
        c2 = (1-cos(rtpsi))/psi;
        c3 = (rtpsi - sin(rtpsi))/rtps3;
    else % hyperbolic case
        c2 = (1-cosh(rtpsi))/psi;
        c3 = (rtpsi - sinh(rtpsi))/rtps3;
    end
else % Parabolic case
    c2 = 1/2;
    c3 = 1/6;
end
