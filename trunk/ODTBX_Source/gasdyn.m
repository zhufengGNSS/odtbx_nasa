function [xDot A Q] = gasdyn(t,x,options)
% GASDYN Simulates comet outgassing forces
%
% [xDot A Q] = gasdyn(t,x,options)
%
% GASDYN provides the forces acting on a spacecraft from comet outgassing.
% The model uses spherical harmonics and the position of the Sun with
% respect to the comet in order to generate resultant forces.
%
% INPUTS
% VARIABLE     SIZE              DESCRIPTION
%    t         1xM               Time steps for analysis
%    x         NxM               Initial 'N' states
%    options   structure         ODTBX options structure
%
% OUTPUTS
%    xDot      NxM               Time derivative of states
%    A         NxNxM             Jacobian
%    Q         NxNxM             Process noise spectral density
%
%
% keyword: dynamic, comet, outgas
% See also POLYDYN
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
%   REVISION HISTORY
%    Author                 Date         	    Comment
%    Benjamin Asher         08/09/2012          Original

%% Initials to be used when actually implemented into ODTBX
epoch = getOdtbxOptions(options,'epoch',[]);
area = getOdtbxOptions(options,'srpArea',[]);
mass = getOdtbxOptions(options,'mass',[]);
SPICE_ID = getOdtbxOptions(options,'SPICE',[]);

et = cspice_str2et(datestr(epoch + t./86400,...
    'dd mmm yyyy HH:MM:SS.FFF'));
Rc = cspice_spkezr(SPICE_ID, et, 'J2000', 'NONE', 'Sun');
lent = length(t);

% These are constants for a 4th degree and order spherical harmonics model
% derived from the comet Tempel-1 at a distance of 2 AU from the Sun.
alpha = [1 0 0 0 0;
    1.7424380 2.4954340e-1 0 0 0;
    8.0649220e-1 1.0767340e-1 1.2179830e-2 0 0;
    -1.6492650e-2 -4.0825580e-3 -1.5942440e-3 -3.3663670e-4 0;
    -6.4582310e-2 -9.6370230e-4 -2.0065580e-3 -9.4770940e-5 -1.2518440e-4
    ];
beta = [0 0 0 0 0;
    0 -3.2874030e-3 0 0 0;
    0 -1.7287340e-3 4.2847190e-4 0 0;
    0 3.1134650e-3 -7.3604510e-4 -1.3549240e-4 0;
    0 4.2770650e-3 -5.6703060e-4 -9.6860660e-5 2.1226330e-5
    ];

%% Create Sun pointing Coordinate System
% Sun Pointing Vector (N vector)
rc = -unit(Rc(1:3,:));             

nc = [Rc(2,:).*Rc(6,:)-Rc(3,:).*Rc(5,:)
    Rc(3,:).*Rc(4,:)-Rc(1,:).*Rc(6,:)
    Rc(1,:).*Rc(5,:)-Rc(2,:).*Rc(4,:)];

% Normal Vector to orbital Plane (B vector)
nc = unit(reshape(nc,size(rc)));    

vc = [rc(2,:).*nc(3,:)-rc(3,:).*nc(2,:)
    rc(3,:).*nc(1,:)-rc(1,:).*nc(3,:)
    rc(1,:).*nc(2,:)-rc(2,:).*nc(1,:)];

% Velocity direction of Comet (T vector)
vc = unit(reshape(vc,size(rc)));            
for ii = lent:-1:1
    DCM(:,:,ii) = [vc(:,ii)';nc(:,ii)';rc(:,ii)'];
    x(1:3,ii) = DCM(:,:,ii)*x(1:3,ii);
end
[xc r] = unit(x(1:3,:));
cone_ang = zeros(1,lent);
clock_ang = zeros(1,lent);
n = 4;                                      % degree and order
Pn = zeros(1,lent);
Pd = 6.960e-5;                              % at 2 AU
I = eye(3);
%% Find xDot, A, and Q
xDot = zeros(6,1,lent);
if nargout > 1
    A = zeros(6,6,lent);
    dPdr = zeros(1,lent);
    dPdcone = zeros(1,lent);
    dPdclock = zeros(1,lent);
end
for j = 1:lent
    %Projection of State (CCI) onto T-B plane
    V = x(1:3,j)-rc(:,j).*x(1:3,j);
    cone_ang(j) = acos(sum(conj(rc(:,j)).*...
        x(1:3,j))/(norm(rc(:,j))*norm(x(1:3,j))));
    clock_dot = sum(conj(vc(:,j)).*V)/(norm(vc(:,j))*norm(V));
    if sum(conj(nc(:,j)).*V) < 0
        clock_ang(j) = 2*pi-acos(clock_dot);
    else
        clock_ang(j) = acos(clock_dot);
    end
    if nargout > 1
        [L dL] = legendre_cos(n,cone_ang(j));
    else
        L = legendre_cos(n,cone_ang(j));
    end
    for ii = 1:n+1
        for jj = 1:ii
            Pn(j) = Pn(j)+Pd/(r(j)^2)*(1+L(jj,ii)*(alpha(ii,jj)*...
                cos((jj-1)*clock_ang(j))+beta(ii,jj)*sin((jj-1)...
                *clock_ang(j))));
            if nargout > 1
                dPdr(j) = dPdr(j)-2*Pd/r(j)^3*(1+L(jj,ii)*(alpha(ii,jj)...
                    *cos((jj-1)*clock_ang(j))+beta(ii,jj)*...
                    sin((jj-1)*clock_ang(j))));
                dPdcone(j) = dPdcone(j)+Pd/r(j)^2*dL(jj,ii)*...
                    (alpha(ii,jj)*cos((jj-1)*clock_ang(j))...
                    +beta(ii,jj)*sin((jj-1)*clock_ang(j)));
                dPdclock(j) = dPdclock(j)+Pd/r(j)^2*L(jj,ii)*(jj-1)*...
                    (-alpha(ii,jj)*sin((jj-1)*clock_ang(j))...
                    +beta(ii,jj)*cos((jj-1)*clock_ang(j)));
            end
        end
    end
    agas = DCM(:,:,j)'*xc(1:3,j).*Pn(j)*area/mass./1000;
    xDot(1:6,:,j) = [x(4:6,j); agas];
    if nargout > 1
        drx = x(1,j)/r(j);
        dry = x(2,j)/r(j);
        drz = x(3,j)/r(j);
        dphix = x(1,j)*x(3,j)/(r(j)^2*sqrt(x(1,j)^2+x(2,j)^2));
        dphiy = x(2,j)*x(3,j)/(r(j)^2*sqrt(x(1,j)^2+x(2,j)^2));
        dphiz = -sqrt(x(1,j)^2+x(2,j)^2)/r(j)^2;
        dlambdax = -x(2,j)^2/(x(1,j)^2+x(2,j)^2);
        dlambday = x(1,j)/(x(1,j)^2+x(2,j)^2);
        dlambdaz = 0;
        dr = [drx dry drz;dphix dphiy dphiz;dlambdax dlambday dlambdaz];
        drdr = 1/r(j)*(I-x(1:3,j)*x(1:3,j)'./r(j)^2);
        dPdR = [dPdr(j) dPdcone(j) dPdclock(j)]*dr;
        dPr = xc(:,j)*dPdR;
        dvdr = area/mass.*(dPr+drdr.*Pn(j))./1000;
        A(1:3,4:6,j) = I;
        A(4:6,1:3,j) = DCM(:,:,j)'*dvdr*DCM(:,:,j);
    end
end
if nargout == 3
    Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9]).^2,[1 1 lent]);
end
end


function [L dL] = legendre_cos(n,ang)
L = zeros(n+1);
if nargout > 1;
    dL = zeros(n+1);
end
for ii = 1:n+1
    for jj = 1:ii
        if ii == 1 && jj == 1
            L(jj,ii) = 1;
            if nargout > 1
                dL(jj,ii) = 0;
            end
        elseif ii == jj+1
            L(jj,ii) = (2*(jj-1)+1)*cos(ang)*L(jj,jj);
            if nargout > 1
                dL(jj,ii) = (2*(jj-1)+1)*cos(ang)*dL(jj,jj)-(2*(jj-1)+1)...
                    *sin(ang)*L(jj,jj);
            end
        elseif ii == jj
            L(jj,ii) = (2*(jj-1)-1)*sin(ang)*L(jj-1,jj-1);
            if nargout > 1
                dL(jj,ii) = (2*(jj-1)-1)*sin(ang)*dL(jj-1,jj-1)...
                    +(2*(jj-1)-1)*cos(ang)*L(jj-1,jj-1);
            end
        else
            L(jj,ii) = 1/(ii-jj)*((2*(ii-1)-1)*cos(ang)*L(jj,ii-1)...
                -(ii+jj-3)*L(jj,ii-2));
            if nargout > 1
                dL(jj,ii) = 1/(ii-jj)*((2*(ii-1)-1)*cos(ang)*dL(jj,ii-1)...
                    -(2*(ii-1)-1)*sin(ang)*L(jj,ii-1)...
                    -(ii+jj-3)*dL(jj,ii-2));
            end
        end
    end
end
end