function [xdot A Q] = polydyn(t,x,options)
% Force model for exterior gravitation of a polyhedron.

% POLYDYN
% Based on the algorithm presented in "Exterior Gravitation of a Polyhedron Derived
% and Compared with Harmonic and Mascon Gravitation Representations of
% Asteroid 4769 Castalia," R. Werner & D. Scheeres, 1996.
%
% xdot = POLYDYN(t,x,options) returns the external gravitation of a
% polyhedron at the field point defined in the state vector, x at time t.
%
% [xdot,A] = POLYDYN(t,x,options) also returns the Jacobian matrix, A.
%
% [xdot,A,Q] = POLYDYN(t,x,options) also returns a process noise power 
% spectral density based on a random walk process noise input.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1x1)	simulation time (secs from epoch)
%      x            (6x1)   Body-centered spacecraft state [pos;vel]
%      options      (1x1)   data structure, see below
%
%   OUTPUTS
%      xdot         (6x1)   derivatives
%      A            (6x6)   jacobian matrix
%      Q            (6x6)   process noise power spectral
%                           density
%
% OPTIONS is an options data structure. The options parameters that are 
% valid for this function are:
%
%   PARAMETER           VALID VALUES            NOTES
%   TR                   TriRep object          See TRIREP
%   RHO                  Scalar>0               Constant density
%   RA                   0< RA <2pi RAD         Right Ascension
%   DEC                  -pi/2< DEC <pi/2 RAD   Declination
%   PRA                  0< PRA <2pi RAD        Prime Meridian Angle
%   W                    Scalar [RAD/SEC]       Angular Rate
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
%
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
%   Kenneth Getzandanner    04/28/2012      Original polydyn.m
%

%% Initialize Variables

% TriRep object
tr = options.tr;

% Density
rho = options.rho;

r = x(1:3,1);
v = x(4:6,1);

% Universal gravitational constant
G = 6.67300e-20;

sumEdge = 0;
sumFace = 0;
sumEdgeA = 0;
sumFaceA = 0;

sumWf = 0;

% Calculate face normals, edges, and vertices
fn = faceNormals(tr);
ic = incenters(tr);
E = edges(tr);
V = tr.X;

% Calculate the DCM from the inertial to the body frame
w = options.w;
PRA = options.PRA;
RA = options.RA;
DEC = options.DEC;

theta = PRA + w*t;

D3w = [cos(-theta) sin(-theta) 0;
      -sin(-theta) cos(-theta) 0;
      0 0 1];
  
D1 = [1 0 0;
      0 cos(DEC-pi/2) sin(DEC-pi/2);
      0 -sin(DEC-pi/2) cos(DEC-pi/2)];
  
D3 = [cos(-pi/2-RA) sin(-pi/2-RA) 0;
      -sin(-pi/2-RA) cos(-pi/2-RA) 0;
      0 0 1];
  
C_IB = D3*D1*D3w;
D = C_IB';

r = D*r;

%% Edges
% Sum the gravitational contributions of each edge
for i=size(E,1):-1:1
    
    % Edge unit vectors
    n12 = (V(E(i,2),:)-V(E(i,1),:))'/norm(V(E(i,2),:)-V(E(i,1),:));
    n21 = -n12;
    
    % Find attached faces
    EA = edgeAttachments(tr,E(i,:));
    fni = fn(EA{1},:);
    ici = ic(EA{1},:);
    
    % Get outward face normals
    na = fni(1,:)';
    nb = fni(2,:)';
    
    % Ensure outward facing normals
%     if ici(1,:)*fni(1,:)'<0
%         na = -na;
%     end
%     
%     if ici(2,:)*fni(2,:)'<0
%         nb = -nb;
%     end
    
    % Calculate the edge normal vectors
    na12 = [0 -n12(3) n12(2); n12(3) 0 -n12(1); -n12(2) n12(1) 0]*na;
    nb21 = [0 -n21(3) n21(2); n21(3) 0 -n21(1); -n21(2) n21(1) 0]*nb;
    
    % Ensure outward-pointing edge normals
    a2b = ici(2,:)-ici(1,:);
    
    if a2b*na12<0
        na12 = -na12;
    end
    
    if -a2b*nb21<0
        nb21 = -nb21;
    end
    
    Ee = na*(na12') + nb*(nb21');
    
    % Vector from the field point to the edge
    re = V(E(i,1),:)'-r;
    
    % Edge length
    e = norm(V(E(i,2),:)-V(E(i,1),:));
    
    % Vectors from the field point to the edge endpoints
    R1 = V(E(i,1),:)'-r;
    R2 = V(E(i,2),:)'-r;
    
    r1 = norm(R1);
    r2 = norm(R2);
    
    % Calculate the logarithm expression, Le
    Le = log((r1+r2+e)/(r1+r2-e));
    
    % Sum the edge gravity contributions
    sumEdge = sumEdge + Ee*re*Le;
    
    % Calculate variational terms (if requested)
    if nargout > 1
        sumEdgeA = sumEdgeA + Ee*Le;
    end
    
end

%% Faces
% Sum the gravitational contributions of each face
for i=size(fn,1):-1:1
   
    Ff = (fn(i,:)')*(fn(i,:)')';
    
    % Get face vertices
    P = tr.X(tr(i,:),:)';
    
    % Calculate vectors from field point to face vertices
    R1 = P(:,1)-r;
    R2 = P(:,2)-r;
    R3 = P(:,3)-r;
    
    r1 = norm(R1);
    r2 = norm(R2);
    r3 = norm(R3);
    
    cR23 = [0 -R2(3) R2(2); R2(3) 0 -R2(1); -R2(2) R2(1) 0]*R3;
    
    % Calculate the solid angle term, wf
    wf = 2*atan2(R1'*cR23,...
        (r1*r2*r3+r1*R2'*R3+r2*R3'*R1+r3*R1'*R2));
    
    % Sum the face gravity contributions
    sumFace = sumFace + Ff*R1*wf;
    sumWf = sumWf + wf;
    
    % Calculate variational terms (if requested)
    if nargout > 1
        sumFaceA = sumFaceA + Ff*wf;
    end
    
end

%% Output Results

% Use the Laplacian to check if the field point is inside the polyhedron
if abs(sumWf-4*pi) < 1e-6
    warning('polydyn:InsidePoly','Field point is inside the polyhedron!')
end

% Acceleration at field point
a = G*rho*(-sumEdge + sumFace);

% State vector derivative
xdot = [v;D'*a];

% Variational terms
if nargout > 1
    A = [zeros(3,3) eye(3,3); 
         D'*G*rho*(sumEdgeA - sumFaceA)*D zeros(3,3)];
end

% Process noise
if nargout == 3
    Q = diag([0 0 0 1e-9 1e-9 1e-9].^2);
end

end