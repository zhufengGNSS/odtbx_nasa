function varargout = errellip(P, varargin)
% ERRELLIP  2-D and 3-D Error Ellipse plotter.
%   ERRELLIP(P) where P is a 2x2 or 3x3 covariance matrix, plots the
%   "1-sigma" 2-D error ellipse or 3-D error ellipsoid, respectively.  This
%   is the boundary defined by all 2-D or 3-D error vectors E that satisify
%   E'*inv(P)*E = K^2, for K = 1.  For the 3-D case, ERRELLIP also plots
%   2-D ellipses that are the projections of the ellipsoid onto the three
%   planes defining the coordinate system of P.
%
%   ERRELLIP(P,K) plots the "k-sigma" ellipsoid.
%
%   ERRELLIP(P,K,DTH) uses the grid spacing DTH (in degrees) for the
%   boundary rather than the default, which is based on the eccentricities
%   of the principal planes of the ellipsoid.  
%
%   H = ERRELLIP(P,...) returns a handle to the plot.
%
%   [X,Y] = ERRELLIP(P,...) and [X12,Y12,X13,Z13,Y23,Z23] = ERRELLIP(P,...)
%   return the data that define the 2-D plots for the 2-D and 3-D cases,
%   respectively.  To get the surface data for the 3-D case, use the handle
%   option above.
%
%keyword: Utilities, 
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

% Russell Carpenter
% NASA Johnson Space Center
% Original version (ellipse.m and fellipse.m) 04/18/91

% Russell Carpenter
% NASA Goddard Space Flight Center
% Date: 2005-02-22 12:07:47 -0500 (Tue, 22 Feb 2005) 

% Error Checking
[n,m] = size(P);
if n~=m,
    error('ERRELLIP:notSquare', 'Input is not square.')
end
if n>3|n<2,
    error('ERRELLIP:dimWrong', 'Input must be 2x2 or 3x3.')
end
[V,D] = eig(P);
lambda = diag(D);
if any(lambda <= 0),
    error('ERRELLIP:notPosDef', 'Input is not positive definite.')
end

% Parse or Define "Sigma-index" of Ellipse
if nargin > 1,
    k = varargin{1};
    if isempty(k),
        k = 1;
    end
else,
    k = 1;
end

% Parse or Define Grid Spacing
if nargin == 3,
    dth = varargin{2}*pi/180;
    if n == 3,
        dth = dth*ones(1,3);
    end
else,
    if n == 3,
        e23 = ecc(lambda(2),lambda(3));
        e13 = ecc(lambda(1),lambda(3));
        dth(3) = gspace(e23);
        dth(2) = gspace(e13);
    end
    e12 = ecc(lambda(1),lambda(2));
    dth(1) = gspace(e12);
end

if n == 2, % Ellipse for n = 2 based on polar coordinates (fellipse.m code)
    a = lambda(1);
    b = lambda(2);
    th = 0:dth:2*pi;
    r = sqrt(a*b*ones(1,length(th))./(a*sin(th).^2 + b*cos(th).^2));
    u = k*r.*cos(th);
    v = k*r.*sin(th);
    tmp = V*[u; v];
    x = tmp(1,:);
    y = tmp(2,:);
    if nargout < 2,
        h = plot(x,y);
        axis equal
        if nargout == 1,
            varargout{1} = h;
        end
    else,
        varargout{1} = x;
        varargout{2} = y;
    end
    
else, % Ellipsoid for n = 3
    % First compute the 2-D projections:
    [x12,y12] = errellip(P(1:2,1:2),k);
    [x13,z13] = errellip(P([1 3],[1 3]),k);
    [y23,z23] = errellip(P(2:3,2:3),k);
    % Only do the (time-consuming) surface computation if asked:
    if nargout < 3,
        for j = 3:-1:1,
            lim = k*sqrt(P(j,j))*[-1,1];
            rng{j} = linspace(lim(1), lim(2), ceil(pi/dth(j)));
        end
        % Create set of data points throughout the space and evaluate the
        % pdf at each of those points:
        [u,v,w] = meshgrid(rng{1},rng{2},rng{3});
        f = pdf3(u,v,w,P);
        % Now evaluate the pdf on one of the "k-sigma" eigenvectors which
        % will define where in the space to draw the surface:
        kvec = V(:,1)*k*sqrt(lambda(1));
        fk = pdf3(kvec(1),kvec(2),kvec(3),P);
        % Create the "k-sigma" isosurface as a patch object; the isonormals
        % function is a refinement to produce a smoother-looking surface:
        h = patch(isosurface(u,v,w,f,fk));
        isonormals(u,v,w,f,h)
        daspect([1,1,1])
        view(126,24)
        % User can change these using the output handle:
        caxis([0 0.5])
        cmap = colormap;
        set(h,'FaceColor',cmap(end,:),'EdgeColor','none','facealpha',0.75)
        camlight headlight, camlight left
        lighting phong
        if nargout == 1,
            varargout{1} = h;
        end
        % Add 2-D projections:
        hold on
        plot3(x12,y12,zeros(length(x12)))
        plot3(x13,zeros(length(x13)),z13)
        plot3(zeros(length(y23)),y23,z23)
        hold off
    else,
        varargout{1} = x12;
        varargout{2} = y12;
        varargout{3} = x13;
        varargout{4} = z13;
        varargout{5} = y23;
        varargout{6} = z23;
    end
end

function e = ecc(a,b)
% Eccentricity.
if a>b,
    e = sqrt((a-b)/a);
elseif b>a,
    e = sqrt((b-a)/b);
else,
    e = 0;
end

function dth = gspace(e)
d2r = pi/180;
% Grid spacing based on eccentricity
if e<0.9,
    dth = 1e01*d2r;
elseif e<0.99,
    dth = 1e00*d2r;
else,
    dth = 1e-1*d2r;
end

%function f = pdf3(u,v,w,lambda)
% 3-D Gaussian PDF
%f = 1/sqrt(2*pi)^3/sqrt(prod(lambda))...
%    * exp(-(u.^2/lambda(1) + v.^2/lambda(2) + w.^2/lambda(3))/2);

function f = pdf3(u,v,w,P)
% 3-D Gaussian PDF
M = inv(P);
f = 1/sqrt(8*pi^3*det(P))*exp(-(M(1,1)*u.^2 + M(2,2)*v.^2 + M(3,3)*w.^2 ...
    + 2*M(1,2)*u.*v + 2*M(1,3)*u.*w + 2*M(2,3)*v.*w)/2);