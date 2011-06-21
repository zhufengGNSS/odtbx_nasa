function varargout = RK4(ode,tspan,y0,options,varargin)
%
% This RK4 function was created for ODEAS validation.
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

X=y0;
delta_time=varargin{1}.OdeSolvOpts.MaxStep;

T(1) = tspan(1);
Y(1,:) = y0';

lastind=ceil(tspan(end)/delta_time)+1;
for i=2:lastind

    DX = feval(ode,tspan(1),X,varargin{:});
    K1 = delta_time * DX;
    XT = X + K1/2;

    DX = feval(ode,tspan(1),XT,varargin{:});
    K2 = delta_time * DX;
    XT = X + K2/2;

    DX = feval(ode,tspan(1),XT,varargin{:});
    K3 = delta_time  * DX;
    XT = X + K3;

    DX = feval(ode,tspan(1),XT,varargin{:});
    K4 = delta_time * DX;

    X = X + (K1 + 2*K2 + 2*K3 + K4)/6;

    T(i) = T(i-1)+delta_time;
    Y(i,:) = X';
end

TOUT = tspan';
YOUT = interp1(T,Y,tspan);

varargout = {TOUT, YOUT};

