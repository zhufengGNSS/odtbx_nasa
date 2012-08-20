function Sig_a = opt_sched(dynfun,datfun,tspan,Xo,options,dynarg,datarg,Pa)
% OPT_SCHED Calculates the sensitivities of a measurement schedule.
%
% Sig_a = opt_sched(dynfun,datfun,tspan,Xo,options,dynarg,datarg,Pa)
%
% OPT_SCHED provides the user a method to iteratively optimize a
% measurement schedule. It provides the sensitivities of the measurements
% at each time step with respect to the final states in a form to be
% displayed in SENSMOS. A higher sensitivity corresponds to a more
% important measurement. See Richard Battin's book "An Introduction to the
% Mathematics and Methods of Astrodynamics" pgs 690-693 for more detailed
% information.
%
% INPUTS
% VARIABLE     SIZE              DESCRIPTION
%    dynfun    function handle   Supplied dynamics function
%    datfun    function handle   Supplied measurement function
%    tspan     1xM               Time steps for analysis
%    Xo        Nx1               Initial 'N' states
%    options   structure         ODTBX options structure
%    dynarg    structure         Options for supplied dynamics function
%    datarg    structure         Options for supplied measurement function
%    Pa        NxN               Final covariance of analysis period
%
% OUTPUTS
%    Sig_a     NxHxM             Sensitivity matrix (N states, H
%                                  measurements, M time steps)
%
%
% keyword: optimize, schedule, measurement
% See also SENSMOS
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
%

[~,Xref,Phi] = integ(dynfun,tspan,Xo,options,dynarg);
Yref = feval(datfun,tspan,Xref,datarg);
[~,H,R] = ominusc(datfun,tspan,Xref,Yref,options,[],datarg);

I1 = eye(size(Xo,1));
I2 = eye(size(R,1));
M = I1;
L = zeros(size(Pa,1),size(Pa,2),size(Phi,3));
L(:,:,end) = -Pa*M*Pa;
for ii = size(Phi,3)-1:-1:1
    Phi_i = Phi(:,:,ii+1)\I1;
    L(:,:,ii) = Phi_i*L(:,:,end)*Phi_i';
end
for jj = size(L,3):-1:1
    Sig_a(:,:,jj) = (L(:,:,jj)+L(:,:,jj)')*(H(:,:,jj)'*(R(:,:,jj)\I2));
end
end