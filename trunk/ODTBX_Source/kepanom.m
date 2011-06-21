function angle = kepanom(S, outstring)
%  KEPANOM  Solves for one angle of anomaly given another using Kepler's Equation
%
%   INPUTS
%   VARIABLE        DESCRIPTION
%      S            Data structure containing the following fields.
%                   Only one angle of anomaly is needed as input.
%           S.ecc   eccentricty (unitless) -required for all conversions
%           S.E     eccentric anomaly (radians)
%           S.M     mean anomaly (radians)
%           S.tran  true anomaly (radians) 
%                   (True anomaly can also be specified in the field S.nu.)
%
%      outstring    String defining the desired angle of anomaly
%                    = 'E';    % eccentric anomaly
%                    = 'M';    % mean anomaly
%                    = 'tran'; % true anomaly
%                    = 'nu';   % true anomaly
%
%   OUTPUTS
%      angle        Angle of the deired anomaly (eccentric, mean, or true)
%                    in radians
%
%   EXAMPLES
%      S.ecc = 0.1;
%      S.M   = 45*pi/180;
%      E     = kepanom(S, 'E');
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Steven Hughes       04/09/2001      Created original subfunctions
%   Kevin Berry         12/15/2008      Created kepanom.m based on provided
%                                       subfunctions

switch upper(outstring)
   case 'E'
      if any(strcmpi('tran',fieldnames(S)) | strcmpi('nu',fieldnames(S)))
          angle = nu2E(S);
      elseif any(strcmpi('M',fieldnames(S)))
          angle = M2E(S);
      end
   case 'M'
      if any(strcmpi('E',fieldnames(S)))
          angle = E2M(S);
      elseif any(strcmpi('tran',fieldnames(S)) | strcmpi('nu',fieldnames(S)))
          angle = nu2M(S);
      end
   case {'TRAN' , 'NU'}
      if any(strcmpi('E',fieldnames(S)))
          angle = E2nu(S);
      elseif any(strcmpi('M',fieldnames(S)))
          angle = M2nu(S);
      end      
   otherwise
      error('Unknown angle requested.')
end

end

%% Subfunctions
function E = nu2E(S)
% Solves for E given nu and ecc using Kepler's Equations
if isfield(S,'tran')
    ta = S.tran;
elseif isfield(S,'nu')
    ta = S.nu;
end
cosnu = cos(ta);
sin_E = sqrt(1 - S.ecc*S.ecc)*sin(ta)/(1+S.ecc*cosnu);       %Vallado pg. 213,  Eq. 4-9
cos_E = (S.ecc+cosnu)/(1+S.ecc*cosnu);
E = mod(atan2(sin_E,cos_E),2*pi);
end

function E = M2E(S)
% Solves Keplers Problem for E given M and ecc using Newton Iteration
if ( -pi < S.M < 0 ) %  pick initial guess for E
   E = S.M - S.ecc;
else
   E = S.M + S.ecc;
end
tol = 1e-8;%  Iterate until tolerance is met
diff = 1;
while ( diff > tol )
   f = S.M - E + S.ecc*sin(E);
   fprime = 1 - S.ecc*cos(E);
   E_new = E + f/fprime;
   diff = abs(E_new - E);
   E = E_new;
end
end

function M = E2M(S)
% Solves for M given E and ecc using Kepler's Equations
M = S.E - S.ecc* sin(S.E);
end

function M = nu2M(S)
% Solves for M given nu and ecc using Kepler's Equations
S.E = nu2E(S);
M = E2M(S);
end

function nu = E2nu(S)
% Solves for nu given E and ecc using Kepler's Equations
%  -uses atan2 to insure the proper quadrant
cos_E = cos(S.E);
sin_nu = sqrt(1 - S.ecc*S.ecc)*sin(S.E)/(1-S.ecc*cos_E);       %Vallado pg. 214,  Eq. 4-10
cos_nu = (cos_E - S.ecc)/(1 - S.ecc*cos_E);              %Vallado pg. 214,  Eq. 4-12
nu = mod(atan2(sin_nu,cos_nu),2*pi);
nu = mod(nu,2*pi);
end  

function nu = M2nu(S)
% Solves for nu given M and ecc using Kepler's Equations
S.E = M2E(S);
nu = E2nu(S);
end