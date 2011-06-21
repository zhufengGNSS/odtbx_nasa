function [xDot,A,Q] = jatForces_km(t,x,jatWorld)

% JATFORCES_KM  Returns derivatives of JAT Force models (state in km).
%
%   xDot = jatForces_km(t,x,jatWorld) returns the derivatives of the JAT
%   force models specified in the given jatWorld java object for an 
%   Earth-centric orbit. jatWorld is created by calling CREATEJATWORLD once
%   prior to simulation/integration.
%
%   This function is a wrapper function for jatForces.  All inputs and 
%   outputs are in kilometers and seconds while jatForces expects inputs in 
%   meters and seconds like JAT.
%
%   [xDot,A] = jatForces(t,x,jatWorld) also returns the Jacobian (the state
%   derivatives of the equations of motion) if possible; otherwise an []
%   empty matrix is returned. Partials are only returned from JAT for
%   atmospheric drag and gravity. The gravitational partials are either
%   the 2-body partials (if the input gravity model is '2body') or the J2 
%   partials. Gravitational partials associated with higher order effects 
%   are not included.  If the model used in computing xDot is of higher
%   order, than A can potentially have significant errors.  In this case, 
%   better filter performance can be obtained by forcing numerical 
%   computation of A by setting A=[], an empty matrix, either in this 
%   function or in a user-provided wrapper function that calls this 
%   function.
%
%   [xDot,A,Q] = jatForces_km(t,x,jatWorld) also returns the process noise
%   matrix that is hardcoded within this file.
%
%   INPUTS
%   VARIABLE        SIZE    	DESCRIPTION (Optional/Default)
%       t           (1x1)       Time since start of simulation (secs)
%       x           (6x1)       Input state (x,y,z position in km 
%                               followed by x,y,z velocity in km/s)
%       jatWorld    (1x1)       java object created by createJATWorld
%
%   OUTPUTS
%       xDot        (6x1)       Derivatives of state (km, seconds)
%       A           (6x6)       State transition matrix (km, seconds)
%       Q           (6x6)       Process noise (km^2, seconds)
%
%   keyword: JAT Forces
%   See also JATFORCES, SETODTBXOPTIONS, GETODTBXOPTIONS, CREATEJATWORLD
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Allen Brown         02/09/2009      Created
%   Sun Hur-Diaz        04/16/2009      Changed comment on xDot dimension
%                                       from Nx1 to 6x1
%   Allen Brown         11/23/2010      Updated comments.

% Call jatForces with a state in m, get answers back in m.
% No units conversion required on the state transition matrix itself.
if( nargout == 1 )
    xDotm = jatForces(t,x*1000,jatWorld);
elseif( nargout == 2 )
    [xDotm,A] = jatForces(t,x*1000,jatWorld);
else
    [xDotm,A,Qm] = jatForces(t,x*1000,jatWorld);
end

% Convert the returned state derivatives wrt time into km:
xDot = xDotm/1000;

% Manipulation of the state transition matrix is not required for a change
% in units.

% Convert the process noise into km
if exist('Qm','var')
    Q = Qm/(1000*1000);
end
