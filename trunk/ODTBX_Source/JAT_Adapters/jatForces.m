function [xDot,A,Q] = jatForces(t,x,jatWorld)

% JATFORCES  Returns derivatives of JAT Force models (state in meters).
%
%   xDot = jatForces(t,x,jatWorld) returns the derivatives of the JAT
%   force models specified in the given jatWorld java object for an 
%   Earth-centric orbit. jatWorld is created by calling CREATEJATWORLD once
%   prior to simulation/integration.
%
%   JAT expects the input state to be in units of m and m/s. This file
%   expects the same input units and keeps the same units for the outputs.
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
%   [xDot,A,Q] = jatForces(t,x,jatWorld) also returns the process noise
%   matrix that is hardcoded within this file.
%
%   INPUTS
%   VARIABLE        SIZE    	DESCRIPTION (Optional/Default)
%       t           (1xn)       Time since start of simulation (secs)
%       x           (6xn)       Input state (x,y,z position in meters 
%                               followed by x,y,z velocity in m/s)
%       jatWorld    (1x1)       java object created by createJATWorld
%
%   OUTPUTS
%       xDot        (6xn)       Derivatives of state (meters, seconds)
%       A           (6x6xn)     State transition matrix (meters, seconds)
%       Q           (6x6xn)     Process noise (meters, seconds)
%
%   keyword: JAT Forces
%   See also SETODTBXOPTIONS, GETODTBXOPTIONS, CREATEJATWORLD, JATFORCES_KM
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
%   Derek Surka         06/26/2007   	Original
%   (missing some revision entries)
%   Allen Brown         01/30/2009      Changed to meters & m/s only
%   Sun Hur-Diaz        08/03/2009      Add vectorized capability
%   Allen Brown         11/23/2010      Updated comments.

[nx,nt]=size(x);
xDot=zeros(nx,nt);
for i=1:nt
    xDot(:,i) = jatWorld.derivs(t(i),x(:,i));
end

if( nargout > 1 )
    A = [];
    useGravPartial = ( jatWorld.get_j2_gravity_flag() || jatWorld.get_2body_gravity_flag() );
    useDragPartial = jatWorld.get_drag_flag();
    if( useGravPartial || useDragPartial )
        A = zeros(6,6,nt);
        for i=1:nt
            A(1:3,4:6,i) = eye(3);
            if( useGravPartial )
                A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.gravityPartials(x(1:3,i));
            end
            if( useDragPartial )
                A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.dragPartials();
            end
        end
    end
end

if nargout == 3,
    Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2),[1 1 nt]);
end

end
