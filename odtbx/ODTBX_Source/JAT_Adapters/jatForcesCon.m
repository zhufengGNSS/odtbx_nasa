function [xDot,A,Q] = jatForcesCon(obj,t,x,jatWorld)
% JATFORCESCON  Returns derivatives of JAT Force models (state in meters).
%
%   xDot = jatForcesCon(obj,t,x,jatWorld) returns the derivatives of the JAT
%   force models specified in the given jatWorld java object for an 
%   Earth-centric orbit. jatWorld is created by calling CREATEJATWORLD once
%   prior to simulation/integration.
%
%   JAT expects the input state to be in units of m and m/s. This file
%   expects the same input units and keeps the same units for the outputs.
%
%   [xDot,A] = jatForcesCon(obj,t,x,jatWorld) also returns the Jacobian (the state
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
%   [xDot,A,Q] = jatForcesCon(obj,t,x,jatWorld) also returns the process noise
%   matrix that is hardcoded within this file.
%
%   Consider parameters currently accepted: 'EARTH-GM', 'SOLAR-GM',
%   'LUNAR-GM', 'SOLRAD 8'
%
%
%   INPUTS
%   VARIABLE        SIZE    	DESCRIPTION (Optional/Default)
%       obj         object      Object containing solve for and consider
%                               parameters (solve_consider.m)
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
%   Phillip Anderson    6/7/2013        Incorporated solve for and consider
%                                       functionality from John Gaebler

    [nx,nt] = size(x);
    % Preallocate arrays
    dyn_max = length(obj.solve.user_order) + length(obj.dyn_cons.user_order);
    xDot = zeros(nx,nt);
    A = zeros(nx,dyn_max,nt);
    Q = zeros(nx,nx,nt);

    [Fei,Fsi,Fmi,Fsri]=deal([]);

    % Make the code more readable and reduce hash lookups
    if isKey(obj.param_order, 'EARTH-GM')
        Fei = obj.param_order('EARTH-GM');
    end
    if isKey(obj.param_order, 'SOLAR-GM')
        Fsi = obj.param_order('SOLAR-GM');
    end
    if isKey(obj.param_order, 'LUNAR-GM')
        Fmi = obj.param_order('LUNAR-GM');
    end
    if isKey(obj.param_order, 'SOLRAD 8')
        Fsri = obj.param_order('SOLRAD 8');
    end
    % To add parameters, introduce a new key here:
%     if isKey(obj.param_order, 'MYNEWKEY')
%         Fmynewkey = obj.param_order('MYNEWKEY');
%     end

    x=x(1:6,:)*1000; % conversion km -> m

    % Check the units
    for i=1:nt
        xDot(1:6,i) = jatWorld.derivs(t(i),x(:,i))/1000;  % conversion m -> km
    end

    SRPjat=nan(3,nt);
    acc_earth=nan(3,nt);
    acc_sun=nan(3,nt);
    acc_moon=nan(3,nt);

    if( nargout > 1 )
%                 A = [];
        useGravPartial = ( jatWorld.get_j2_gravity_flag() || jatWorld.get_2body_gravity_flag() );
        useDragPartial = jatWorld.get_drag_flag();

        for i=1:nt
            if useGravPartial
                A = zeros(nx,nx,nt);    
                A(1:3,4:6,i) = eye(3);
                if( useGravPartial )
                    A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.gravityPartials(x(1:3,i));
                end
            end

            % JAT Force List order:
            % 0 Earth Gravity
            % 1 Solar Gravity
            % 2 Lunar Gravity
            % 3 Drag
            % 4 Solar Radiation Force
            %
            % If a force is not added, all higher slots are shifted.
            % Put additional parameters after these
            cur = 0;
            
            % See ODEAS spec for more information as to these derivations
            if Fei % earth gravity error
                % need spacecraft position and earth GM
                acc_earth(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fei,i) = A(4:6,Fei,i) + acc_earth(:,i);
                cur=cur+1;
            end
            if Fsi % solar gravity error
                % need spacecraft inertial position, sun inertial position,
                % and sun GM
                acc_sun(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fsi,i) = A(4:6,Fsi,i) + acc_sun(:,i);
                cur=cur+1;
            end
            if Fmi % lunar gravity error
                % need spacecraft inertial position, lunar inertial
                % position, and lunar GM
                acc_moon(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fmi,i) = A(4:6,Fmi,i) + acc_moon(:,i);
                cur=cur+1;
            end
            if useDragPartial % Drag
                A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.dragPartials();
                cur=cur+1;
            end
            if Fsri % Solar radiation
                SRPjat(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray;
                A(4:6,Fsri,i) = A(4:6,Fsri,i) + (SRPjat(:,i))/1000;% m -> km          +[-3.05073589;6.75748831;2.03244886]*1e-9
                cur=cur+1;
            end

            % If there are any other consider parameters, their data
            % collection should go here.
            % Reference ODEAS spec to find dForce/dError
%             if Fmynewkey % New key as defined above
%                 Accel(:,i) = ...
%                     jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray;
%                 A(4:6,Fmynewkey,i) = A(4:6,Fmynewkey,i) + (Accel(:,i))/1000;% m -> km          +[-3.05073589;6.75748831;2.03244886]*1e-9
%                 cur=cur+1;
%             end
        end
    end

    if nargout == 3,
        Q = zeros(nx,nx,nt);
        Q(4:6,4:6,:) = repmat(diag([1e-9 1e-9 1e-9].^2),[1 1 nt])/1000^2;
        Q(Fei,Fei,end) = 0;
        Q(Fsi,Fsi,end) = 0;
        Q(Fmi,Fmi,end) = 0;
        Q(Fsri,Fsri,end) = 0;
    end
end

