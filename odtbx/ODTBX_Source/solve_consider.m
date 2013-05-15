classdef solve_consider
    %solve_consider Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solve = struct('param', {''}, 'user_order', {''});
        dyn_cons = struct('param', {''}, 'user_order', {''});
        loc_cons = struct('param', {''}, 'user_order', {''});
        param_order = containers.Map;
        external_func = 'jat';
    end
    
    
    methods
        function obj = solve_consider(varargin)
            % Provides the opportunity to define solve for/consider
            % params when the class is created.
            if nargin >= 1 % First variable is always the solve param
                obj.solve.param = varargin{1};
            end
            if nargin == 2 % If length is two, short state, second is external function
                obj.external_func = varargin(2);
            end
            if nargin > 2 % If length is bigger than two, we'll get a full state
                obj.dyn_cons.param = varargin{2}; 
            end
            if nargin >= 3 % Full state
                obj.loc_cons.param = varargin{3};
            end
            if nargin >= 4 % Full state
                obj.external_func = varargin(4);
            end
            
            map_params(obj);
%             obj.order_size = max(cell2mat(obj.param_order.values()));
        end
        
        
        function map_params(obj)
            % Maps all of the consider functions to values representing
            % the order they should be listed in the A matrix.
            if (~isempty(obj.solve.param))
                for order = 1:length(obj.solve.param)
                    % Assign the value
                    obj.solve.user_order{order,1} = order;
                    % Add to map
                    obj.param_order(obj.solve.param{order}) = ...
                        obj.solve.user_order{order,1};
                end
            end  
            current_len = length(obj.solve.param);
            if (~isempty(obj.dyn_cons.param))
                for order = 1:length(obj.dyn_cons.param)
                    % Assign the value
                    obj.dyn_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    obj.param_order(obj.dyn_cons.param{order}) = ...
                        obj.dyn_cons.user_order{order,1};
                end
            end
            current_len = current_len + length(obj.dyn_cons.param);
            if (~isempty(obj.loc_cons.param))
                for order = 1:length(obj.loc_cons.param)
                    % Assign the value
                    obj.loc_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    obj.param_order(obj.loc_cons.param{order}) = ...
                        obj.loc_cons.user_order{order,1};
                end
            end
            obj.solve.user_order
            obj.dyn_cons.user_order
            obj.loc_cons.user_order
        end
        
        
        function unmap_params(obj)
            remove(obj.param_order, obj.param_order.keys());
        end
        
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            % containers.Map objects are treated as handles, so the hash is
            % shared between different instances of the solve_consider
            % class. We need to make sure we're using the correct hash when
            % we're running the external force models.
            
            unmap_params(obj);
            map_params(obj);
            
            % Interface maintains compatibility with jatForces?
            if (strcmpi(obj.external_func, 'jat'))
                [xDot,A,Q] = obj.jatForces(t,x,jatWorld);
            elseif (strcmpi(obj.external_func, 'gmat'))
                [xDot,A,Q] = obj.gmatForces(t,x,jatWorld);
            else
                disp 'External function not recognized.'
            end
        end
        
        
        function [xDot,A,Q] = jatForces(obj,t,x,jatWorld)
%             x
            % John Gaebler's Code, revised
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

            x=x(1:6,:)*1000; % conversion km -> m

            for i=1:nt
                xDot(1:6,i) = jatWorld.derivs(t(i),x(:,i))/1000;  % conversion m -> km
            end

            SRPjat=nan(3,nt);
            acc_earth=nan(3,nt);
            acc_sun=nan(3,nt);
            acc_moon=nan(3,nt);

            if( nargout > 1 )
                A = [];
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
                    cur = 0;
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
                    if useDragPartial 
                        A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.dragPartials();
                        cur=cur+1;
                    end
                    if Fsri
                        SRPjat(:,i) = ...
                            jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray;
                        A(4:6,Fsri,i) = A(4:6,Fsri,i) + (SRPjat(:,i))/1000;% m -> km          +[-3.05073589;6.75748831;2.03244886]*1e-9
                        cur=cur+1;
                    end
                    
                    % If there are any more consider parameters, they
                    % should go here
                    
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
        
        
        function [xDot,A,Q] = gmatForces(obj,t,x,jatWorld)
            % Need GMAT interface info
        end
        
    end
    
end

