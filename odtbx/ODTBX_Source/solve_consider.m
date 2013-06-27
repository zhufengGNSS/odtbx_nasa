classdef solve_consider
    %solve_consider Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solve = struct('param', {''}, 'user_order', {''});
        dyn_cons = struct('param', {''}, 'user_order', {''});
        loc_cons = struct('param', {''}, 'user_order', {''});
        param_order = containers.Map;
        external_force = 'jat';
        external_meas = 'gsmeas';
    end
    
    
    methods
        %% Constructor
        function obj = solve_consider(varargin)
            % Provides the opportunity to define solve for/consider
            % params when the class is created.
            if nargin >= 1 % First variable is always the solve param
                obj.solve.param = varargin{1};
            end
            if nargin == 2 % If length is two, short state, second is external function
                obj.external_force = varargin{2};
            end
            if nargin > 2 % If length is bigger than two, we'll get a full state
                obj.dyn_cons.param = varargin{2}; 
            end
            if nargin >= 3 % Full state
                obj.loc_cons.param = varargin{3};
            end
            if nargin >= 4 % Full state
                obj.external_force = varargin{4};
            end
            if nargin >= 5 % Full state
                obj.external_meas = varargin{5};
            end
            
            % containers.Map objects are treated as handles class, so every
            % instance of the solve_consider class point back to the same
            % Map. 
            % Because there is no easy way to do deep copy (clone) of a
            % Map, we have to use a little trick inspired by:
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/278249
            % It's possible to force Matlab to recognize a containers.Map
            % object as value-based if it's returned from a function.
            [obj, obj.param_order] = map_params(obj);
            
%             length(obj.solve.param)
%             obj.solve.user_order
%             length(obj.dyn_cons.param)
%             obj.dyn_cons.user_order
%             length(obj.loc_cons.param)
%             obj.loc_cons.user_order
        end
        
        
        %% Parameter Mapping
        function [obj,temp_map] = map_params(obj)
            temp_map = containers.Map();
            % Maps all of the consider functions to values representing
            % the order they should be listed in the A matrix.
            if (~isempty(obj.solve.param))
                for order = 1:length(obj.solve.param)
                    % Assign the value
                    obj.solve.user_order{order,1} = order;
                    % Add to map
                    temp_map(obj.solve.param{order}) = ...
                        obj.solve.user_order{order,1};
                end
            end  
            current_len = length(obj.solve.param);
            if (~isempty(obj.dyn_cons.param))
                for order = 1:length(obj.dyn_cons.param)
                    % Assign the value
                    obj.dyn_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    temp_map(obj.dyn_cons.param{order}) = ...
                        obj.dyn_cons.user_order{order,1};
                end
            end
            current_len = current_len + length(obj.dyn_cons.param);
            if (~isempty(obj.loc_cons.param))
                for order = 1:length(obj.loc_cons.param)
                    % Assign the value
                    obj.loc_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    temp_map(obj.loc_cons.param{order}) = ...
                        obj.loc_cons.user_order{order,1};
                end
            end

            % temp_map is returned from this function to overcome
            % shortcomings of containers.Map objects.

        end
        
        
        %% External Force Models
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            
            % Interface maintains compatibility with jatForces?
            if (strcmpi(obj.external_force, 'jat'))
                [xDot,A,Q] = jatForcesCon(obj,t,x,jatWorld);
            elseif (strcmpi(obj.external_force, 'gmat'))
                [xDot,A,Q] = gmatForces(obj,t,x,jatWorld);
            else
                disp 'Force function not recognized.'
            end
        end
        
        
        function [xDot,A,Q] = gmatForces(obj,t,x,jatWorld)
            % Need GMAT interface info
        end
        
        
        %% Measurement Models
        
        function [y,H,R] = meas(obj,t,x,options)
            cons_iono_flag = any(strncmpi(obj.loc_cons.param,'ION',3));
            cons_trop_flag = any(strncmpi(obj.loc_cons.param,'TRP',3));
            
            if cons_iono_flag == true
                options=setOdtbxOptions(options,'useIonosphere',true);
            end
            if cons_trop_flag == true
                options=setOdtbxOptions(options,'useTroposphere',true);
            end
            
            
            if (strcmpi(obj.external_meas, 'gsmeas'))
                [y,H,R] = gsmeasCon(obj,t,x,options);
            else
                disp 'Measurement function not recognized.'
            end
        end
    
    end
    
end

