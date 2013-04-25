classdef solve_consider
    %solve_consider Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solve = struct('param', [], 'user_order', [], 'func_order', []);
        dyn_cons = struct('param', [], 'user_order', [], 'func_order', []);
        loc_cons = struct('param', [], 'user_order', [], 'func_order', []);
        external_func = 'jat';
    end
    
    
    methods
        function obj = solve_consider(solve_in, dyn_cons_in, loc_cons_in, external_func_in)
            if nargin > 0
                % Provides the opportunity to define solve for/consider
                % params when the class is created.
                obj.solve.param = solve_in;
                obj.dyn_cons.param = dyn_cons_in;
                obj.loc_cons.param = loc_cons_in;
                obj.external_func = external_func_in;
            end
            
            if (~isempty(obj.solve.param))
                obj.solve.user_order = 1:length(obj.solve.param);
            end
            
            if (~isempty(obj.dyn_cons.param))
                obj.dyn_cons.user_order = 1:length(obj.dyn_cons.param);
            end
            
            if (~isempty(obj.loc_cons.param))
                obj.loc_cons.user_order = 1:length(obj.loc_cons.param);
            end
        end
        
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            % Interface maintains compatibility with jatForces?
            if (obj.external_func == 'jat')
                
            elseif (obj.external_func == 'gmat')
                
            else
                disp 'External function not recognized.'
            end
        end
        
    end
    
end

