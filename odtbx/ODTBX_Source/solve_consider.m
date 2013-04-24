classdef solve_consider
    %solve_consider Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solve_param = [];
        solve_user_order = [];
        solve_func_order = [];
        dyn_cons_param = [];
        dyn_cons_user_order = [];
        dyn_cons_func_order = [];
        loc_cons_param = [];
        loc_cons_user_order = [];
        loc_cons_func_order = [];
    end
    
    methods
        function obj = solve_consider(solve, dyn_cons, loc_cons)
            if nargin > 0
                obj.solve_param = solve;
                obj.dyn_cons_param = dyn_cons;
                obj.loc_cons_param = loc_cons;
            else
                disp 'Need solve for and consider parameters'
            end
            
            if (length(solve_param) > 0)
                solve_user_order = 1:length(solve_param);
            end
            
            if (length(dyn_cons_param) > 0)
                dyn_cons_user_order = 1:length(dyn_cons_param);
            end
            
            if (length(loc_cons_param) > 0)
                loc_cons_user_order = 1:length(loc_cons_param);
            end
        end
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            
        end
    end
    
end

