classdef est_control < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        est_type
        myest
        
    end
    
    methods
        
        function obj = est_control(varargin)
            if nargin >= 1,
                obj.est_type = varargin{1};
                obj.myest = feval(obj.est_type, varargin{2:end});
            else
                disp "Must provide inputs!";
            end
        end
        
        
        function set_controllers(obj,varargin)
            if nargin >= 2,
                obj.myest.events_fcn = varargin{1};
                obj.myest.control_events_fcn = varargin{2};
            else
                obj.myest.events_fcn = @obj.events;
                obj.myest.control_events_fcn = @obj.control_events;
            end
        end
        
        
        function varargout = run_sim(obj)
            [t,Xhat,Phat,e,Y] = obj.myest.run_estimator();
            if nargout >= 3,
                varargout{1} = t;
                varargout{2} = Xhat;
                varargout{3} = Phat;
            end
            if nargout >= 4,
                varargout{4} = e;
            end
            if nargout >= 5,
                varargout{5} = Y;
            end
        end
        
        function [value,isterminal,direction] = events(obj,t,X,varargin)
            % See header in integev.m for details on event function formats
%             disp "Controller"
            % Consider using functions for conditions

            % Event 1:
            condition1 = X(1); %  We can change this to be anything related to t or X
            terminal1 = 0;
            direction1 = 0;
            
            % Event 2:
            condition2 = t - 295; %  We can change this to be anything related to t or X
            terminal2 = 1;
            direction2 = 0;
                        
            % Put all the values together to be returned
            value = [condition1; condition2]; % Condition to look for
            isterminal = [terminal1; terminal2]; % Do we need to end the integration at this point?
            direction = [direction1; direction2]; % Is there a direction involved?
        end % events
        
        
        function [X_state_mod, Phi_mod] = control_events(obj,t,X,Phi,varargin)
%             disp "Controller"
            X_state_mod = X;
            Phi_mod = Phi;
            
        end
        
    end
    
end

