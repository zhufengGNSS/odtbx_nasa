classdef est_control < handle
    %EST_CONTROL Handles the events and controls for a user-specified specified estimator.
    %   Acts as a wrapper for a user-defined estimator. Users will define
    %   all of the events and controls in this class.
    
    properties
        % Estimator
        est_type
        myest
    end
    
    methods
        %% Estimator wrapper functions
        function obj = est_control(varargin)
            % Input Parsing and Setup            
            if nargin >= 1,
                obj.est_type = varargin{1};
                obj.myest = feval(obj.est_type, varargin{2:end});
            else
                disp "Must provide inputs!";
            end
        end
        
        
        function set_controllers(obj,varargin)
            % This function needs to be run to substitute the
            % events/controls functions from this class (or other specified functions)
            % for the defaults located in the estimator.
            if nargin >= 2,
                obj.myest.events_fcn = varargin{1};
                obj.myest.control_events_fcn = varargin{2};
            else
                obj.myest.events_fcn = @obj.events;
                obj.myest.control_events_fcn = @obj.control_events;
            end
        end
        
        
        function varargout = run_sim(obj)
            % Run estimator
            [t,Xhat,Phat,e,Y] = obj.myest.run_estimator();
            
            % Output results
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
        
        
        %% Event conditions
        function [value,isterminal,direction] = events(obj,t,X,varargin)
            % This function is used to kick the integrator out of its loop
            % at a certain point (as determined by time or state
            % conditions) in order to perform an action desginated by
            % control_events.
            
            % See header in integev.m for details on event function formats
            
            % Event 1: Position
            [condition1,terminal1,direction1] = obj.event1(t,X);            
            
            % Event 2: Time
            [condition2,terminal2,direction2] = obj.event2(t,X);
                        
            % Put all the values together to be returned
            value = [condition1; condition2]; % Conditions to look for
            isterminal = [terminal1; terminal2]; % Do we need to end the integration at this point?
            direction = [direction1; direction2]; % Is there a direction involved?
        end % events
        
        
        function [value,isterminal,direction] = event1(obj,t,X)
            % This function contains a condition that could halt the
            % propagation of the system being tested. No action is taken to
            % implement any controls in this function. 
            % This function must be called from the events function (or
            % nothing will happen).

            % We can change this to be anything related to t or X
            value = X(1); % Condition that will trigger the event when equal to zero
            isterminal = 1; % Whether the condition will halt propagation
            direction = 0; % If the event is purely zero-finding or if it is directional
            
        end % event2
        
        
        function [value,isterminal,direction] = event2(obj,t,X)
            % This function contains a condition that could halt the
            % propagation of the system being tested. No action is taken to
            % implement any controls in this function. 
            % This function must be called from the events function (or
            % nothing will happen).

            % We can change this to be anything related to t or X
            value = t - 290; % Condition that will trigger the event when equal to zero
            isterminal = 1; % Whether the condition will halt propagation
            direction = 0; % If the event is purely zero-finding or if it is directional
            
        end % event2
        
        
        %% Controls
        function [X_state_mod, P_mod] = control_events(obj,t,X,P,varargin)
            % This function is used to change the state/covariance once a
            % condition has been detected.
            
            % This function currently only operates on the most recent
            % event that halts the integration (isterminal = 1). At this
            % time, non-terminal events are ignored.
            
            % Event1 control
            [value,~,~] = obj.event1(t,X);
            if (value == 0)
                
            end
            
            % Event2 control
            [value,~,~] = obj.event2(t,X);
            if (value == 0)
                 X(11) = X(11).*1.5;
            end
            
            X_state_mod = X;
            P_mod = P;
            
        end % control_events
        
    end % methods
    
end % est_control

