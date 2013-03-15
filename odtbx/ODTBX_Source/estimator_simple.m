classdef estimator_simple < handle
    % ESTIMATOR Estimator base class
    %   This class is used to define the properties common to all ODTBX
    %   estimators. It will promote uniformity in implementations and
    %   coding.
    
    properties
        % In
        dynfun
        dynarg
        datfun
        datarg
        tspan
        Xo
        Xbaro
        Po
        Pbaro
        options
        events_fcn
        control_events_fcn
        
        % Out
        t
        X
        Xhat
        Phat
        e
        Y
    end
    
    methods
        function obj = estimator(varargin)
            obj.dynfun = [];
            obj.dynarg = [];
            obj.datfun = [];
            obj.datarg = [];
            obj.tspan = [];
            obj.Xo = [];
            obj.Xbaro = [];
            obj.Po = [];
            obj.Pbaro = [];
            obj.options = [];
            
            obj.events_fcn = @obj.events_default;
            obj.control_events_fcn = @obj.control_events_default;

        end
        
        
        function varargout = run_estimator(obj)
            % This is a dummy function. Subclasses of estimator should
            % overwrite this function with the function that will run the
            % estimator.
        end
        
        
        function [value,isterminal,direction] = events_default(obj,t,X,varargin)
            % This function is used to kick the integrator out of its loop
            % at a certain point (as determined by time or state
            % conditions) in order to perform an action desginated by
            % control_events.
            
            % See header in integev.m for details on event function formats
            
            % Put all the values together to be returned
            value = [0]; % Condition to look for
            isterminal = [0]; % Do we need to end the integration at this point?
            direction = [0]; % Is there a direction involved?
        end % events
        
        
        function [X_state_mod, P_mod] = control_events_default(obj,t,X,P,varargin)
            % This function is used to change the state/covariance once a
            % condition has been detected.
            X_state_mod = X;
            P_mod = P;
            
        end
    end % Methods
    
end % Class

