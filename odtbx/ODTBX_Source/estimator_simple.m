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
%         X_state
        t
        X
        Xhat
        Phat
        e
        y
        Y
        Pa
        Pv
        Pw
        Phata
        Phatv
        Phatw
        Phatm % new estseq
        Phatt % new estseq
        Sigma_a % estbat
        Sig_sa % estseq
        eflag % estseq
        P
        Pdy
        Pdyt
        Pm % estseq
        restartRecord % estseq
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
            value = [0; 0]; % Condition to look for
            isterminal = [0; 0]; % Do we need to end the integration at this point?
            direction = [0; 0]; % Is there a direction involved?
        end % events
        
        
        function [X_state_mod, Phi_mod] = control_events_default(obj,t,X,Phi,varargin)
            % This function is used to change the state/covariance once a
            % condition has been detected.
            X_state_mod = X;
            Phi_mod = Phi;
            
        end
    end % Methods
    
    methods(Static)
        
        function y = refine(u,refine)
            y = [reshape([u(1:end-1);repmat(u(1:end-1),refine,1)+...
            cumsum(repmat(diff(u)/(refine+1),refine,1),1)],[],1);u(end)]';
        end
        
    end % Methods (Static)
end % Class

