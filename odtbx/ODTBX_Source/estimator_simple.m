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
        function obj = estimator()
            dynfun = [];
            dynarg = [];
            datfun = [];
            datarg = [];
            tspan = [];
            Xo = [];
            Xbaro = [];
            Po = [];
            Pbaro = [];
            options = [];

        end
        
        function varargout = run_estimator(obj)
            % This is a dummy function. Subclasses of estimator should
            % overwrite this function with the function that will run the
            % estimator.
        end
    end % Methods
    
    methods(Static)
        
        function y = refine(u,refine)
            y = [reshape([u(1:end-1);repmat(u(1:end-1),refine,1)+...
            cumsum(repmat(diff(u)/(refine+1),refine,1),1)],[],1);u(end)]';
        end
        
    end % Methods (Static)
end % Class

