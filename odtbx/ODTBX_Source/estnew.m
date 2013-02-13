classdef estnew < estimator_simple
    % ESTNEW A new, simple ODTBX estimator
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = estnew()
%             Initialize inputs (state, covariance, timespan, models)
        end
    
        function obj = run_estimator()

            for time_sim = 1:length(obj.tspan)
                
            end
%             For i=1:length(tspan)
% 
%             Calculate true/estimated measurement at tspan(i)
%             Perform measurement update (the 'kalmup' function will be useful)
% 
%             While ~done
%             Propagate true/estimated state [from tspan(i) to tspan(i+1)] together
%             using dynfun
%                     - use ode event function if supplied by user
%             Check integration time to see if t_end = tspan(i+1)
%                     - if yes, done = true and leave loop
%                     - if not, terminal event was detected, adjust state/covariance based on
%             user-supplied function and try to propagate to tspan(i+1) again.  Keep
%             going until t_end = tspan(i+1)
%             End while loop
% 
%             End for loop
% 
%             Calculate errors, package data, and output results (state errors,
%             covariance, residuals, etc)

        end
        
    end
    
end

