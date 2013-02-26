classdef estnew < estimator_simple
    % ESTNEW A new, simple ODTBX estimator
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = estnew(varargin)
%             Initialize inputs (state, covariance, timespan, models)
            %% Input Parsing and Setup
            % Parse the input list and options structure.  Pre-allocate arrays, using a
            % cell index for the monte carlo cases, which will avoid the need for each
            % case to have time series at common sample times.  Use an extra dimension
            % "on the right" within each monte carlo case to accomodate the time
            % series, which will avoid the need for conversions from cell to double for
            % plotting.  Where it makes sense, use cell indices to partition
            % large matrices into submatrices, to avoid the need for opaque indexing
            % computations.
            %
            % This should be a subfunction, or if there is a lot of commonality with
            % estseq's version, a private function.
            %
            % If there are no input arguments, perform a built-in self-test.  If there
            % are no output arguments, then plot the results of the input self-test
            % number as a demo.

            if nargin >= 4,
                if all(isfield(varargin{1}, {'tru','est'})),
                    obj.dynfun = varargin{1};
                else
                    obj.dynfun.tru = varargin{1};
                    obj.dynfun.est = varargin{1};
                end
                if all(isfield(varargin{2}, {'tru','est'})),
                    obj.datfun = varargin{2};
                else
                    obj.datfun.tru = varargin{2};
                    obj.datfun.est = varargin{2};
                end
                
                obj.tspan = varargin{3};
                
                if isstruct(varargin{4}),
                    obj.Xo = varargin{4}.Xo;
                    obj.Xbaro = varargin{4}.Xbaro;
                else
                    obj.Xo = varargin{4};
                    obj.Xbaro = varargin{4};
                end
            end
            
            if nargin >= 5,
                if isstruct(varargin{5}),
                    obj.Po = varargin{5}.Po;
                    obj.Pbaro = varargin{5}.Pbaro;
                else
                    obj.Po = varargin{5};
                    obj.Pbaro = varargin{5};
                end
                if isempty(obj.Po)
                    obj.Po = 1/eps*eye(size(obj.Xo));
                end
                if isempty(obj.Pbaro)
                    obj.Pbaro = 1/eps*eye(size(obj.Xo));
                end
            elseif nargin >= 4,
                obj.Po = 1/eps*eye(size(obj.Xo));
                obj.Pbaro = obj.Po;
            end
            
            if nargin >=6,
                obj.options = varargin{6};
            else
                obj.options = setOdtbxOptions('OdeSolvOpts',odeset);
            end
            
            if nargin >= 7,
                if all(isfield(varargin{7}, {'tru','est'}))
                    obj.dynarg = varargin{7};
                else
                    obj.dynarg.tru = varargin{7};
                    obj.dynarg.est = varargin{7};
                end
            elseif nargin >= 4,
                obj.dynarg.tru = [];
                obj.dynarg.est = [];
            end
            
            if nargin >= 8,
                if all(isfield(varargin{8}, {'tru','est'}))
                    obj.datarg = varargin{8};
                else
                    obj.datarg.tru = varargin{8};
                    obj.datarg.est = varargin{8};
                end
            elseif nargin >= 4,
                obj.datarg.tru = [];
                obj.datarg.est = [];
            end
        end
    
        function varargout = run_estimator(obj)
            
            % Preallocate variables
            state_component_length = length(obj.Xo);
            num_time_steps = length(obj.tspan);
            
            X_state = NaN(state_component_length*2,1);
            obj.X = NaN(state_component_length,num_time_steps); % X and Xhat could just be consolidated into X_state, but that makes for some long, ugly code whenever you have to access them
            obj.Xhat = NaN(state_component_length,num_time_steps);
            obj.Phat = NaN([size(obj.Pbaro),num_time_steps]);
            obj.y = NaN(state_component_length,num_time_steps);     % Measurement innovations
            obj.eflag = NaN(state_component_length*2,num_time_steps);
%             obj.Y = NaN(state_component_length,num_time_steps);     % True measurements
            obj.Pdy = NaN(state_component_length,state_component_length,num_time_steps); % Measurement innovations covariance
            
%             monteseed = getOdtbxOptions(obj.options, 'MonteCarloSeed', NaN);
%             monteseed_use = NaN(1,1); % Pre-allocate array
%             monteseed_use = monteseed;

%             eratio = getOdtbxOptions(obj.options, 'EditRatio', []) % Default is empty, meaning no meas. editing
%             eflag_set  = getOdtbxOptions(obj.options, 'EditFlag', []) % Default is empty (no meas. editing)


            for sim_time_ind = 1:length(obj.tspan)
                if (sim_time_ind == 1)
                    % Filter initial conditions
                    obj.Xhat(:,1) = obj.Xbaro;
                    obj.X(:,1) = obj.Xo;
                    obj.Phat(:,:,1) = obj.Pbaro;
                    
                    % True covariance
                    obj.Pa(:,:,1) = obj.Po;
                    obj.Pv(:,:,1) = zeros(size(obj.Po));
                    obj.Pw(:,:,1) = zeros(size(obj.Po));
                    obj.Pm(:,:,1) = zeros(size(obj.Po));
                    obj.P(:,:,1) = obj.Pa(:,:,1);
                    % Assumed covariance
                    obj.Phata(:,:,1) = obj.Pbaro;
                    obj.Phatv(:,:,1) = zeros(size(obj.Pbaro));
                    obj.Phatw(:,:,1) = zeros(size(obj.Pbaro));
                    obj.Phatm(:,:,1) = zeros(size(obj.Pbaro));
                    
%                     obj.X(:,1) = obj.Xo + covsmpl(obj.Po, 1, monteseed_use);
                    
                end
                              
                % Calculate true/estimated measurement at tspan(sim_time_ind)
                Yref = ...
                    feval(obj.datfun.tru,obj.tspan(sim_time_ind),obj.X(:,sim_time_ind),obj.datarg.tru);
                Ybar = ...
                    feval(obj.datfun.est,obj.tspan(sim_time_ind),obj.Xhat(:,sim_time_ind),obj.datarg.est);
                [~,Href(:,:),R(:,:)] = ...
                    ominusc(obj.datfun.tru,obj.tspan(sim_time_ind),obj.X(:,sim_time_ind),Yref,obj.options,[],obj.datarg.tru);
                [~,Hsref(:,:),Rhat(:,:)] = ...
                    ominusc(obj.datfun.est,obj.tspan(sim_time_ind),obj.Xhat(:,sim_time_ind),Ybar,obj.options,[],obj.datarg.est);
%                 if (sim_time_ind == 1) % Make Y the right size
%                     obj.Y = NaN(size(Yref),num_time_steps);
%                 end
%                 covsmpl(R(:,:,sim_time_ind))
                obj.Y(:,sim_time_ind) = Yref + covsmpl(R(:,:));
                num_measurements = size(obj.Y(:,1));
                isel = 1:num_measurements;
                
                % Perform measurement update
%                 [obj.Xhat(:,sim_time_ind),obj.Phat(:,:,sim_time_ind),obj.eflag(isel,sim_time_ind),obj.y(isel,sim_time_ind),obj.Pdy(isel,isel,sim_time_ind),~] = ...
%                     kalmup(obj.datfun.est, obj.tspan(sim_time_ind),obj.Xhat(:,sim_time_ind),obj.Phat(:,:,sim_time_ind),obj.Y(:,sim_time_ind),...
%                     obj.options,eflag_set,eratio,obj.datarg.est,isel); 
%                 [obj.Xhat(:,sim_time_ind),obj.Phat(:,:,sim_time_ind)] = ...
                    kalmup(obj.datfun.est,obj.tspan(sim_time_ind),obj.Xhat(:,sim_time_ind),obj.Phat(:,:,sim_time_ind),obj.Y(:,sim_time_ind))
                           
                % Prepare for propagation
                done = false;
                opts = odeset('Event',@obj.events);
                
                % Define the time range and starting state for propagation
                prop_begin_time = obj.tspan(sim_time_ind);
                prop_end_time = obj.tspan(sim_time_ind+1);
                X_state(:,1) = [obj.Xhat(:,sim_time_ind); obj.X(:,sim_time_ind)];
                X_state_begin = X_state(:,1);
                
                % Propagation loop
                while ~done
                    % Propagation
                    time_span = [prop_begin_time, prop_end_time];
                    [time_prop, X_state_prop] = integ(obj.wrapperdyn,time_span,X_state_begin,opts,obj.dynarg);

                    % Check for full propagation
                    if (time_prop(end) == time_span(end))
                        % Save the final propagated state back to the state
                        % (will be saved to time sim_time_ind+1)
                        X_state(:,1) = X_state_prop(:,end);
                        done = true;
                    else
                        % Adjust state/covariance based on user-supplied function

                        % Repropagate from where the loop ended (where
                        % the event occurred)
                        prop_begin_time = time_prop(end);
                        X_state_begin = X_state_prop(end);
                    end

                end
                
                % Pull the variables out of the state
                obj.X = X_state(state_component_length+1:state_component_length*2,sim_time_ind+1);
                obj.Xhat = X_state(1:state_component_length,sim_time_ind+1);
                
            end

%             Calculate errors, package data, and output results (state errors,
%             covariance, residuals, etc)

            %% Output results
            if nargout >= 3,
                varargout{1} = obj.t;
                varargout{2} = obj.Xhat;
                varargout{3} = obj.Phat;
            end
            if nargout >= 4,
                varargout{4} = obj.e;
            end
            if nargout >= 5,
                varargout{5} = obj.y;
            end
            if nargout >= 6,
                varargout{6} = obj.Pa;
                varargout{7} = obj.Pv;
                varargout{8} = obj.Pw;
                varargout{9} = obj.Phata;
                varargout{10} = obj.Phatv;
                varargout{11} = obj.Phatw;
            end
            if nargout >= 12,
                varargout{12} = obj.Sig_sa;
            end
            if(nargout >= 13)
                varargout{13} = obj.eflag;
            end
            if nargout >= 14,
                varargout{14} = obj.Pdy;
            end
            if nargout >= 15,
                varargout{15} = obj.Pdyt;
            end
            if nargout >= 16
                varargout{16} = obj.Pm;
                varargout{17} = obj.Phatm;
            end


        end % run_estimator
        
        
        function [xdot,A,Q] = wrapperdyn(obj,t,X,opts)

            [xdot1,A1,Q1] = feval(obj.dynfun.est,t,X(1:6),opts.est);
            [xdot2,A2,Q2] = feval(obj.dynfun.tru,t,X(7:12),opts.tru);

            xdot = [xdot1;xdot2];
            A = blkdiag(A1,A2);
            Q = blkdiag(Q1,Q2);

        end % wrapperdyn
        
        
        function [value,isterminal,direction] = events(obj,t,X)
            % Consider using functions for conditions
            
            % Event 1:
            condition1 = X; %  We chan change this to be anything related to t or y
            terminal1 = 0;
            direction1 = 0;
            
            % Event 2:
            condition2 = t; %  We chan change this to be anything related to t or y
            terminal2 = 0;
            direction2 = 0;
                        
            % Put all the values together to be returned
            value = [condition1; condition2]; % Condition to look for
            isterminal = [terminal1; terminal2]; % Do we need to end the integration at this point?
            direction = [direction1; direction2]; % Is there a direction involved?
        end % events
        
    end % Methods
    
end % Class

