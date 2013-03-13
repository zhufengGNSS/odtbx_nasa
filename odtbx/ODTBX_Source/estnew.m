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
            
            
            obj.events_fcn = @obj.events_default;
            obj.control_events_fcn = @obj.control_events_default;
        end
    
        
        function varargout = run_estimator(obj)
            
            % Preallocate variables
            state_component_length = length(obj.Xo);
            num_time_steps = length(obj.tspan);
            
            % X and Xhat could just be consolidated into X_state, 
            % but that makes for some long, ugly code whenever you have 
            % to access them. We use a hybrid approach: Inputs and outputs
            % to the class will be broken down
            X_state = NaN(state_component_length*2,1);
            obj.X = NaN(state_component_length,num_time_steps); 
            obj.Xhat = NaN(state_component_length,num_time_steps);
            obj.Phat = NaN([size(obj.Pbaro),num_time_steps]);
                        
            monteseed = getOdtbxOptions(obj.options, 'MonteCarloSeed', NaN);
            monteseed_use = monteseed;
            
            eratio = getOdtbxOptions(obj.options, 'EditRatio', []); % Default is empty, meaning no meas. editing
            eflag_set  = getOdtbxOptions(obj.options, 'EditFlag', []); % Default is empty (no meas. editing)

            for sim_time_ind = 1:length(obj.tspan)-1
                if (sim_time_ind == 1)
                    % Filter initial conditions
                    obj.Xhat(:,1) = obj.Xbaro;
                    obj.X(:,1) = obj.Xo + covsmpl(obj.Po, 1, monteseed_use);
                    obj.Phat(:,:,1) = obj.Pbaro;
                    
                    X_state = [obj.Xhat(:,sim_time_ind); obj.X(:,sim_time_ind)];
                    
                    % Find size of measurements
                    Y_size = length(feval(obj.datfun.tru,obj.tspan(sim_time_ind), ...
                        X_state(:,sim_time_ind),obj.datarg.tru));
                    R = NaN(2*Y_size);
                end
                           
                % Prepare for propagation
                done = false;
                
                if (sim_time_ind ~= length(obj.tspan)) % Only propagates when there's something left to propagate

                    % Define the time range and starting state for propagation
                    prop_begin_time = obj.tspan(sim_time_ind);
                    prop_end_time = obj.tspan(sim_time_ind+1);
                    X_state(:,1) = [obj.Xhat(:,sim_time_ind); obj.X(:,sim_time_ind)];
                    X_state_begin = X_state(:,1);

                    % Propagation loop
                    while ~done
                        % Propagation
                        time_span = [prop_begin_time, prop_end_time];
                        [time_prop, X_state_prop, time_event, X_event, Phi_state_prop, Phi_event, S_state_prop, S_event] = integev(@obj.wrapperdyn,time_span,X_state_begin,[],obj.dynarg,@obj.events_fcn);
             
                        % Check for full propagation
                        if (time_prop(end) == time_span(end))
                            % Save the final propagated state back to the state
                            % (will be saved to time sim_time_ind+1)
                            X_state(:,1) = X_state_prop(:,end);
                            Phi_state(:,:,1) = Phi_state_prop(1:state_component_length,1:state_component_length,end);
                            S_state(:,:,1) = S_state_prop(1:state_component_length,1:state_component_length,end);
                            done = true;
                        else
                            S_event(:,:,1) = (S_event(:,:,end) + S_event(:,:,end)')/2;
                            Phi_event(:,:,end) = Phi_event(:,:,end)*obj.Phat(:,:,sim_time_ind)*Phi_event(:,:,end)' + S_event(:,:,end);
                            Phi_event(:,:,end) = (Phi_event(:,:,end) + Phi_event(:,:,end)')/2;
                            
                            % Adjust state/covariance based on user-supplied function
                            [X_new, P_new] = feval(@obj.control_events_fcn, time_event(:,end), X_event(:,end), Phi_event(:,:,end));

                            % What do we need to do with the new P? Use it
                            % to update Xhat?
                            
                            % Repropagate from where the loop ended (where
                            % the event occurred)
                            prop_begin_time = time_event(end);
                            X_state_begin = X_new;
                            
                        end

                    end

                    % Pull the variables out of the state to save to object
                    obj.X(:,sim_time_ind+1) = X_state(state_component_length+1:state_component_length*2,1);
                    obj.Xhat(:,sim_time_ind+1) = X_state(1:state_component_length,1);
                    obj.t(sim_time_ind+1,1) = time_prop(end);
                    
                    S_state(:,:,1) = (S_state(:,:,1) + S_state(:,:,1)')/2;
                    obj.Phat(:,:,sim_time_ind+1) = Phi_state(:,:,1)*obj.Phat(:,:,sim_time_ind)*Phi_state(:,:,1)' + S_state(:,:,1);
                    obj.Phat(:,:,sim_time_ind+1) = (obj.Phat(:,:,sim_time_ind+1) + obj.Phat(:,:,sim_time_ind+1)')/2;
                    
                end
                
                % Calculate true/estimated measurement at tspan(sim_time_ind)
                [Y_state(1:Y_size*2,1),R(:,:)] = ...
                    obj.wrappermeas(obj.datfun,obj.tspan(sim_time_ind+1), ...
                    X_state(:,1),obj.datarg,obj.options);
                
                % Apply measurement errors
                obj.Y(:,sim_time_ind+1) = Y_state(Y_size+1:Y_size*2,1) + ...
                    covsmpl(R(Y_size+1:Y_size*2,Y_size+1:Y_size*2));
                
                % Perform measurement update
                [obj.Xhat(:,sim_time_ind+1),obj.Phat(:,:,sim_time_ind+1), eflag_set] = ...
                    kalmup(obj.datfun.est,obj.tspan(sim_time_ind+1),...
                    obj.Xhat(:,sim_time_ind+1),obj.Phat(:,:,sim_time_ind+1),...
                    obj.Y(:,sim_time_ind+1),obj.options, eflag_set, eratio,...
                    obj.datarg.est);
            end

%             Calculate errors, package data, and output results (state errors,
%             covariance, residuals, etc)
            obj.e = obj.Xhat(:,:) - obj.X(:,:);


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
                varargout{5} = obj.Y;
            end

        end % run_estimator
        
        
        function [xdot,A,Q] = wrapperdyn(obj,time,X,opts)
%             
%             [xdot1,A1,Q1] = feval(obj.dynfun.est,time,X(1:6,1),opts.est);
%             [xdot2,A2,Q2] = feval(obj.dynfun.tru,time,X(7:12,1),opts.tru);
            
            [xdot1,A1,Q1] = feval(obj.dynfun.est,time,X,opts.est);
            [xdot2,A2,Q2] = feval(obj.dynfun.tru,time,X,opts.tru);

            xdot = [xdot1;xdot2];
            A = blkdiag(A1,A2);
            Q = blkdiag(Q1,Q2);

        end % wrapperdyn
      
        
        function [Y,R] = wrappermeas(obj, datfun, time, X_state, datarg, options)
            state_size = length(X_state);
            % Ybar
            Ybar = feval(datfun.est,time,X_state(1:state_size/2),datarg.est);
            
            % Yref
            Yref = feval(datfun.tru,time,X_state(state_size/2+1:state_size),datarg.tru);
            
            [~,Hsref(:,:),Rhat(:,:)] = ominusc(datfun.est,time,...
                X_state(1:state_size/2),Ybar,options,[],datarg.est);
            
            [~,Href(:,:),Rref(:,:)] = ominusc(datfun.tru,time,...
                X_state(state_size/2+1:state_size),Yref,options,[],datarg.tru);
            
            % Combine
            Y = [Ybar;Yref];
            R = blkdiag(Rhat,Rref);
            
        end % wrappermeas
        
        
        %% Examples
        % These are what event functions and control functions should look
        % like. 
        
        function [value,isterminal,direction] = events_default(obj,t,X,varargin)
            % This function is used to kick the integrator out of its loop
            % at a certain point (as determined by time or state
            % conditions) in order to perform an action desginated by
            % control_events.
            
            % See header in integev.m for details on event function formats
            
            % Consider using functions for conditions

%             % Event 1:
%             condition1 = X(1); %  We can change this to be anything related to t or X
%             terminal1 = 0;
%             direction1 = 0;
%             
%             % Event 2:
%             condition2 = t - 295; %  We can change this to be anything related to t or X
%             terminal2 = 1;
%             direction2 = 0;
                        
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
    
end % Class

