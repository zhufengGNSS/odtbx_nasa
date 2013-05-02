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
        function obj = solve_consider(varargin)
            % Provides the opportunity to define solve for/consider
            % params when the class is created.
            if nargin >= 1
                obj.solve.param = varargin(1);
            end
            if nargin >= 2
                obj.dyn_cons.param = varargin(2); 
            end
            if nargin >= 3
                obj.loc_cons.param = varargin(3);
            end
            if nargin >= 4
                obj.external_func = varargin(4);
            end
            
            % Even if variables weren't passed in, these will occur if
            % there is data in the variables
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
                [xDot,A,Q] = obj.jatForces(obj,t,x,jatWorld);
            elseif (obj.external_func == 'gmat')
                [xDot,A,Q] = obj.gmatForces(obj,t,x,jatWorld);
            else
                disp 'External function not recognized.'
            end
        end
        
        
        function [xDot,A,Q] = obj.jatForces(obj,t,x,jatWorld)
            % John Gaebler's Code, revised
            [nx,nt] = size(x);
            xDot = zeros(nx,nt);
            A = zeros(nx,nx,nt);
            Q = zeros(nx,nx,nt);

            [Fei,Fsi,Fmi,Fsri]=deal([]);

            if exist('dyn_cons','var')
                solv = length(solve);
                Fei =solv+find(strcmpi(dyn_cons,'EARTH-GM'));
                Fsi =solv+find(strcmpi(dyn_cons,'SOLAR-GM'));
                Fmi =solv+find(strcmpi(dyn_cons,'LUNAR-GM'));
                Fsri=solv+find(strcmpi(dyn_cons,'SOLRAD 8'));

                Fearth=x(Fei,1);
                Fsun=x(Fsi,1);
                Fmoon=x(Fmi,1);
                Fsrp=x(Fsri,1); % Fraction of Cr
            end

            x=x(1:6,:)*1000; % conversion km -> m

            for i=1:nt
                xDot(1:6,i) = jatWorld.derivs(t(i),x(:,i))/1000;  % conversion m -> km
            end

            SRPjat=nan(3,nt);
            acc_earth=nan(3,nt);
            acc_sun=nan(3,nt);
            acc_moon=nan(3,nt);

            if( nargout > 1 )
                A = [];

                useGravPartial = ( jatWorld.get_j2_gravity_flag() || jatWorld.get_2body_gravity_flag() );
                useDragPartial = jatWorld.get_drag_flag();

                for i=1:nt
                    if useGravPartial
                        A = zeros(nx,nx,nt);    
                        A(1:3,4:6,i) = eye(3);
                        if( useGravPartial )
                            A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.gravityPartials(x(1:3,i));
                        end
                    end

                    % JAT Force List order:
                    % 0 Earth Gravity
                    % 1 Solar Gravity
                    % 2 Lunar Gravity
                    % 3 Drag
                    % 4 Solar Radiation Force
                    %
                    % If a force is not added, all higher slots are shifted.
                    cur = 0;
                    if Fei % earth gravity error
                        % need spacecraft position and earth GM
                        acc_earth(:,i)=jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                        A(4:6,Fei,i)=A(4:6,Fei,i) + acc_earth(:,i);
                        cur=cur+1;
                    end
                    if Fsi % solar gravity error
                        % need spacecraft inertial position, sun inertial position,
                        % and sun GM
                        acc_sun(:,i)=jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                        A(4:6,Fsi,i)=A(4:6,Fsi,i) + acc_sun(:,i);
                        cur=cur+1;
                    end
                    if Fmi % lunar gravity error
                        % need spacecraft inertial position, lunar inertial
                        % position, and lunar GM
                        acc_moon(:,i)=jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                        A(4:6,Fmi,i)=A(4:6,Fmi,i) + acc_moon(:,i);
                        cur=cur+1;
                    end
                    if useDragPartial 
                        A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.dragPartials();
                        cur=cur+1;
                    end
                    if Fsri
                        SRPjat(:,i)=jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray;
                        A(4:6,Fsri,i)=A(4:6,Fsri,i) + (SRPjat(:,i))/1000;% m -> km          +[-3.05073589;6.75748831;2.03244886]*1e-9
                        cur=cur+1;
                    end
                end
            end

            if nargout == 3,
                Q = zeros(nx,nx,nt);
                Q(4:6,4:6,:) = repmat(diag([1e-9 1e-9 1e-9].^2),[1 1 nt])/1000^2;
                Q(Fei,Fei,end) = 0;
                Q(Fsi,Fsi,end) = 0;
                Q(Fmi,Fmi,end) = 0;
                Q(Fsri,Fsri,end) = 0;
            end
        end
        
        
        function [xDot,A,Q] = obj.gmatForces(obj,t,x,jatWorld)
            % Need GMAT interface info
        end
        
    end
    
end

