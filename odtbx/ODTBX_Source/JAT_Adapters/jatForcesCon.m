function [xDot,A,Q] = jatForcesCon(obj,t,x,jatWorld)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% John Gaebler's Code, revised
    [nx,nt] = size(x);
    % Preallocate arrays
    dyn_max = length(obj.solve.user_order) + length(obj.dyn_cons.user_order);
    xDot = zeros(nx,nt);
    A = zeros(nx,dyn_max,nt);
    Q = zeros(nx,nx,nt);

    [Fei,Fsi,Fmi,Fsri]=deal([]);

    % Make the code more readable and reduce hash lookups
    if isKey(obj.param_order, 'EARTH-GM')
        Fei = obj.param_order('EARTH-GM');
    end
    if isKey(obj.param_order, 'SOLAR-GM')
        Fsi = obj.param_order('SOLAR-GM');
    end
    if isKey(obj.param_order, 'LUNAR-GM')
        Fmi = obj.param_order('LUNAR-GM');
    end
    if isKey(obj.param_order, 'SOLRAD 8')
        Fsri = obj.param_order('SOLRAD 8');
    end

    x=x(1:6,:)*1000; % conversion km -> m

    % Check the units
    for i=1:nt
        xDot(1:6,i) = jatWorld.derivs(t(i),x(:,i))/1000;  % conversion m -> km
    end

    SRPjat=nan(3,nt);
    acc_earth=nan(3,nt);
    acc_sun=nan(3,nt);
    acc_moon=nan(3,nt);

    if( nargout > 1 )
%                 A = [];
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
                acc_earth(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fei,i) = A(4:6,Fei,i) + acc_earth(:,i);
                cur=cur+1;
            end
            if Fsi % solar gravity error
                % need spacecraft inertial position, sun inertial position,
                % and sun GM
                acc_sun(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fsi,i) = A(4:6,Fsi,i) + acc_sun(:,i);
                cur=cur+1;
            end
            if Fmi % lunar gravity error
                % need spacecraft inertial position, lunar inertial
                % position, and lunar GM
                acc_moon(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray/1000;% m -> km
                A(4:6,Fmi,i) = A(4:6,Fmi,i) + acc_moon(:,i);
                cur=cur+1;
            end
            if useDragPartial 
                A(4:6,1:3,i) = A(4:6,1:3,i) + jatWorld.dragPartials();
                cur=cur+1;
            end
            if Fsri
                SRPjat(:,i) = ...
                    jatWorld.spacetime.getForce(cur).acceleration(jatWorld.spacetime.time,jatWorld.spacetime.earthRef,jatWorld.sc).getArray;
                A(4:6,Fsri,i) = A(4:6,Fsri,i) + (SRPjat(:,i))/1000;% m -> km          +[-3.05073589;6.75748831;2.03244886]*1e-9
                cur=cur+1;
            end

            % If there are any more consider parameters, they
            % should go here

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

