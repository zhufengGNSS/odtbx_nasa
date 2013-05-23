classdef solve_consider
    %solve_consider Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solve = struct('param', {''}, 'user_order', {''});
        dyn_cons = struct('param', {''}, 'user_order', {''});
        loc_cons = struct('param', {''}, 'user_order', {''});
        param_order = containers.Map;
        external_func = 'jat';
    end
    
    
    methods
        function obj = solve_consider(varargin)
            % Provides the opportunity to define solve for/consider
            % params when the class is created.
            if nargin >= 1 % First variable is always the solve param
                obj.solve.param = varargin{1};
            end
            if nargin == 2 % If length is two, short state, second is external function
                obj.external_func = varargin(2);
            end
            if nargin > 2 % If length is bigger than two, we'll get a full state
                obj.dyn_cons.param = varargin{2}; 
            end
            if nargin >= 3 % Full state
                obj.loc_cons.param = varargin{3};
            end
            if nargin >= 4 % Full state
                obj.external_func = varargin(4);
            end
            
            % containers.Map objects are treated as handles class, so every
            % instance of the solve_consider class point back to the same
            % Map. 
            % Because there is no easy way to do deep copy (clone) of a
            % Map, we have to use a little trick inspired by:
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/278249
            % It's possible to force Matlab to recognize a containers.Map
            % object as value-based if it's returned from a function.
            obj.param_order = map_params(obj);
        end
        
        
        function temp_map = map_params(obj)
            temp_map = containers.Map();
            % Maps all of the consider functions to values representing
            % the order they should be listed in the A matrix.
            if (~isempty(obj.solve.param))
                for order = 1:length(obj.solve.param)
                    % Assign the value
                    obj.solve.user_order{order,1} = order;
                    % Add to map
                    temp_map(obj.solve.param{order}) = ...
                        obj.solve.user_order{order,1};
                end
            end  
            current_len = length(obj.solve.param);
            if (~isempty(obj.dyn_cons.param))
                for order = 1:length(obj.dyn_cons.param)
                    % Assign the value
                    obj.dyn_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    temp_map(obj.dyn_cons.param{order}) = ...
                        obj.dyn_cons.user_order{order,1};
                end
            end
            current_len = current_len + length(obj.dyn_cons.param);
            if (~isempty(obj.loc_cons.param))
                for order = 1:length(obj.loc_cons.param)
                    % Assign the value
                    obj.loc_cons.user_order{order,1} = order + current_len;
                    % Add to map
                    temp_map(obj.loc_cons.param{order}) = ...
                        obj.loc_cons.user_order{order,1};
                end
            end

            % temp_map is returned from this function to overcome
            % shortcomings of containers.Map objects.
        end
        
        % ----------------------------------------------------------------
        % The following functions currently are written into the class
        % definition, but there's a case to be made for removing them and
        % just making them separate files. As long as they have access to
        % the class variables, they should be fine. That would change the
        % inputs and outputs of the files, though, as the object would have
        % to be explicitly passed.
        % ----------------------------------------------------------------
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            
            % Interface maintains compatibility with jatForces?
            if (strcmpi(obj.external_func, 'jat'))
                [xDot,A,Q] = obj.jatForces(t,x,jatWorld);
            elseif (strcmpi(obj.external_func, 'gmat'))
                [xDot,A,Q] = obj.gmatForces(t,x,jatWorld);
            else
                disp 'External function not recognized.'
            end
        end
        
        
        function [xDot,A,Q] = jatForces(obj,t,x,jatWorld)
%             x
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
        
        
        function [xDot,A,Q] = gmatForces(obj,t,x,jatWorld)
            % Need GMAT interface info
        end
        
        
        function [y,H,R] = gsmeas(obj,t,x,options)
            %% Get values from options
            gsID         = getOdtbxOptions(options, 'gsID', [] );
            gsList       = getOdtbxOptions(options, 'gsList', []);
            gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
            epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
            elMin        = getOdtbxOptions(options, 'gsElevationConstraint', 10)*pi/180; % convert from degs to rads
            uselt        = getOdtbxOptions(options, 'useLightTime', false);
            useRange     = getOdtbxOptions(options, 'useRange', true );
            useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
            useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
            useUnit      = getOdtbxOptions(options, 'useUnit', false );
            useAngles    = getOdtbxOptions(options, 'useAngles', false );
            Sched        = getOdtbxOptions(options, 'Schedule',[]); %Tracking Schedule
            numtypes     = useRange + useRangeRate + useDoppler+3*useUnit+2*useAngles;
            Type         = getOdtbxOptions(options, 'rangeType','2way');
%             solve        = getOdtbxOptions(options, 'solvefor',[]);
%             dyn_cons     = getOdtbxOptions(options, 'dynamicConsider',[]);
%             loc_cons     = getOdtbxOptions(options, 'localConsider',[]);

            num_sf = length(obj.solve.param);
            num_dc = length(obj.dyn_cons.param);
            num_lc = length(obj.loc_cons.param);
            if isempty(num_dc),num_dc = 0;end
            if isempty(num_lc),num_lc = 0;end
            ind_iono = num_sf + num_dc + find(strncmpi(loc_cons,'ION',3));
            ind_trop = num_sf + num_dc + find(strncmpi(loc_cons,'TRP',3));
            if size(x,1)<max(ind_iono)
                ind_iono=[];
            end
            if size(x,1)<max(ind_trop)
                ind_trop=[];
            end

            numtypes     = useRange + useRangeRate + useDoppler+3*useUnit;

            if isempty(gsECEF)
                if( isempty(gsList) && ~isempty(gsID) ); gsList = createGroundStationList(); end
                gsECEF = zeros(3,length(gsID));
                for n=1:length(gsID)
                        gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
                end
            end
            M            = size(gsECEF,2) * numtypes;
            N            = length(t);
            if size(t,1)==N, t=t'; end

            if isnan(epoch); error('An epoch must be set in the options structure.'); end


            %% call rrdot with all the options passing straight through
            y    = nan(M,N);
            H    = zeros(M,size(x,1),N);  
            if uselt
                %x2 needs states before and after each state in x1 at time steps
                %comparable to the light time delay in order for the interpolation
                %within the lightTimeCorrection function to have a high level of
                %accuracy.
                c    = JATConstant('c')/1000;
                ltDT = sqrt(sum(x(1:3,:).^2))/c; %
                t2   =  unique([t-ltDT, t, t+ltDT]);

                % Get the ground station positions at times t2 (see subfunction below)
                x2 = getGSstates(gsECEF,epoch,t2);

                % Run rrdotlt for each ground station
                for n=1:size(gsECEF,2)
                    % Check schedule for times when tracking is done
                    if ~isempty(Sched)
                        gSch = Sched(n==Sched(:,1),2:3);
                        tind = [];
                        for m=1:size(gSch,1)
                            tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
                        end
                    else
                        tind = 1:length(t);
                    end
                    t1  = t(tind);
                    x1  = x(:,tind);

                    if ~isempty(t1)

                        [y1,H1,R,t2_lt,x2_lt] = rrdotlt(t1,x1,t2,x2(:,:,n),options);

                        % apply the elevation constraint
                        Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
                        Ephem.SatCoords   = 'ECI';
                        Ephem.Epoch       = epoch+t2_lt{1}/86400; %UTC
                        Ephem.StationInfo = 'ECEF';
                        Ephem.staPos      = gsECEF(:,n)*1000;
                        [~,el]            = jatStaAzEl(Ephem);
                        index0            = find( el < elMin );
                        if length(x2_lt)==2 %then it was a 2way measurement
                            Ephem.Epoch = epoch+t2_lt{2}/86400; %UTC
                            [~,el]      = jatStaAzEl(Ephem);
                            index1      = find( el < elMin );
                            index0      = union(index0,index1);
                        end
                        y1(:,index0) = NaN;

                        % combine with results from previous stations
                        indstart                     = 1 + numtypes*(n-1);
                        indstop                      = numtypes*n;
                        y(indstart:indstop, tind)    = y1;
                        H(indstart:indstop, :, tind) = H1;
                    end
                end
                clear R;
            else
                % Get the ground station positions at times t (see subfunction below)
                gx = getGSstates(gsECEF,epoch,t);

                % Run rrdot for each ground station
                for n=1:size(gsECEF,2)
                    % Check schedule for times when tracking is done
                    if ~isempty(Sched)
                        gSch = Sched(n==Sched(:,1),2:3);
                        tind = [];
                        for m=1:size(gSch,1)
                            tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
                        end
                    else
                        tind = 1:length(t);
                    end
                    if isempty(tind),continue,end
                    t1 = t(tind);
                    x1 = x(1:6,tind);
                    x2 = gx(:,tind,n);

                    [y1,H1] = rrdotang(t1,x1,x2,options);

                    % apply the elevation constraint
                    Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
                    Ephem.SatCoords   = 'ECI';
                    Ephem.Epoch       = epoch+t1/86400;%UTC
                    Ephem.StationInfo = 'ECEF';
                    Ephem.staPos      = gsECEF(:,n)*1000;
                    [~,el]            = jatStaAzEl(Ephem);        
                    y1(:,el<elMin)    = NaN;

                    % combine with results from previous stations
                    indstart                     = 1 + numtypes*(n-1);
                    indstop                      = numtypes*n;
                    y(indstart:indstop, tind)    = y1;
                    if ~isempty(ind_iono) && ~isempty(ind_trop)
                        H(indstart:indstop, [1:num_sf ind_iono(n) ind_trop(n)], tind) = H1; 
                    elseif ~isempty(ind_iono)
                        H(indstart:indstop, [1:num_sf ind_iono(n)], tind) = H1; 
                    elseif ~isempty(ind_trop)
                        H(indstart:indstop, [1:num_sf ind_trop(n)], tind) = H1; 
                    else
                        H(indstart:indstop, 1:num_sf, tind) = H1(:,1:num_sf,:); 
                    end
                end
            end

            %% double for 2way measurements
            if strcmpi(Type,'2way')
                y=2*y;
                H(:,1:num_sf,:)=2*H(:,1:num_sf,:);
            end

            %% Add consider parameter data to H
            bias = strncmpi(loc_cons,'MEASBI',6);
            numbias=sum(bias);
            mbi = find(bias)+ num_sf + num_dc;

            if max(mbi)<=size(x,1) % NOT ROBUST SOLUTION!!!!!!!!!!!!!!
                y(1:2:end,:) = y(1:2:end,:) + x(mbi,:);
                Hbias = zeros(size(H,1),3);
                Hbias(1:2:5,:) = diag(ones(1,numbias));
                H(:,mbi,:)=Hbias(:,:,ones(1,size(H,3)));
            end



            %% Set the measurement covariance output
            if nargout > 2,
                %get sigma out of the options
                sigmaDefault = ones(1,size(y,1))*1e-3;
                sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
                if( length(sigma)~=length(sigmaDefault) )
                    disp('WARNING: In gsmeas, length(rSigma) does not match total number of measurements');
                    disp('   Setting rSigma = default (1e-3) for all measurements');
                    sigma = sigmaDefault;
                end
                sigma           = diag(sigma);
                R = repmat(sigma.^2,[1,1,N]);
            end


            %% Subfunction for getting ECI states of ground stations at times t
            function gx = getGSstates(gsECEF,epoch,t)

                D = jatDCM('ecef2eci', epoch+t/86400);
                w = [0;0;JATConstant('wEarth')];

                gx = zeros(6,length(t),size(gsECEF,2));
                for nt = 1:length(t);
                    M          = rotransf(-D(:,:,nt)*w,D(:,:,nt));
                    gx(:,nt,:) = M(1:6,1:6)*[gsECEF;zeros(size(gsECEF))];
                end
            end
        end
        
    end
    
end

