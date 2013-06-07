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
        %% Constructor
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
        
        
        %% Parameter Mapping
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
        
        
        %% External Force Models
        
        function [xDot,A,Q] = extForces(obj,t,x,jatWorld)
            
            % Interface maintains compatibility with jatForces?
            if (strcmpi(obj.external_func, 'jat'))
                [xDot,A,Q] = jatForcesCon(obj,t,x,jatWorld);
            elseif (strcmpi(obj.external_func, 'gmat'))
                [xDot,A,Q] = gmatForces(obj,t,x,jatWorld);
            else
                disp 'External function not recognized.'
            end
        end
        
        
        function [xDot,A,Q] = gmatForces(obj,t,x,jatWorld)
            % Need GMAT interface info
        end
        
        
        %% Measurement Models
%         function [y,H,R] = gsmeas(obj,t,x,options)
%             %% Get values from options
%             gsID         = getOdtbxOptions(options, 'gsID', [] );
%             gsList       = getOdtbxOptions(options, 'gsList', []);
%             gsECEF       = getOdtbxOptions(options, 'gsECEF', []);
%             epoch        = getOdtbxOptions(options, 'epoch', NaN ); %UTC
%             elMin        = getOdtbxOptions(options, 'gsElevationConstraint', 10)*pi/180; % convert from degs to rads
%             uselt        = getOdtbxOptions(options, 'useLightTime', false);
%             useRange     = getOdtbxOptions(options, 'useRange', true );
%             useRangeRate = getOdtbxOptions(options, 'useRangeRate', true );
%             useDoppler   = getOdtbxOptions(options, 'useDoppler', false );
%             useUnit      = getOdtbxOptions(options, 'useUnit', false );
%             useAngles    = getOdtbxOptions(options, 'useAngles', false );
%             Sched        = getOdtbxOptions(options, 'Schedule',[]); %Tracking Schedule
%             numtypes     = useRange + useRangeRate + useDoppler+3*useUnit+2*useAngles;
%             Type         = getOdtbxOptions(options, 'rangeType','2way');
% %             solve        = getOdtbxOptions(options, 'solvefor',[]);
% %             dyn_cons     = getOdtbxOptions(options, 'dynamicConsider',[]);
% %             loc_cons     = getOdtbxOptions(options, 'localConsider',[]);
% 
%             num_sf = length(obj.solve.param);
%             num_dc = length(obj.dyn_cons.param);
%             num_lc = length(obj.loc_cons.param);
%             if isempty(num_dc),num_dc = 0;end
%             if isempty(num_lc),num_lc = 0;end
%             ind_iono = num_sf + num_dc + find(strncmpi(obj.loc_cons.param,'ION',3));
%             ind_trop = num_sf + num_dc + find(strncmpi(obj.loc_cons.param,'TRP',3));
%             if size(x,1)<max(ind_iono)
%                 ind_iono=[];
%             end
%             if size(x,1)<max(ind_trop)
%                 ind_trop=[];
%             end
% 
%             numtypes = useRange + useRangeRate + useDoppler+3*useUnit;
% 
%             if isempty(gsECEF)
%                 if( isempty(gsList) && ~isempty(gsID) ); gsList = createGroundStationList(); end
%                 gsECEF = zeros(3,length(gsID));
%                 for n=1:length(gsID)
%                         gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
%                 end
%             end
%             
%             M = size(gsECEF,2) * numtypes;
%             N = length(t);
%             if size(t,1)==N, t=t'; end
% 
%             if isnan(epoch); error('An epoch must be set in the options structure.'); end
% 
% 
%             %% call rrdot with all the options passing straight through
%             y = nan(M,N);
%             H = zeros(M,size(x,1),N);  
%             if uselt
%                 %x2 needs states before and after each state in x1 at time steps
%                 %comparable to the light time delay in order for the interpolation
%                 %within the lightTimeCorrection function to have a high level of
%                 %accuracy.
%                 c    = JATConstant('c')/1000;
%                 ltDT = sqrt(sum(x(1:3,:).^2))/c; %
%                 t2   =  unique([t-ltDT, t, t+ltDT]);
% 
%                 % Get the ground station positions at times t2 (see subfunction below)
%                 x2 = getGSstates(gsECEF,epoch,t2);
% 
%                 % Run rrdotlt for each ground station
%                 for n=1:size(gsECEF,2)
%                     % Check schedule for times when tracking is done
%                     if ~isempty(Sched)
%                         gSch = Sched(n==Sched(:,1),2:3);
%                         tind = [];
%                         for m=1:size(gSch,1)
%                             tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
%                         end
%                     else
%                         tind = 1:length(t);
%                     end
%                     t1  = t(tind);
%                     x1  = x(:,tind);
% 
%                     if ~isempty(t1)
% 
%                         [y1,H1,R,t2_lt,x2_lt] = rrdotlt(t1,x1,t2,x2(:,:,n),options);
% 
%                         % apply the elevation constraint
%                         Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
%                         Ephem.SatCoords   = 'ECI';
%                         Ephem.Epoch       = epoch+t2_lt{1}/86400; %UTC
%                         Ephem.StationInfo = 'ECEF';
%                         Ephem.staPos      = gsECEF(:,n)*1000;
%                         [~,el]            = jatStaAzEl(Ephem);
%                         index0            = find( el < elMin );
%                         if length(x2_lt)==2 %then it was a 2way measurement
%                             Ephem.Epoch = epoch+t2_lt{2}/86400; %UTC
%                             [~,el]      = jatStaAzEl(Ephem);
%                             index1      = find( el < elMin );
%                             index0      = union(index0,index1);
%                         end
%                         y1(:,index0) = NaN;
% 
%                         % combine with results from previous stations
%                         indstart                     = 1 + numtypes*(n-1);
%                         indstop                      = numtypes*n;
%                         y(indstart:indstop, tind)    = y1;
%                         H(indstart:indstop, :, tind) = H1;
%                     end
%                 end
%                 clear R;
%             else
%                 % Get the ground station positions at times t (see subfunction below)
%                 gx = getGSstates(gsECEF,epoch,t);
% 
%                 % Run rrdot for each ground station
%                 for n=1:size(gsECEF,2)
%                     % Check schedule for times when tracking is done
%                     if ~isempty(Sched)
%                         gSch = Sched(n==Sched(:,1),2:3);
%                         tind = [];
%                         for m=1:size(gSch,1)
%                             tind = union(tind, find( gSch(m,1)<=t & t<=gSch(m,2) ));
%                         end
%                     else
%                         tind = 1:length(t);
%                     end
%                     if isempty(tind),continue,end
%                     t1 = t(tind);
%                     x1 = x(1:6,tind);
%                     x2 = gx(:,tind,n);
% 
%                     [y1,H1] = rrdotang(t1,x1,x2,options);
% 
%                     % apply the elevation constraint
%                     Ephem.satPos      = x1(1:3,:)*1000; %ECI satellite coordinates (m)
%                     Ephem.SatCoords   = 'ECI';
%                     Ephem.Epoch       = epoch+t1/86400;%UTC
%                     Ephem.StationInfo = 'ECEF';
%                     Ephem.staPos      = gsECEF(:,n)*1000;
%                     [~,el]            = jatStaAzEl(Ephem);        
%                     y1(:,el<elMin)    = NaN;
% 
%                     % combine with results from previous stations
%                     indstart                     = 1 + numtypes*(n-1);
%                     indstop                      = numtypes*n;
%                     y(indstart:indstop, tind)    = y1;
%                     
%                     indstart
%                     indstop
%                     num_sf
% %                     ind_iono(n)
%                     ind_trop(n)
%                     size(H1)
%                     size(H)
%                     
%                     if ~isempty(ind_iono) && ~isempty(ind_trop)
%                         H(indstart:indstop, [1:num_sf ind_iono(n) ind_trop(n)], tind) = H1; 
%                     elseif ~isempty(ind_iono)
%                         H(indstart:indstop, [1:num_sf ind_iono(n)], tind) = H1; 
%                     elseif ~isempty(ind_trop)
%                         H(indstart:indstop, [1:num_sf ind_trop(n)], tind) = H1; 
%                     else
%                         H(indstart:indstop, 1:num_sf, tind) = H1(:,1:num_sf,:); 
%                     end
%                 end
%             end
% 
%             %% double for 2way measurements
%             if strcmpi(Type,'2way')
%                 y=2*y;
%                 H(:,1:num_sf,:)=2*H(:,1:num_sf,:);
%             end
% 
%             %% Add consider parameter data to H
%             bias = strncmpi(obj.loc_cons.param,'MEASBI',6);
%             numbias=sum(bias);
%             mbi = find(bias)+ num_sf + num_dc;
% 
%             if max(mbi)<=size(x,1) % NOT ROBUST SOLUTION!!!!!!!!!!!!!!
%                 y(1:2:end,:) = y(1:2:end,:) + x(mbi,:);
%                 Hbias = zeros(size(H,1),3);
%                 Hbias(1:2:5,:) = diag(ones(1,numbias));
%                 H(:,mbi,:)=Hbias(:,:,ones(1,size(H,3)));
%             end
% 
% 
% 
%             %% Set the measurement covariance output
%             if nargout > 2,
%                 %get sigma out of the options
%                 sigmaDefault = ones(1,size(y,1))*1e-3;
%                 sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
%                 if( length(sigma)~=length(sigmaDefault) )
%                     disp('WARNING: In gsmeas, length(rSigma) does not match total number of measurements');
%                     disp('   Setting rSigma = default (1e-3) for all measurements');
%                     sigma = sigmaDefault;
%                 end
%                 sigma           = diag(sigma);
%                 R = repmat(sigma.^2,[1,1,N]);
%             end
% 
% 
%             %% Subfunction for getting ECI states of ground stations at times t
%             function gx = getGSstates(gsECEF,epoch,t)
% 
%                 D = jatDCM('ecef2eci', epoch+t/86400);
%                 w = [0;0;JATConstant('wEarth')];
% 
%                 gx = zeros(6,length(t),size(gsECEF,2));
%                 for nt = 1:length(t);
%                     M          = rotransf(-D(:,:,nt)*w,D(:,:,nt));
%                     gx(:,nt,:) = M(1:6,1:6)*[gsECEF;zeros(size(gsECEF))];
%                 end
%             end
%         end
%         
    end
    
end

