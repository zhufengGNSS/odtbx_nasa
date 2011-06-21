classdef Attitude < handle
    %ATTITUDE Class to store spacecraft attitude information
    %
    % (This file is part of ODTBX, The Orbit Determination Toolbox, and is
    %  distributed under the NASA Open Source Agreement.  See file source for
    %  more details.)
    
    % ODTBX: Orbit Determination Toolbox
    %
    % Copyright (c) 2003-2011 United States Government as represented by the
    % administrator of the National Aeronautics and Space Administration. All
    % Other Rights Reserved.
    %
    % This file is distributed "as is", without any warranty, as part of the
    % ODTBX. ODTBX is free software; you can redistribute it and/or modify it
    % under the terms of the NASA Open Source Agreement, version 1.3 or later.
    %
    % You should have received a copy of the NASA Open Source Agreement along
    % with this program (in a file named License.txt); if not, write to the
    % NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.
    
    %  REVISION HISTORY
    %   Author                  Date         	Comment
    %   Kenneth Getzandanner    03/29/2011      Original Attitude.m    

    %% Properties
    properties
        tspan;  % Timespan
        C_CI;   % Rotation Matrix
    end
    
    %% Methods
    methods
        
        %%
        % *Attitude(tspan)* [constructor]
        function obj = Attitude(tspan)
            obj.tspan = tspan;
            nt = length(tspan);
            obj.C_CI = repmat(eye(3,3),[1,1,nt]);
        end
        
        %%
        % *nadir(obj,t,r)* Simple function to assign nadir pointing
        % attitude at time t given a position vector r
        function nadir(obj,t,r)
            
            % Determine number of timesteps
            nt = length(t);
            
            % Define a reference vector
            ref = repmat([0;0;1],1,nt);
            
            % Assign unit vectors
            l = unit(-r);
            m = unit(cross(ref,l));
            n = unit(cross(l,m));
            
            % Create DCM
            DCM(1,:,:) = m;
            DCM(2,:,:) = n;
            DCM(3,:,:) = l;
            
            % Assign C_CI the value of DCM at each element in t
            for i=1:nt
                obj.C_CI(:,:,obj.tspan == t(i)) = DCM(:,:,i);
            end 
        end 
    end
end