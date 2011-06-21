classdef Camera < handle
    %CAMERA Class to model an onboard imaging device used for optical navigation
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
    %   Kenneth Getzandanner    03/29/2011      Original Camera.m

    %% Properties
    properties
        body;   % Target celestial body
        R;      % Inertial position
        C_CI;   % DCM from inertial to camera frame
        f;      % Camera focal length
        FOV;    % Camera Field of View [rad]
        shadow; % Shadow landmark visibility {true|false}
    end
    
    %% Methods
    methods
        
        %%
        % *Camera(body,R,C_CI,f,FOV,shadow)* [constructor]
        function obj = Camera(body,R,C_CI,f,FOV,shadow)
            obj.body = body;
            obj.R = R;
            obj.C_CI = C_CI;
            obj.f = f;
            obj.FOV = FOV;
            obj.shadow = shadow;
        end
        
        %% 
        % *getCameraFrame(obj,dt)* Returns the position of landmarks in the 
        % camera focal plane at time dt
        function [x y] = getCameraFrame(obj,dt)
            
            % Get inertial landmark positions
            xi = getInertialLmk(obj.body,dt);
            
            % Get inertial landmark normal vector
            ni = getInertialLmkNormal(obj.body,dt);
            
            for i=size(xi,2):-1:1
                
                % Rotate landmark position into camera frame
                A = obj.C_CI*(xi(:,i) - obj.R);
                
                % Calculate landmark position in the focal plane
                x(i) = obj.f*A(1)/A(3);
                y(i) = obj.f*A(2)/A(3);
            end
            
            % Determine landmark visibility
            vis = isVisible(obj,x,y,xi,ni,dt);
            x(~vis) = NaN;
            y(~vis) = NaN;
        end
        
        %%
        % *getCameraPartials(obj,dt)* Returns the measurement partials
        % matrix H for measurements generated at time dt
        function H = getCameraPartials(obj,dt)
            
            O = zeros(1,3);
            
            % Get inertial landmark positions
            xi = getInertialLmk(obj.body,dt);
            
            % Calculate the measurement partials for each landmark
            for i=2*size(xi,2):-2:2
                A = obj.C_CI*(xi(:,i/2) - obj.R);
                dAdr = -obj.C_CI;
                H((i-1):i,:) = obj.f/A(3)*[[1 0 -A(1)/A(3)]*dAdr O;
                                           [0 1 -A(2)/A(3)]*dAdr O];
            end
        end
        
        %%
        % *isVisible(obj,x,y,xi,ni,dt)* Determine if a landmark is visible
        % based on target body occultation and shadowing
        function vis = isVisible(obj,x,y,xi,ni,dt)
            
            % Loop through each measurements
            for i=size(xi,2):-1:1
                
                % Calculate the line-of-sight vector
                LOS = obj.R-xi(:,i);
                
                % Calculate maximum visible distance in focal plane
                max = obj.f*tan(obj.FOV/2);
                
                % Determine visibility based on occultation
                vis(i) = dot(ni(:,i),LOS)>0 & abs(x(i))<max & abs(y(i))<max;
                
                % If selected, calculate shadow visibility
                if obj.shadow
                    sunVec = getState(obj.body,dt,'J2000','SUN');
                    sunVec = -unit(sunVec(1:3));
                    sunVis(i) = dot(ni(:,i),sunVec)>0;
                end
            end
            
            % Combine visibility
            if obj.shadow
                vis = sunVis & vis;
            end
        end
        
        %%
        % *plotCameraFrame(obj,dt)* Simulate OpNav image by plotting
        % landmark positions in the focal plane
        function plotCameraFrame(obj,dt)
            [x y] = getCameraFrame(obj,dt);
            max = obj.f*tan(obj.FOV/2);
            plot(x,y,'*','MarkerSize',12);
            axis([-max max -max max])
            set(gca,'xdir','reverse')
        end
    end
    
end

