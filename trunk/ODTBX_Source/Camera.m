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
        xmax;   % Maximum visible distance in the focal plane (X)
        ymax;   % Maximum visible distance in the focal plane (Y)
        theta;  % Attitude Error [rad] (uses small angle approximation)
        bias;   % Optical bias parameter
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
            
            % Calculate maximum visible distance in focal plane
            obj.xmax = f*tan(FOV/2);
            obj.ymax = obj.xmax;
            
            % Assume zero attitude error & bias
            obj.theta = [0;0;0];
            obj.bias = [0;0];
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
                Z = obj.C_CI*(xi(:,i) - obj.R);
                E = [1 obj.theta(3) -obj.theta(2);
                     -obj.theta(3) 1 obj.theta(1);
                     obj.theta(2) -obj.theta(1) 1];
                A = E*Z;
                
                % Calculate landmark position in the focal plane
                x(i) = obj.f*A(1)/A(3) + obj.bias(1);
                y(i) = obj.f*A(2)/A(3) + obj.bias(2);
            end
            
            % Determine landmark visibility
            vis = isVisible(obj,x,y,xi,ni,dt);
            x(~vis) = NaN;
            y(~vis) = NaN;
        end
        
        %%
        % *getCameraPartials(obj,dt)* Returns the measurement partials
        % matrix H for measurements generated at time dt
        function [H Hf Ha Hb] = getCameraPartials(obj,dt)
            
            O = zeros(1,3);
            
            % Get inertial landmark positions
            xi = getInertialLmk(obj.body,dt);
            
            % Calculate the measurement partials for each landmark
            for i=2*size(xi,2):-2:2
                Z = obj.C_CI*(xi(:,i/2) - obj.R);
                E = [1 obj.theta(3) -obj.theta(2);
                     -obj.theta(3) 1 obj.theta(1);
                     obj.theta(2) -obj.theta(1) 1];
                A = E*Z;
                dAdr = -E*obj.C_CI;
                H((i-1):i,:) = obj.f/A(3)*[[1 0 -A(1)/A(3)]*dAdr O;
                                           [0 1 -A(2)/A(3)]*dAdr O];
                
                % Calculate partials with respect to focal length
                if nargout > 1
                    Hf((i-1):i,:) = 1/A(3)*[A(1);
                                            A(2)];
                end
                
                % Calculate partials with respect to attitude error
                if nargout > 2
                    dAdth = [0 -Z(3) Z(2);
                             Z(3) 0 -Z(1);
                             -Z(2) Z(1) 0];
                    Ha((i-1):i,:) = obj.f/A(3)*[1 0 -A(1)/A(3);
                                                0 1 -A(2)/A(3)]*dAdth;
                end
                
                if nargout > 3
                    Hb((i-1):i,:) = eye(2,2);
                end
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
                
                % Determine visibility based on occultation
                vis(i) = dot(ni(:,i),LOS)>0 & abs(x(i))<obj.xmax & ...
                    abs(y(i))<obj.ymax;
                
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

