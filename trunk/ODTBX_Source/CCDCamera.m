classdef CCDCamera < Camera
    %CCDCAMERA Class extends the Camera onboard imaging device model to
    %include CCD pixel mapping and distortion effects.
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
    %   Kenneth Getzandanner    07/15/2011      Original CCDCamera.m
    
    %% Properties
    properties
        Rx;     % CCD Resolution (X)
        Ry;     % CCD Resolution (Y)
        Er;     % Distortion Coefficient of Magnitude (R)
        Em;     % Distortion Coefficient of Magnitude (M)
        En;     % Distortion Coefficient of Magnitude (N)
        Kx;     % Pixel Mapping (X)
        Ky;     % Pixel Mapping (Y)
        s0;     % Zero pixel offset (pixel)
        l0;     % Zero pixel offset (line)
    end
    
    %% Methods
    methods
        
        %%
        % *Camera(body,R,C_CI,f,FOV,shadow,E,Res)* [constructor]
        function obj = CCDCamera(body,R,C_CI,f,FOV,shadow,E,Res)
            obj = obj@Camera(body,R,C_CI,f,FOV,shadow);
            obj.Er = E(1);
            obj.Em = E(2);
            obj.En = E(3);
            
            obj.Rx = Res(1);
            obj.Ry = Res(2);
            
            obj.s0 = Res(1)/2;
            obj.l0 = Res(2)/2;
            
            obj.Kx = obj.s0/(2*obj.xmax);
            obj.Ky = -obj.l0/(2*obj.ymax);
        end
        
        %% 
        % *getCameraFrame(obj,dt)* Returns the position of landmarks on the 
        % camera CCD at time dt
        function [s l x y] = getCameraFrame(obj,dt)
            
            % Calculate landmark positions in the camera focal plane
            [x y] = getCameraFrame@Camera(obj,dt);
            
            % Convert from xy positions to pixel/line locations for each
            % landmark
            for i = length(x):-1:1
                A = [x(i)*(x(i)^2+y(i)^2) x(i)*y(i) x(i)^2;
                    y(i)*(x(i)^2+y(i)^2) y(i)^2 x(i)*y(i)];
                
                % Add distortion effects
                P = [x(i);y(i)] + A*[obj.Er;obj.Em;obj.En];
                
                % Map to pixel/line locations
                S = [obj.s0;obj.l0] + [obj.Kx 0;0 obj.Ky]*P;
                
                s(i) = S(1);
                l(i) = S(2);
            end
            
        end
        
        %%
        % *getCameraPartials(obj,dt)* Returns the measurement partials
        % matrix H for measurements generated at time dt
        function H = getCameraPartials(obj,dt)
            
            % Calculate measurement partials with respect to xy
            H = getCameraPartials@Camera(obj,dt);
            [~,~,x y] = getCameraFrame(obj,dt);
            
            % Use the chain rule to calculate measurement partials with
            % respect to pixel/line locations
            for i = length(x):-1:1
                Hp(:,:,i) = [2*x(i)*obj.En+y(i)*obj.Em+2*x(i)^2*obj.Er+obj.Er*(x(i)^2+y(i)^2)+1 ...
                    x(i)*obj.Em+2*x(i)*y(i)*obj.Er;
                    obj.En*y(i)+2*x(i)*y(i)*obj.Er x(i)*obj.En+2*y(i)*obj.Em+2*obj.Er*y(i)^2+...
                    obj.Er*(x(i)^2+y(i)^2)+1];
                
            end
            for i = (1:length(x))*2
                H(i-1:i,:) = diag([obj.Kx obj.Ky])*Hp(:,:,i/2)*H(i-1:i,:);
            end
        end
        
        %%
        % *plotCameraFrame(obj,dt)* Simulate OpNav image by plotting
        % landmark positions on the camera CCD
        function plotCameraFrame(obj,dt)
            [s l] = getCameraFrame(obj,dt);
            plot(s,l,'*','MarkerSize',12);
            axis([0 obj.Rx 0 obj.Ry])
            set(gca,'ydir','reverse')
            set(gca,'XAxisLocation','top')
            xlabel('Pixel [s]')
            ylabel('Line [l]')
        end
    end
end
