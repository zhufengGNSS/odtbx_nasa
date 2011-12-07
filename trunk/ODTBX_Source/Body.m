classdef Body < handle
    %BODY Class representing a simple solar system body (i.e. planet, asteroid, comet, etc.)
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
    %   Kenneth Getzandanner    03/29/2011      Original Body.m
    
    %% Properties
    properties
        a;      % semi-axis (a)
        b;      % semi-axis (b)
        c;      % semi-axis (c)
        RA;     % Right Ascension (Spin Axis)
        DEC;    % Declination (Spin Axis)
        PRA;    % Prime Meridian Angle at Epoch
        w;      % Spin Rate
        lmk;    % Landmark Locations [nx2]
        SPID;   % SPICE ID
        GM;     % Gravitational Parameter
        epoch;  % Epoch    
    end
    
    %% Methods
    methods
        %%
        % *Body(varargin)* [constructor]
        function obj = Body(SA,SPIN,LMK,EPOCH,SPID,GM)
            obj.a = SA(1);
            obj.b = SA(2);
            obj.c = SA(3);
            obj.RA = SPIN(1);
            obj.DEC = SPIN(2);
            obj.PRA = SPIN(3);
            obj.w = SPIN(4);
            obj.lmk = LMK;
            obj.epoch = EPOCH;
            obj.SPID = SPID;
            obj.GM = GM;
        end
        
        %%
        % *getBodyLmk(obj)* Get the cartesian position of each landmark in
        % the body fixed frame
        function [lmk] = getBodyLmk(obj)
            
            % Loop through each landmark
            for i=size(obj.lmk,1):-1:1
                
                % Get latitude and longitude of landmark
                B = obj.lmk(i,1);
                L = obj.lmk(i,2);
                
                % Convert from lat/long to body fixed cartesian (XYZ) 
                x = obj.a*cos(B)*cos(L);
                y = obj.b*cos(B)*sin(L);
                z = obj.c*sin(B);
                
                % Return cartesian vector
                lmk(:,i) = [x;y;z];
            end
        end
        
        %%
        % *getBodyLmkNormal(obj)* Get the landmark normal vector in the
        % body fixed frame
        function [n] = getBodyLmkNormal(obj)
            
            % Get body-fixed landmark locations
            xb = getBodyLmk(obj);
            
            % Get body-fixed landmark normal vectors
            for i=size(obj.lmk,1):-1:1
                n(:,i) = unit([2*xb(1,i)/obj.a^2;...
                    2*xb(2,i)/obj.b^2; ...
                    2*xb(3,i)/obj.c^2]);
            end
        end
        
        %%
        % *getInertialLmk(obj,dt)* Get the cartesian position of each
        % landmark at the inertial frame dt seconds from epoch
        function xi = getInertialLmk(obj,dt)
            
            % Calculate pole rotation at time dt
            theta = obj.PRA + obj.w*dt;
            
            % Generate DCM to rotate from body-fixed to inertial
            % coordinates using the body spin state
            D3w = [cos(-theta) sin(-theta) 0;
                -sin(-theta) cos(-theta) 0;
                0 0 1];
            
            D1 = [1 0 0;
                0 cos(obj.DEC-pi/2) sin(obj.DEC-pi/2);
                0 -sin(obj.DEC-pi/2) cos(obj.DEC-pi/2)];
            
            D3 = [cos(-pi/2-obj.RA) sin(-pi/2-obj.RA) 0;
                -sin(-pi/2-obj.RA) cos(-pi/2-obj.RA) 0;
                0 0 1];
            
            D = D3*D1*D3w;
            
            xb = getBodyLmk(obj);
            
            % Rotate landmark coordinates into inertial frame using the DCM
            for i=size(obj.lmk,1):-1:1
                xi(:,i) = D*xb(:,i);
            end
        end
        
        %%
        % *getInertialLmkNormal(obj,dt)* Get the landmark normal vector in 
        % the inertial frame
        function ni = getInertialLmkNormal(obj,dt)
            
            % Calculate pole rotation at time dt
            theta = obj.PRA + obj.w*dt;
            
            % Generate DCM to rotate from body-fixed to inertial
            % coordinates using the body spin state
            D3 = dcm('ax3',obj.RA);
            D2 = dcm('ax2',obj.DEC);
            D3w = dcm('ax3',theta);
            D = D3w*D2*D3;
            
            % Get landmark normal vectors in the body frame
            nb = getBodyLmkNormal(obj);
            
            % Rotate landmark normals into inertial frame using the DCM
            for i=size(obj.lmk,1):-1:1
                ni(:,i) = D'*nb(:,i);
            end
        end
        
        %%
        % *plotBody(obj)* Generates a 3D model of the body using *surf*
        function plotBody(obj)
            
            % Generate and plot an ellipsoid based on body dimensions
            [x y z] = ellipsoid(0,0,0,obj.a,obj.b,obj.c);
            surf(x,y,z);
            axis equal
            
            % Change body color to gray
            kd = get(gcf,'children');
            kd = get(kd,'children');
            
            sp = kd(1);
            set(sp,'EdgeColor',[0.5 0.5 0.5])
            set(sp,'FaceColor',[229 229 229]/255)
        end
        
        %%
        % *plotLmk(obj)* Generates a 3D model of the body with landmarks
        function plotLmk(obj)
            
            % Get body fixed landmark locations and normals
            xb = getBodyLmk(obj);
            nb = getBodyLmkNormal(obj);
            
            % Calculate scale based on mean radius
            mr = (obj.a+obj.b+obj.c)/3;
            s = 0.05;
            
            % Plot the body and add landmarks (red circles)
            hold on
            plotBody(obj)
            for i=1:size(xb,2)
                h = obj.plotCircle3D(xb(:,i)',nb(:,i)',mr*s);
                set(h,'LineWidth',2)
            end
            hold off
            
        end
        
        %%
        % *getState(obj,dt,FRAME,CB)* Returns the inertial position of the
        % body at dt seconds past the epoch in the reference frame FRAME
        % with respect to the central body CB.
        function X = getState(obj,dt,FRAME,CB)
            et = obj.epoch + dt;
            X = cspice_spkezr(obj.SPID, et,FRAME,'NONE',CB);
        end
    end
    
    %% Methods (Static)
    methods (Static)
        
        %%
        % *plotCircle3D(center,normal,radius)* Plot a circle in 3D given a
        % center location, radius, and normal vector
        function h = plotCircle3D(center,normal,radius)
            theta=0:0.01:2*pi;
            v=null(normal);
            points=repmat(center',1,size(theta,2))+...
                radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
            h = plot3(points(1,:),points(2,:),points(3,:),'r-');
            
        end
    end
end

