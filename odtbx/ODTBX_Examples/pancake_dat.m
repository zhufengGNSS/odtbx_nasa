function [y H R] = pancake_dat(t,x,~)
% Get measurement data for the pancake tutorial.

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

% Find length of timespan
lent = length(t);

% Initialize output variables
y = nan(1,lent);
H = nan(1,7,lent);
R = ones(1,1,lent);

% Calculate measurments and partials at each time step
for i=length(t):-1:1
    
    % Spacecraft/Ground Station position
    Rs = x(1:2,i);
    Rg = x(6:7,i);
    
    % Check visibility
    if dot(Rg,Rs-Rg)>=0
        
        % Calculate range
        y(1,i) = norm(Rs-Rg);
        
        % Calculate partials
        %H(:,:,i) = [unit(Rs-Rg)' zeros(1,5)];
        H(:,:,i) = [unit(Rs-Rg)' zeros(1,3) -unit(Rs-Rg)'];
    end

end