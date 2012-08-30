function [y H R] = pancake_dat(t,x,~)

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
        H(:,:,i) = [unit(Rs-Rg)' zeros(1,5)];
    end

end