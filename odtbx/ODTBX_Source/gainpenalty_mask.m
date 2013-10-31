function [link_budget] = gainpenalty_mask(antenna, link_budget, alpha, beta) 
% Check and apply additional mask angle penalty to:
% Ar, RP, CN0
% Note, the linkbudget already applies a penalty for those 
% angles outside the pattern so don't doubly-apply a penalty.

  % Apply the receiver gain penalty for the user-defined mask angle
    if beta < (max(antenna.pattern(:,1)))*pi/180
        % Check and apply additional mask angle penalty to:
        % Ar, RP, CN0
        % Note, the gpslinkbudget already applies a penalty for those 
        % angles outside the pattern so don't doubly-apply a penalty.
        mask_ind = (alpha > beta) & (link_budget.Ar ~= -100.0);
        if sum(sum((mask_ind))) > 0
            link_budget.Ar(mask_ind) = -100; % change Ar
            % alter dependent values:
            link_budget.RP(mask_ind) = link_budget.RP(mask_ind) - 100;
            link_budget.CN0(mask_ind) = link_budget.CN0(mask_ind) - 100;
        end
    end

end