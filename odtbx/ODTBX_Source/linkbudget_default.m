function [linkbudget] = linkbudget_default(linkbudget, fieldname, default_val)
    % LINKBUDGET_DEFAULT Sets a required linkbudget value to a given
    % default if the value doesn't already exist in the link budget
    % structure
    
    if nargin < 2
        error('MATLAB:linkbudget_default:NotEnoughInputs','Not enough input arguments.');
    end
    if nargin < 3
        default = [];
    end

    if ~isempty(linkbudget) && ~isa(linkbudget,'struct')
        error('MATLAB:linkbudget_default:Arg1Notlinkbudgetstruct',...
            'First argument must be a link budget structure.');
    end
    
    if ~isfield(linkbudget, fieldname) || isempty(linkbudget.(fieldname))
        linkbudget.(fieldname) = default_val;
        warning('%s not set, setting default value of %s.', fieldname, num2str(default_val));
    end
end