function [linkbudget] = linkbudget_default(linkbudget, fieldname, default_val)
    % LINKBUDGET_DEFAULT Sets a required linkbudget value to a default if
    % the value doesn't already exist in the link budget structure. Warns
    % user when a value is set as a default.
    
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
        if iscell(default_val)
            defaultname = mat2str(cell2mat(default_val));
        else
            defaultname = num2str(default_val);
        end
        warning('%s not set, setting default value of %s.', fieldname, defaultname);
    end
end