function [] = linkbudget_default(linkbudget, fieldname, default_val)
    % LINKBUDGET_DEFAULT Sets a required linkbudget value to a given
    % default if the value doesn't already exist in the link budget
    % structure
    
    if ~isfield(linkbudget, fieldname) || isempty(linkbudget.(fieldname))
        linkbudget.(fieldname) = default_val;
    end
end