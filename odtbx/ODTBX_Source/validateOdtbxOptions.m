function [valid, options] = validateOdtbxOptions(options)

% VALIDATEODTBXOPTIONS  Checks to see if the given structure is valid.
%
%   VALID = VALIDATEODTBXOPTIONS(OPTIONS) returns TRUE if the input options
%   structure is a valid ODTBX options structure. Extra field names in the
%   input options structure are ignored.
%
%   [VALID, OPTIONS] = VALIDATEODTBXOPTIONS(OPTIONS) returns an updated
%   options structure containing fields missing from the input structure.
%   New fields have default values as given by ODTBXOPTIONS.
%
%   keyword: options
%   See also ODTBXOPTIONS, GETODTBXOPTIONS, SETODTBXOPTIONS
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
%   Author         Date            Comment
%                 (MM/DD/YYYY)
%   Derek Surka    07/20/2007      Original validateJATOptions
%   Derek Surka    09/07/2007      Renamed and revised
%   Allen Brown    02/26/2009      Removed undocumented usage in favor
%                                  of getOdtbxOptions.
%   Ravi Mathur    06/29/2014      Add missing fields if desired

% First do a quick sanity check on the input
if ~isa(options,'struct') || ~isfield(options,'optionsType')
    valid = false;

% Now do a full check on input compared to a valid OdtbxOptions structure
else
    validOptions = odtbxOptions(options.optionsType); % Valid OdtbxOptions struct of same type
    validFields = fieldnames(validOptions);

    % All fields from valid OdtbxOptions struct must be in input struct
    hasFields = isfield(options,validFields);
    valid = all(hasFields);
    
    % Return immediately if input is valid or shouldn't be fixed
    if(valid || (nargout < 2))
        return;
    end
    
    % Add default value for any missing field in input struct
    missingInd = find(~hasFields);
    for i = 1:length(missingInd)
        currField = validFields{missingInd(i)};
        options.(currField) = validOptions.(currField);
    end
    
    % Display number of fields added to the input struct
    warning('ODTBX:validateOdtbxOptions:AddedOptions', 'added %d fields to input struct', length(missingInd));
    valid = true; % Input is now valid
end
