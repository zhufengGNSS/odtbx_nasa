function [valid] = validateOdtbxOptions(options)

% VALIDATEODTBXOPTIONS  Checks to see if the given structure is valid.
%
%   VALID = VALIDATEODTBXOPTIONS(OPTIONS) returns TRUE if the input options
%   structure is a valid ODTBX options structure. Extra field names in the
%   input options structure are ignored.
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

if ~isa(options,'struct') || ~isfield(options,'optionsType')
    valid = false;
else
    validStruct = odtbxOptions(options.optionsType);
    fields = fieldnames(validStruct);

    validFields = isfield(options,fields);
    if ~all(validFields)
        valid=false;
    else
        valid=true;
    end
end
