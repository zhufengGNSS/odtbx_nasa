function y = getOdtbxOptions(options, name, default)

% GETODTBXOPTIONS Get Adaptor OPTIONS parameters.
%
%   y = getOdtbxOptions(OPTIONS, NAME, DEFAULT) 
%
%   Extracts the value of the NAME property from the structure 
%   OPTIONS, returning the DEFAULT value if the property value is not 
%   specified in OPTIONS. Case is ignored for property names.  If DEFAULT
%   is not provided [] is assumed.
%
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property.  However this use is discouraged because it
%   can make software maintenance more difficult.
%
%   INPUTS
%   VARIABLE    SIZE    		DESCRIPTION (Optional/Default)
%      options     structure    Options structure from ODTBXOPTIONS
%      name        string       name of field in structure
%      default     (1x1)     	(optional) default value for the field that
%                               has not been set, [] is an allowable
%                               default value for an argument
%
%   OUTPUTS 
%      y           (1x1)		Value of the specified field, or the given
%                               default argument
%
%
%    keyword: option, 
%    See also SETODTBXOPTIONS, ODTBXOPTIONS, VALIDATEODTBXOPTIONS, ODEGET
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
%   Author     			Date         	Comment
%               		(MM/DD/YYYY)
%   Derek Surka          07/20/2007     Original getJATOptions
%   Derek Surka          09/07/2007     Renamed and revised
%   Allen Brown          02/26/2009     Updated documentation
%   Ravi Mathur          06/29/2014     Uses standardized validation

if ((nargin < 2) || (nargin > 3))
  error('ODTBX:getOdtbxOptions:IncorrectInputs','Incorrect number of input arguments.');
elseif (nargin == 2)
  default = [];
end

% Fast validity check of input options
validateOdtbxOptions(options);

% Get the requested field's full name
[~, fullname] = getIndex(name,fieldnames(options),false); 

% Get the requested field value, default as needed
if isempty(fullname)
    y = default;
else
    y = options.(fullname);
    if isempty(y)
        y = default;
    end
end
