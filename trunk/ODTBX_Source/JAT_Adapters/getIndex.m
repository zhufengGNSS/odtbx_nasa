function [index,fullname] = getIndex(str,strArray,requireMatch,errmsg)

% GETINDEX Returns the index within strArray of the entry that starts with str.
%
%   [index,fullname] = getIndex(str,strarray,requireMatch,errmsg)
%
%   Returns the index within strArray of the entry that starts with
%   str. It is sufficient to type only the leading characters that
%   uniquely identify the property. Case is ignored for property names. If
%   a second output argument is specified, the full name of the matching
%   string will be returned.
%
%   If requireMatch is true (default), an error message will be displayed
%   if a match is not found. The error message can be specified in errmsg,
%   otherwise, a default message will be used.
%
%   INPUTS
%   VARIABLE            SIZE    		DESCRIPTION (Optional/Default)
%      str             string     	Input string to be matched
%      strArray         (1xN)       Cell array of possible matching strings
%      requireMatch     (1x1)     	true = error if no match found
%                                   false = return empty index but no error
%      errmsg          string       Optional error message to be displayed if a match is not found 
%
%   OUTPUTS
%      index            (1x1)		Index into strArray of matching string
%      fullname        string       Full matching string 
%
%
%    keyword: options
%    See also ODTBXOPTIONS, GETODTBXOPTIONS, SETODTBXOPTIONS
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
%   Author      	Date         	Comment
%               	MM/DD/YYYY
%   Derek Surka     06/29/2007     Generalized from Mathworks' getOptions
   
% JAG 6/23/10 moved creation of errmsg so that it is not created unless
% needed.

if( nargin < 3 )
    requireMatch = true;
end

j=find(strcmpi(str,strArray));
% This should compare the whole entries, so the whole name must be entered
% There should not be more than one output.

% IF there is a partial match for an option name the code will be slow.
if isempty(j)
    lowStr      = lower(str);
    lowStrArray = lower(strArray);
    j           = strmatch(lowStr,lowStrArray);
end

if isempty(j)               % if no matches
    if( requireMatch )
        if( nargin < 4 )
            errmsg = ['MATLAB:getIndex: Unrecognized property ' str];
        end
        error(errmsg);
    else
        fullname = [];
    end
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = strmatch(lowStr,lowStrArray,'exact');
    if length(k) == 1
        j        = k;
        fullname = strArray{j};
    else
        msg = ['Ambiguous property: ' str '. Possible matches {' strArray{j(1)}];
        for k = j(2:length(j))
            msg = [msg ', ' strArray{k}];
        end
        msg = [msg '}'];
        error(['MATLAB:getIndex: ', msg]);
    end
else
    fullname = strArray{j};
end

index = j;
