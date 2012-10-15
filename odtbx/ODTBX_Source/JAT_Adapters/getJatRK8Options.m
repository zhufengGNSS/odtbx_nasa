function o = getJatRK8Options(options, name, default)

% getJatRK8Options Get adaptor options for jatRK8.
%
%  This function operates on a different structure compared to 
%  GETODTBXOPTIONS and is applicable only to jatRK8.
%
%   VAL = getJatRK8Options(OPTIONS,'NAME') 
%
%   Extracts the value of the named property from integrator options structure 
%   OPTIONS, returning an empty matrix if the property value is not specified 
%   in OPTIONS. It is sufficient to type only the leading characters that 
%   uniquely identify the property. Case is ignored for property names. [] 
%   is a valid OPTIONS argument.
%
%   INPUTS
%   VARIABLE    SIZE    		DESCRIPTION (Optional/Default)
%      options     strucutre     	Input to function, Options structure
%                                   from setJatRK8Options
%      Name        string      	'', name of field in structure
%      default     (1x1)     		'', default value for the field that has not been set
%
%   OUTPUTS 
%      o    	   (1x1)			Output from function, Value of the field                                     
%
%   VAL = getJatRK8Options(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS. For example
%   returns val = 1.2 if the cr property is not specified in opts.
%
%   See setJatRK8Options for options details.
%
%
%    keyword: JAT Adaptor, options, RK8
%    See also jatRK8, setJatRK8Options
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
%   Author                Date          Comment
%                         (MM/DD/YYYY)
%   Mark W. Reichelt and
%   Lawrence F. Shampine  3/1/94        Copyright 1984-2003 The MathWorks, Inc.
%   									Revision: 1.37.4.3  Date: 2004/04/16 22:05:31
%   Kathryn Bradley      04/17/2006   	Ammended original to work for adaptor options
%   Emergent Space Technologies
% undocumented usage for fast access with no error checking
% (for additional changes, see the svn repository)
%   Allen Brown          02/11/2009     renamed to getJatRK8Options,
%                                       updated documentation

if (nargin == 4) && isequal(flag,'fast')
   o = getknownfield(options,name,default);
   return
end

if nargin < 2
  error('MATLAB:odeget:NotEnoughInputs','Not enough input arguments.');
end
if nargin < 3
  default = [];
end

if ~isempty(options) && ~isa(options,'struct')
  error('MATLAB:odeget:Arg1NotODESETstruct',...
        'First argument must be an options structure created with ODESET.');
end

if isempty(options)
  o = default;
  return;
end

 

Names = [
    'stepSize '
    'coefD    '
    'cr       '  
    'mass     '
    'cArea    '
    'mjd_utc  '
    'JGMOrder '
    'JGMDegree'];

names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error('MATLAB:odeget:InvalidPropName',...
        ['Unrecognized property name ''%s''.  ' ...
         'See ODESET for possibilities.'], name);
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    msg = sprintf('%s).', msg);
    error('MATLAB:odeget:AmbiguousPropName', msg);
  end
end

if any(strcmp(fieldnames(options),deblank(Names(j,:))))
  o = options.(deblank(Names(j,:)));
  if isempty(o)
    o = default;
  end
else
  o = default;
end

% --------------------------------------------------------------------------
function v = getknownfield(s, f, d)
%GETKNOWNFIELD  Get field f from struct s, or else yield default d.

if isfield(s,f)   % s could be empty.
  v = subsref(s, struct('type','.','subs',f));
  if isempty(v)
    v = d;
  end
else
  v = d;
end
