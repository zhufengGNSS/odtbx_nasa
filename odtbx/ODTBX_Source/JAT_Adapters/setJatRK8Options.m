function options = setJatRK8Options(varargin)

% SETJATRK8OPTIONS Set adaptor options for jatRK8.
%
%  This function operates on a different structure compared to 
%  SETODTBXOPTIONS and is applicable only to jatRK8.
%
%   OPTIONS = setJatRK8Options('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   OPTIONS = setJatRK8Options(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = setJatRK8Options(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%   
%   setJatRK8Options with no input arguments displays all property names and their
%   possible values.
%
%   INPUTS
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%   stepSize    (1x1)		Input to function, This is the stepSize for the integrator in seconds.
%   coefD       (1x1)		 '', This is the coefficient of drag for the spacecraft in orbit.
%   cr          (1x1)		 '', Coefficient of reflectivity for the spacecraft in orbit.
%   mass        (1x1)		 '', The mass of the spacecraft in orbit in kilograms.
%   cArea       (1x1)		 '', The cross-sectional area of the spacecraft in orbit in m^2.
%   mjd_utc     (1x1)		 '', The Modified Julian Date in UTC time.
%   JGMOrder    (1X1)        '', Order for the JGM Universe
%   JGMDegree   (1x1)        '', Degree for the JGM Universe
%
%   OUTPUTS 
%   options     structure    Output from function, Options structure                                    
%
%   NOTE: 
%     Some of the properties available through setJatRK8Options have changed in MATLAB 6.0. 
%     Although we still support the v5 properties when used with the v5 syntax 
%     of the ODE solvers, any new functionality will be available only with the 
%     new syntax. To see the properties available in v5, type in the command line  
%     more on, type setJatRK8Options, more off
%
%   NOTE:
%     This portion describes the properties available in v5
%
%    keyword: JAT Adaptor, options, RK8
%    See also jatRK8, getJatRK8Options
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
%   Mark W. Reichelt and Lawrence F. Shampine 	5/6/94		Copyright 1984-2003 The MathWorks, Inc.
%   											Revision: 1.37.4.3   Date: 2004/04/16 22:05:31
%   Kathryn Bradley      04/17/2006     Ammended original to work for 
%   Emergent Space Technologies         adaptor options
% (for additional changes, see the svn repository)
%   Allen Brown          02/11/2009     renamed to getJatRK8Options,
%                                       updated documentation

Names = [
    'stepSize ' 
    'coefD    '
    'cr       '  
    'mass     '
    'cArea    '
    'mjd_utc  '
    'JGMOrder '
    'JGMDegree'];

m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in setJatRK8Options(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error('MATLAB:setJatRK8Options:NoPropNameOrStruct',...
            ['Expected argument %d to be a string property name ' ...
                     'or an options structure\ncreated with setJatRK8Options.'], i);
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('MATLAB:setJatRK8Options:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error('MATLAB:setJatRK8Options:NoPropName',...
            'Expected argument %d to be a string property name.', i);
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error('MATLAB:setJatRK8Options:InvalidPropName',...
            'Unrecognized property name ''%s''.', arg);
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error('MATLAB:setJatRK8Options:AmbiguousPropName', msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error('MATLAB:setJatRK8Options:NoValueForProp',...
        'Expected value for property ''%s''.', arg);
end
