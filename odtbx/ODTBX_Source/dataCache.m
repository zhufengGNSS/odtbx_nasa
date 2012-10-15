function [varargout] = dataCache(varargin)
%
% Creates and manipulates a data cache.
%
% dataCache is the function that creates, manipulates, and removes data
% from a data cache.  This function does not retain the data cache itself 
% so that multiple data caches can be managed.  Clients are expected to
% not directly manipulate the data cache produced by this function.
% This dataCache is similar to an unordered map, where all keys are unique.
% No restrictions are placed on the stored data other than the data should
% be able to be stored in a cell array element.  Adding data and getting
% data are designed to be efficient.  However deleting data is more
% heavyweight.
%
%   INPUTS:
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%   command     1 (string)  The command string for manipulating the cache:
%                           must be one of the following:
%                           'create', 'add', 'delete', 'get',
%                           'length', 'isempty'
%   cache       1 (struct)  The cache struct to be operated upon (not 
%                           required when command is 'create', reqired for
%                           all other commands)
%   key         1 (string)  The name/key identifier that is associated with
%                           data in the cache, zero-length key strings are
%                           not allowed (optional for some commands)
%   data        any         The data to be stored with the given key (only
%                           required for adding data)
%
%   OUTPUTS:
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%   cache       1 (struct)  The manipulated cache (optional)
%   OR
%   length      1 double    The number of separate data items stored in the
%                           cache, 0 if the cache is empty (optional)
%   OR
%   empty       1 double    Boolean indicating the cache is empty (1) or
%                           it contains data (0)
%   data        any         The data to returned with the given key, if no
%                           key matches or if the data stored is empty then
%                           [] is returned (optional)
%
% EXAMPLES:
% [cache] = dataCache('create');
%   Creates a cache struct for the client to hold.
%
% [cache] = dataCache('add', cache, key, data);
%   Adds data to the cache with the given key.  If the key already exists
%   in the cache the data is overwritten.
%
% [cache] = dataCache('delete', cache, key);
%   Removes data from the cache using the given key.  If the key does not
%   exist in the cache then the cache is not altered.
%
% [data] = dataCache('get', cache, key);
%   Returns data from the cache using the given key.  If the key does not
%   exist in the cache then the returned data is []. (Does not alter cache)
%
% [length] = dataCache('length', cache);
%   Returns the length of the cache, i.e. the number of keys currently in
%   the map.  Note, this does not check if any data is associated with the
%   keys.  Returns 0 if the cache is empty. (Does not alter cache)
%
% [empty] = dataCache('isempty', cache);
%   Returns 1 if the cache is empty or 0 if the cache has some data. (Does
%   not alter cache)
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

% Implementation:
% The data cache itself is mplemented as a struct 
% with two cell arrays, one that holds the identifying field for the
% cached data and the other that holds the data itself.  
% cache.keys and cache.data, both nx1,
% The sizes and layout of cache.data and cache.keys are always kept
% consistent.  However, they are be sized 0x1 at creation.

if nargin < 1 || ~ischar(varargin{1})
    error('dataCache.m: Client did not specify correct input arguments.\nHere is the help output:\n%s', help('dataCache'));
end
if nargout < 1
    error('dataCache.m: Client did not specify correct output arguments.\nHere is the help output:\n%s', help('dataCache'));
end

varargout = cell(1,1);

% run through the commands in order of most likely first:
switch lower(varargin{1})
    case 'get',
        if nargin ~= 3 || nargout ~= 1 || ~isstruct(varargin{2}) || ...
                ~ischar(varargin{3}) || ~isfield(varargin{2},'keys') || ...
                (length(varargin{3}) < 1)
            error('dataCache.m: Bad inputs/outputs for get command.');
        end
                
        found = 0;
        
        if ~isempty(varargin{2})
            % search for the key
            for i = 1:length(varargin{2}.keys)
                if strcmp(varargin{3},varargin{2}.keys{i})
                    found = i;
                    break;
                end
            end
        end
        
        if found == 0
            varargout{1} = [];
        else
            varargout{1} = varargin{2}.data{i};
        end            
        
    case 'add',
        if nargin ~= 4 || nargout ~= 1 || ~isstruct(varargin{2}) || ...
                ~ischar(varargin{3}) || ~isfield(varargin{2},'keys') || ...
                (length(varargin{3}) < 1)
            error('dataCache.m: Bad inputs/outputs for add command.');
        end
        
        overwrite = 0;
        
        if ~isempty(varargin{2})
            % search for the key
            for i = 1:length(varargin{2}.keys)
                if strcmp(varargin{3},varargin{2}.keys{i})
                    overwrite = i;
                    break;
                end
            end
        end
        
        varargout{1} = varargin{2};
        if overwrite > 0
            % overwrite at this index:
            varargout{1}.data{overwrite,1} = varargin{4};
        else
            if isempty(varargin{2})
                % ah, the first one to add
                varargout{1}(1).keys = cell(1,1);
                varargout{1}.keys{1,1} = varargin{3};
                varargout{1}(1).data = cell(1,1);
                varargout{1}.data{1,1} = varargin{4};
            else
                % don't overwrite, append to end:
                varargout{1}.keys{end+1,1} = varargin{3};
                varargout{1}.data{end+1,1} = varargin{4};
            end
        end
        
    case 'length',
        if nargin ~= 2 || nargout ~= 1 || ~isstruct(varargin{2}) || ...
                ~isfield(varargin{2},'keys')
            error('dataCache.m: Bad inputs/outputs for length command.');
        end
        
        if isempty(varargin{2})
            varargout{1} = 0;
        else
            varargout{1} = length(varargin{2}.keys);
        end
        
    case 'isempty',
        if nargin ~= 2 || nargout ~= 1|| ~isstruct(varargin{2}) || ...
                ~isfield(varargin{2},'keys')
            error('dataCache.m: Bad inputs/outputs for isempty command.');
        end
        
        varargout{1} = isempty(varargin{2}) || isempty(varargin{2}.keys);
        
    case 'create',
        if nargin ~= 1 || nargout ~= 1
            error('dataCache.m: Bad inputs/outputs for create command.');
        end
        varargout{1} = struct('keys',cell(0,0),'data',cell(0,0));

    case 'delete',
        if nargin ~= 3 || nargout ~= 1 || ~isstruct(varargin{2}) || ...
                ~ischar(varargin{3}) || ~isfield(varargin{2},'keys') || ...
                (length(varargin{3}) < 1)
            error('dataCache.m: Bad inputs/outputs for add command.');
        end
        
        found = 0;
        if ~isempty(varargin{2})
            % search for the key
            for i = 1:length(varargin{2}.keys)
                if strcmp(varargin{3},varargin{2}.keys{i})
                    found = i;
                    break;
                end
            end
        end

        if found > 0
            % We found data to be deleted
            newlength = length(varargin{2}.keys)-1;
            varargout{1}(1).keys = cell(newlength,1);
            varargout{1}(1).data = cell(newlength,1);
            
            if found > 1
                % all items below found are the same:
                for i = 1:found-1
                    varargout{1}.keys{i} = varargin{2}.keys{i};
                    varargout{1}.data{i} = varargin{2}.data{i};
                end
            end
            if found < length(varargin{2}.keys)
                % Shift any the higher data down by one
                for i = found:newlength
                    varargout{1}.keys{i} = varargin{2}.keys{i+1};
                    varargout{1}.data{i} = varargin{2}.data{i+1};
                end
            end
            
        else
            varargout{1} = varargin{2}; % not found, no change
        end
        
    otherwise,
        error('dataCache.m: Client did not specify correct input arguments.\nHere is the help output:\n%s', help('dataCache'));
end

