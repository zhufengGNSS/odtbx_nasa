classdef FilterGpsMulti < FilterGpsData
    % A FilterGpsData class that uses multiple FilterGpsData filters on GPS measurement data.
    %
    % This filter stores and applies multiple user-supplied FilterGpsData
    % filters to the same GPS measurement data at the same time as they
    % were a single filter.  Data is preserved if all sub-filters agree.
    % A FilterGpsMulti with no added filters does not filter out any data.
    %
    % Note, like all FilterGpsData classes, this is a MATLAB Value Class.
    % It is designed so that you can create copies of it and manipulate
    % each separately.  However with value classes some functions, such as
    % add(), create new objects rather than updating the current object.
    %
    % Usage:
    % Construct an object of this class and add several other filters to it,
    % example: Construct and add multiple filters:
    %
    %     "f1 = FilterCnoThresh(22);"
    %     "f2 = FilterGpsMulti('include',[1 2 4 13 21]);"
    %     "f3 = FilterTimeWindow(time1,time2); % (times previously defined)"
    %
    %     "f = FilterGpsMulti(f1, f2, f3); % added during construction"
    %       or
    %     "f = FilterGpsMulti; % added after construction"
    %     "f = f.add(f1);" % returns a different object, but re-assigns to f
    %     "f = f.add(f3);"
    %     "f = f.add(f2);"
    %
    % Then call the filter method of this object,
    %     "filtered_data = f.filter(input_data)".
    %
    % See also: makeGpsData, FilterGpsData
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
    
    %% Public Properties
    properties (GetAccess=public, SetAccess=protected)
        filters = cell(1,0); % a cell array of FilterGpsData filters
    end
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterGpsMulti(varargin)
            % Constructor that can optionally take other filters to use.
            %
            % One or more FilterGpsData objecs can be added as a variable
            % number of arguments.  An empty argument list means that no
            % filters will be applied.
            %
            % Throws error if provided an argument that is not a FilterGpsData
            % object.
            if nargin > 0
                obj = obj.add(varargin{:});
            end
        end
        
        function obj = add(obj, varargin)
            % Function to add one or more filters to use.
            %
            % Note that with value classes this function creates new objects
            % rather than updating the current object.  To add a filter to an
            % existing set, use:
            % "obj = obj.add(...)"
            %
            % One or more FilterGpsData objecs can be added as a variable
            % number of arguments.
            %
            % Throws error if provided an argument that is not a FilterGpsData
            % object.
            
            if nargin > 1 % arg 1 is always the obj
                for i = 2:nargin
                    if isa(varargin{i-1},'FilterGpsData')
                        obj.filters{end+1} = varargin{i-1};
                    else
                        error('FilterGpsMulti passed an argument that is not a FilterGpsData type.');
                    end
                end
            end
        end
        
        function ind = getInd(obj, gps_data, prn)
            % Filter the data of a specific PRN in the GPS dataset.
            %
            % Outputs: 1xN set of selected indices to retain from a filter
            %       instance from the gps_data and prn selection.
            %
            % Inputs:
            %   obj: this instance
            %   gps_data: (single) struct of GPS data defined by makeGpsData, it
            %       can contain zero, one, or many sets of GPS PRN data.
            %   prn: (scalar) index into the gps_data.PRN_data to filter
            %
            % Throws an error if any argument is invalid or incorrectly
            % constructed.
            if isempty(gps_data)
                error('getInd passed empty gps_data argument');
            end
            if isempty(prn)
                error('getInd passed empty prn argument');
            end
            if (prn > length(gps_data.GPS_PRN)) || (prn > length(gps_data.PRN_data))
                error('getInd passed out-of-range PRN index (passed: %d, max: %d)', length(gps_data.PRN_data));
            end
            
            ind = []; % default to all filtered out
            
            if isempty(obj.filters)
                % no filters means preserve all data
                ind = 1:length(gps_data.PRN_data{prn}.epoch);
                return;
            end
            
            res = cell(size(obj.filters)); % to hold sub-filter results
            
            % run through and apply each filter
            for i = 1:length(obj.filters)
                res{i} = obj.filters{i}.getInd(gps_data, prn);
                
                % bail early if any filter removes all data
                if isempty(res{i})
                    return; % empty ind
                end
            end
            
            % only one filter?, we're done
            if length(res) == 1
                ind = res{1};
                return;
            end
            
            % Since intersection usually downsizes the population,
            % intersect the data in order from fewest returned indices to
            % greatest number of returned indices for efficiency.
            stats = zeros(1,length(res)); % # ind returned from each filter
            for i = 1:length(res)
                stats(i) = length(res{i});
            end
            [~, sind] = sort(stats); % default ascending sort
            ind = res{sind(1)};
            for i = 2:length(res)
                ind = intersect(ind,res{sind(i)});
            end
            
        end
    end % methods
    
end % classdef
