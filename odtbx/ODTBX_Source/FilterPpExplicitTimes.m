classdef FilterPpExplicitTimes < FilterPpData
    % A FilterPpData class for filtering GPS physical parameter by time.
    %
    % This filter will preserve any data in the dataset that match the epoch
    % times provided when the class is created, the filter_times property.
    %
    % The time data uses the GPS physical parameter convention, a MATLAB 
    % serial date in UTC (1xN double array).
    %
    % Note, this comparison uses an exact bitmatch with no tolerance
    % between the given times and the GPS physical parameter times.
    %
    % Usage:
    % Construct an object of this class with a given set of times, as in:
    %
    %   f = FilterPpExplicitTimes(time_data);
    %
    % Then call the filter method of
    % this object, "filtered_data = f.filter(input_data)".
    %
    % See also: FilterPpData, gps_phys_params
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
        filter_times = []; % MATLAB serial date, UTC
    end
    
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterPpExplicitTimes(arg)
            % Constructor that takes an array of times as the argument.
            %
            % The time data uses the GPS physical parameter convention, a 
            % MATLAB serial date in UTC (1xN double array).
            % The filter_times defaults to an empty array if not supplied.
            %
            % throws error if given an improperly-sized, non-empty, array.
            if nargin > 0
                if ~isempty(arg) && size(arg,1) ~= 1
                    error('FilterPpExplicitTimes constructor given incorrectly sized argument (expected 1xN array)');
                end
                obj.filter_times = arg;
            end
        end
        
        function ind = getInd(obj, pp_data)
            % Function to filter the data in the GPS physical parameter dataset.
            %
            % Outputs: 1xN set of selected indices to retain from a filter
            %       instance from the pp_data
            % Inputs:
            %   obj: this instance
            %   pp_data: (single) struct of GPS physical parameter data
            %   populated by gps_phys_params.
            %
            % Throws an error if any argument is invalid or incorrectly
            % constructed.
            
            if ~isstruct(pp_data) || ~isfield(pp_data,'epoch')
                error('filter given arg that is not a GPS physical parameters struct');
            end
            
            [~,~,ind] = intersect(obj.filter_times,pp_data.epoch);
            
        end % function
        
    end % methods
    
end % classdef