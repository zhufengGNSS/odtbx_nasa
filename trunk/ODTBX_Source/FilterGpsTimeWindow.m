classdef FilterGpsTimeWindow < FilterGpsData
    % A FilterGpsData class for filtering GPS measurement data by selecting
    % a window of time.
    %
    % Any data for any PRN in the dataset will be preserved if the epoch is
    % the same as, or in between, two epochs provided when the class is created,
    % the early & late properties.  The time data uses the ODTBX convention
    % of datenum with convertTime GPS epoch.
    %
    % Usage:
    % Construct an object of this class with two times as the argument,
    % "f = FilterGpsTimeWindow(time1,time2)".  Then call the filter
    % method of this object, "filtered_data = f.filter(input_data)".
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
        early = 0; % the earlier time of the window
        late = 0; % the later time of the window
    end
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterGpsTimeWindow(arg1, arg2)
            % Constructor that takes two times as the arguments.
            %
            % The time data uses the ODTBX convention
            % of datenum with convertTime GPS epoch.  The arguments are sorted
            % in time order and placed into the early and late properties.
            % The filter_times defaults to zero if not supplied.
            %
            % throws error if not provided zero or two arguments or if the
            % arguments are incorrectly specified
            if nargin > 0
                if nargin ~= 2
                    error('FilterGpsTimeWindow constructor given incorrect number of arguments (expects 2)');
                end
                if isempty(arg1) || length(arg1) ~= 1
                    error('FilterGpsTimeWindow constructor given incorrectly sized argument for arg 1');
                end
                if isempty(arg2) || length(arg2) ~= 1
                    error('FilterGpsTimeWindow constructor given incorrectly sized argument for arg 2');
                end
                if arg1 <= arg2
                    obj.early = arg1;
                    obj.late = arg2;
                else
                    obj.early = arg2;
                    obj.late = arg1;
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
            ind = find( gps_data.PRN_data{prn}.epoch >= obj.early & ...
                gps_data.PRN_data{prn}.epoch <= obj.late );
        end
    end
    
end % classdef
