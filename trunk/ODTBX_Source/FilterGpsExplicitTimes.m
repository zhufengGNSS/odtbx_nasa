classdef FilterGpsExplicitTimes < FilterGpsData
    % A FilterGpsData class for filtering GPS measurement data by selecting
    % explicit time points.
    %
    % Any data for any PRN in the dataset will be preserved if the epoch is
    % the same as values provided when the class is created,
    % the filter_times property.  The time data uses the ODTBX convention
    % of datenum with convertTime GPS epoch (1xN double array).
    %
    % Note, frequently direct time comparisons of floating point values can
    % fail due to rounding or truncation error in calculations.  Therefore
    % a tolerance can be applied to the provided times during filtering.
    % The tolerance, if non-zero, is applied to preserve a time point
    % if (for each point) 
    %
    %   filter_times-TOL <= t_data <= filter_times+TOL,
    %       where t_data are the GPS time points being filtered
    %
    % The tolerance defaults to zero for exact matching where the time data
    % has not been manipulated and a direct bit match is acceptable.
    %
    % Usage:
    % Construct an object of this class with an array of times as the
    % argument, 
    %   "f = FilterGpsExplicitTimes(gps_times)".
    % If a tolerance should be applied to handle inexact matches, specify
    % that at construction:
    %   "f = FilterGpsExplicitTimes(gps_times, tol)".
    % 
    % Then call the filter
    % method of this object, "filtered_data = f.filter(input_data)".
    %
    % See also: makeGpsData, FilterGpsData, compValsTol
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
        filter_times = []; % the times to compare against
        tol = 0.0; % comparison tolerance, units are same as the timescale units
    end
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterGpsExplicitTimes(arg1, arg2)
            % Constructor that takes an array of times as the argument, and
            % an optional comparison tolerance.
            %
            % The time data uses the ODTBX convention
            % of datenum with convertTime GPS epoch (1xN double array).
            % The filter_times defaults to an empty array if not supplied.
            %
            % The tolerance defaults to zero if not specified.  Units are
            % the same timescale as the times (1.0 = 1 day in length).
            % Allowable values: zero or positive real values.
            %
            % Throws error if given an improperly-sized, non-empty, array,
            % or if the tolerance is negative.
            if nargin > 0
                if ~isempty(arg1) && size(arg1,1) ~= 1
                    error('FilterGpsExplicitTimes constructor given incorrectly sized argument (expected 1xN array)');
                end
                obj.filter_times = arg1;
                if nargin > 1
                    if length(arg2) ~= 1 || arg2 < 0
                        error('FilterGpsExplicitTimes constructor given improper tolerance');
                    end
                    obj.tol = arg2;
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
            if obj.tol == 0.0
                % no tolerance to apply, use exact bitmatches
                [~,~,ind] = intersect(obj.filter_times,gps_data.PRN_data{prn}.epoch);
            else
                % apply a tolerance for comparison
                [ind] = compValsTol(gps_data.PRN_data{prn}.epoch, obj.filter_times, obj.tol);
            end
        end
    end
    
end % classdef
