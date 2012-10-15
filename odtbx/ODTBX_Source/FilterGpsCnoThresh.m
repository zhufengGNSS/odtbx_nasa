classdef FilterGpsCnoThresh < FilterGpsData
    % A FilterGpsData class for filtering GPS measurement data via a C/No
    % threshold.
    %
    % Any data for any PRN in the dataset will be preserved if its C/No is
    % greater than or equal to the value provided when the class is created,
    % the cno_thresh property.
    %
    % Usage:
    % Construct an object of this class with a C/No threshold as the
    % argument, "f = FilterGpsCnoThresh(22)".  Then call the filter method of
    % this object, "filtered_data = f.filter(input_data)".
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
        cno_thresh = 0.0; % threshold value to apply
    end
    
    %% Methods
    methods (Access=public)
        
        % Constructor that takes a C/No threshold value as an argument.
        % The threshold defaults to zero if not supplied.
        function obj = FilterGpsCnoThresh(threshold)
            if nargin > 0
                obj.cno_thresh = threshold;
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
            ind = find(gps_data.PRN_data{prn}.raw_SNR >= obj.cno_thresh);
        end
    end
    
end % classdef
