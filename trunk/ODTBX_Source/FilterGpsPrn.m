classdef FilterGpsPrn < FilterGpsData
    % A FilterGpsData class for filtering GPS measurement data by PRN
    % values.
    %
    % This filter preserves all data for a PRN in the dataset by comparing
    % it to a list of given PRNs at constrution.  This class can either
    % include or exclude data compared to the set PRNs.
    %
    % Usage:
    % Construct an object of this class with an array of PRNs and whether
    % to treat them as PRNs to be "included" or "excluded", examples:
    % ex: To include data only from any of the given PRNs:
    %     "f = FilterGpsPrn('include',[1 2 4 13 21])".
    % ex: To exclude data from any of the given PRNs (and include all others):
    %     "f = FilterGpsPrn('exclude',[1 2 4 13 21])".
    % Then call the filter method of this object,
    % "filtered_data = f.filter(input_data)".
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
        PRN_list = []; % the list of PRNs to preserve
        
        % true=include data from only PRNs that match the PRN_list,
        % false=include data from any PRN other than those on the PRN_list
        include = true;
    end
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterGpsPrn(flag, PRNs)
            % Constructor that takes an include/exclude flag and a list of PRNs
            % as arguments.
            %
            % flag is either 'include' or 'exclude', no other value is allowed,
            % defaults to 'include' if empty.
            %
            % The PRNs argument is a 1xN list of PRN values.
            % The PRN_list property defaults to [] if not supplied.
            %
            % throws error if not provided zero or two arguments or if the
            % arguments are incorrectly specified
            
            if nargin > 0
                if nargin ~= 2
                    error('FilterGpsPrn constructor given incorrect number of arguments (expects 2)');
                end
                flagi = strcmpi(flag,'include');
                flage = strcmpi(flag,'exclude');
                if ~flagi && ~flage
                    error('FilterGpsPrn constructor given incorrect flag, expects ''include'' or ''exclude'' (given ''%s'')',flag);
                end
                if isempty(PRNs) || size(PRNs,1) ~= 1
                    error('FilterGpsPrn constructor given incorrectly sized argument for PRNs');
                end
                obj.include = flagi;
                obj.PRN_list = sort(PRNs);
                if any(diff(obj.PRN_list) == 0)
                    error('FilterGpsPrn constructor given duplicate PRNs.');
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
            
            ind = [];
            if obj.isSelected(gps_data, prn)
                ind = 1:length(gps_data.PRN_data{prn}.epoch);
            end
        end
    end % methods
    
    
    %% Methods
    methods (Access=protected)
        function sel = isSelected(obj, gps_data, prn)
            % If the specified prn in the gps_data is "selected" for inclusion
            %
            % False if the prn should not be included in the outptut 
            % data.
            %
            % Inputs:
            %   obj: this instance
            %   gps_data: (single) struct of GPS data defined by makeGpsData, it
            %       can contain zero, one, or many sets of GPS PRN data.
            %   prn: (scalar) index into the gps_data.PRN_data to filter
            %
            % Throws an error if any argument is invalid or incorrectly
            % constructed.
            
            sel = false;
            
            if obj.include
                if ~isempty(obj.PRN_list) && any(obj.PRN_list == gps_data.GPS_PRN(prn))
                    sel = true;
                end
            else
                if isempty(obj.PRN_list) || ~any(obj.PRN_list == gps_data.GPS_PRN(prn))
                    sel = true;
                end
            end
        end
        
    end % methods
    
end % classdef
