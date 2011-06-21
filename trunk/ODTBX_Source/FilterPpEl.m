classdef FilterPpEl < FilterPpData
    % A FilterPpData class for filtering GPS physical parameter data on antenna
    % boresight (elevation) angle.
    %
    % This filter will preserve any data in the dataset when the elevation
    % angle is less than or equal to the value provided when the class
    % is created, the ang_thresh property.  Either the transmit elevation
    % angle or the receive elevation angle can be specified at construction,
    % the angtype property.
    %
    % Usage:
    % Construct an object of this class with an angle threshold and whether
    % to apply it to the transmit ('tx') or receive ('rx') antenna, as in:
    %
    %   % only preserve transmit elevation angles of pi/2 rad or higher:
    %   f = FilterPpEl(pi/2, 'tx')
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
        ang_thresh = 0.0; % angle threshold for comparison (rad)
        angtype = 1; % apply threshold to: 1=transmit, 2=receive antenna angles
    end
    
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterPpEl(ang, type)
            % Constructor that takes an angle threshold value and antenna designator.
            %
            % ang: elevation angle threshold (rad), defaults to 0.0 rad
            % type: antenna designation: transmit ('tx') or receive ('rx')
            %   antenna, no other value is allowed
            %
            % Throws exception if the type is incorrect.
            if nargin > 0
                obj.ang_thresh = ang;
            end
            
            if nargin > 1
                flagr = strcmpi(type,'rx');
                flagt = strcmpi(type,'tx');
                if ~flagr && ~flagt
                    error('FilterPpEl constructor given incorrect type, expects ''rx'' or ''tx'' (given ''%s'')',type);
                end
                
                if flagr
                    obj.angtype = 2;
                else
                    obj.angtype = 1;
                end
            end
        end % function
        
        function ind = getInd(obj, pp_data)
            % Filter the indices in the GPS physical parameter dataset.
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
            
            if ~isstruct(pp_data) || ~isfield(pp_data,'TX_el') || ~isfield(pp_data,'RX_el')
                error('filter given arg that is not a GPS physical parameters struct');
            end
            
            if obj.angtype == 1
                % filter via tx angle
                ind = find(pp_data.TX_el <= obj.ang_thresh);
            else
                % filter via rx angle
                ind = find(pp_data.RX_el <= obj.ang_thresh);
            end
            
        end % function
        
    end % methods
    
end % classdef