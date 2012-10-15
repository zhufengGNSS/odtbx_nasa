classdef FilterEstCnoThresh
    % Filters GPS data structs against an estimated C/No threshold.
    %
    % Any data for any PRN in the dataset will be preserved if its C/No is
    % greater than or equal to the value provided when the class is created,
    % the cno_thresh property.  This class can operate on the CN0 output
    % from gps_est_cno, as well as the corresponding GPS physical
    % parameters data, and the original GPS measurement data.
    %
    % Usage:
    % Construct an object of this class with a C/No threshold as the
    % argument:
    %   "f = FilterEstCnoThresh(22)".
    %
    % Then call the filter method of this object on the estimated C/No
    % output data from
    %   "CN0 = gps_est_cno(...);"
    %   "fil_CN0 = f.filter(CN0)"
    %
    % or filter both the estimated C/No and the GPS physical parameters it
    % is based on:
    %   "CN0 = gps_est_cno(phys_param, ...);"
    %   "[fil_CN0, fil_physparam] = f.filter(CN0, phys_param)"
    %
    % or additionally filter the original GPS measurement data:
    %   "phys_param = gps_phys_params(gps_meas, ...);"
    %   "CN0 = gps_est_cno(phys_param, ...);"
    %   "[fil_CN0, fil_physparam, fil_meas] = f.filter(CN0, phys_param, gps_meas)"
    %
    % See also: gps_est_cno, gps_phys_params, convertTime
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
        function obj = FilterEstCnoThresh(threshold)
            if nargin > 0
                if ~isnumeric(threshold) || length(threshold) > 1
                    error('FilterEstCnoThresh passed an invalid argument, a scalar is required.');
                end
                obj.cno_thresh = threshold;
            end
        end
        
        function [varargout] = filter(varargin)
            % Filter C/No data, and optional associated data structs.
            %
            % Filter the C/No data against the threshold provided at
            % construction.  If GPS physical parameters are provided then
            % they will be filtered against the C/No filter result.  If GPS
            % measurement data is provided then it will also be filtered
            % against the C/No filter result. The metadata is of both
            % GPS structs are preserved.
            %
            % This filter assumes the datasets of any provided data sets
            % are actually properly related.  If the times do not properly
            % match then unexpected results may occur.
            %
            % Inputs:
            %   arg 1: this instance, object (required)
            %   arg 2: 2xN double array of the estimated C/No from
            %       gps_est_cno.  (required)  
            %       Description:
            %       (1,:) time, utc epoch of the measurement,
            %       as a matlab serial date, see convertTime
            %       (2,:) signal carrier to noise ratio [db]
            %   arg 3: (optional) struct of GPS physical parameter data
            %       defined by gps_phys_params.
            %   arg 4: (optional) struct of GPS data defined by
            %       makeGpsData.
            %
            % Outputs:
            %   out 1: filtered C/No, same description as input arg 2, may
            %       be empty
            %   out 2: (optional) the filtered GPS physical parameters, if
            %       provided in input arg 3
            %   out 3: (optional) the filtered GPS measurement data, if
            %       provided in input arg 4
            %
            % Throws an exception if the any argument is not an
            % appropriate type, or if the underlying metadata in the
            % GPS structs do not match.
            
            varargout = {[] [] []};
            if nargin < 2 || nargout == 0
                return;
            end
            
            obj = varargin{1};
            if ~isa(obj,'FilterEstCnoThresh')
                error('Incorrect object passed, expected FilterEstCnoThresh type.');
            end
            
            cno = varargin{2};
            if ~isnumeric(cno)
                error('Incorrect type passed for C/No, expected numeric array.');
            end
            
            if isempty(cno)
                return;
            end
            
            if size(cno,1) ~= 2
                error('Incorrect dimensions on C/No argument, expected 2xN (got %dx%d).',size(cno));
            end
            
            cnoind = (cno(2,:) >= obj.cno_thresh);
            if isempty(cnoind)
                varargout{1} = zeros(2,0);
            else
                varargout{1} = cno(:,cnoind);
            end
            
            if nargin > 2 && nargout > 1
                % filter phys params, and gps meas if we have them too
                
                % arg check on phys params
                pp_data = varargin{3};
                if ~isstruct(pp_data)
                    error('Incorrect type on physical param argument, expected struct');
                end
                
                % construct a physical params filter based on the filtered
                % C/No times (should be an exact match from the physical
                % params)
                f = FilterPpExplicitTimes(cno(1,cnoind));
                
                if nargin > 3 && nargout > 2
                    % also filter the GPS data at the same time
                    
                    % arg check on gps meas struct
                    gps_data = varargin{4};
                    if ~isstruct(gps_data)
                        error('Incorrect type on GPS measurement argument, expected struct');
                    end
                    
                    [varargout{2}, varargout{3}] = f.filterMeas(pp_data, gps_data);
                else
                    % only apply the filter to the phys param data
                    varargout{2} = f.filter(pp_data);
                end
            end
            
        end % function
    end
    
end % classdef
