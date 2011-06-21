classdef FilterPpData
    % An abstract class for filtering GPS physical parameter data from gps_phys_params.
    %
    % This is the top-level class for a Template design pattern that can
    % filter GPS physical parameter data.  More types of filters can easily
    % be created in the future without altering existing classes.
    % Sub-classes only need to implement the getInd method.
    %
    % Note, his is a MATLAB Value Class.  It is designed so that you can
    % create copies of it and manipulate each separately.
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

    %% Abstract Methods
    methods (Abstract,Access=public)
        ind = getInd(obj, pp_data); 
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
    end
    
    %% Methods
    methods (Access=public)
        function data_out = filter(obj, pp_data)
            % Filter a GPS physical parameter dataset struct. 
            %
            % The metadata is preserved.  Returns a blank GPS physical 
            % parameters struct if pp_data is [].
            %
            % Outputs: filtered GPS physical parameters (could be empty of
            %   all data points or even contain all of the original data)
            % Inputs:
            %   obj: this instance
            %   pp_data: (single) struct of GPS physical parameter data
            %       defined by gps_phys_params.
            %
            % Throws an exception if the argument is not a GPS physical
            % parameters struct.
            
            data_out = makePpData; % always output a valid struct
            
            if isempty(pp_data)
                return;
            end
            
            % arg check:
            if ~isstruct(pp_data) || ~isfield(pp_data,'epoch')
                error('filter given arg that is not a GPS physical parameters struct');
            end
            
            data_out.meta = pp_data.meta; % copy the metadata
            
            if isempty(pp_data.epoch)
                return; % no usable data
            end
            
            % filter the dataset
            ind = obj.getInd(pp_data);
            if ~isempty(ind)
                % copy filtered data
                data_out.epoch = pp_data.epoch(ind);
                data_out.TX_az = pp_data.TX_az(ind);
                data_out.TX_el = pp_data.TX_el(ind);
                data_out.RX_az = pp_data.RX_az(ind);
                data_out.RX_el = pp_data.RX_el(ind);
                data_out.range = pp_data.range(ind);
                data_out.range_rate = pp_data.range_rate(ind);
                data_out.GPS_yaw = pp_data.GPS_yaw(ind);
            end
        end % function
        
        
        function [data_out, gps_out] = filterMeas(obj, pp_data, gps_data)
            % Filter GPS physical parameters and GPS measurement data.
            %
            % Filtering occurs on directly on the GPS physical parameters 
            % and the GPS measurement data is filtered to reflect only the 
            % contents of the filtered physical parameters.
            % The metadata is of both is preserved.
            %
            % Outputs:
            % data_out: filtered GPS physical parameters (could be empty of
            %   all data points or even contain all of the original data)
            % gps_out: filtered GPS measurement data (could be empty of
            %   all data points or even contain all of the original data),
            %   empty, [], if no gps_data is provided.
            %
            % Inputs:
            %   obj: this instance
            %   pp_data: (single) struct of GPS physical parameter data
            %       defined by gps_phys_params.
            %   gps_data: (single) struct of GPS data defined by
            %       makeGpsData, it will only contain the GPS PRN data for
            %       the pp_data.  (Optional, this function behaves as
            %       obj.filter(pp_data) if not provided.)
            %
            % Throws an exception if the argument either argument is not an
            % appropriate struct type, or if the RX_ID and GPS PRN do not
            % match between the arguments (these are from separate
            % datasets).
            
            gps_out = []; % default
            
            if isempty(pp_data)
                data_out = makePpData;
                return;
            end
            
            % first, filter by the physical parameters
            data_out = obj.filter(pp_data);
            
            if nargin > 2
                % also filter gps_data by the physical param filter results
                
                % check that these two datasets are consistent on receiver
                % ID & GPS PRN - they must be from the same dataset!
                if pp_data.meta.RX_meta.RX_ID ~= gps_data.RX_meta.RX_ID
                    error('Receiver ID (RX_ID) in the GPS measurements doesn''t match the physical params, Aborting.');
                end
                
                if ~any(gps_data.GPS_PRN == pp_data.meta.TX_ID.GPS_PRN)
                    error('GPS PRNs in GPS measurements don''t match the physical params, Aborting.');
                end
                
                % We'll filter by explicit times and by GPS PRN since the
                % filtered phys params only pertain to one PRN
                f1 = FilterGpsPrn('include', data_out.meta.TX_ID.GPS_PRN);
                
                % call convertTime from UTC to GPS here
                gpstimes = convertTime('GPS','UTC',data_out.epoch); % gps epoch
                
                % filter the times in the GPS epoch, use a 1msec tolerance
                % to handle convertTime calc errors.
                f2 = FilterGpsExplicitTimes(gpstimes, 0.001/86400);
                
                f = FilterGpsMulti(f1,f2);
                gps_out = f.filter(gps_data);
            end
        end % function
        
    end % methods
    
end % classdef