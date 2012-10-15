classdef FilterGpsData
    % An abstract class for filtering GPS measurement data from makeGpsData.
    %
    % This is the top-level class for a Template design pattern that can
    % filter GPS data.  More types of filters can easily be created in the
    % future without altering existing classes.  Sub-classes only need to
    % implement the getInd method.
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
        ind = getInd(obj, gps_data, prn);
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
        % Throws an error if any argument is invalid or incorrectly
        % constructed.
    end
    
    %% Methods
    methods (Access=public)
        function data_out = filter(obj, gps_data)
            % Filter the GPS dataset struct.
            %
            % Filtering is performed on each GPS PRN in the data struct.  
            % The metadata is preserved.
            %
            % Outputs: filtered GPS measurement data (could be empty of 
            %   all data points or even contain all of the original data)
            %
            % Inputs:
            %   obj: this instance
            %   gps_data: (single) struct of GPS data defined by makeGpsData, it
            %       can contain zero, one, or many sets of GPS PRN data.
            
            if isempty(gps_data)
                gps_data = makeGpsData; % always output a valid struct
            end
            data_out = copyGpsMetaData(gps_data);
            
            if ~FilterGpsData.hasData(gps_data)
                return; % no usable data
            end
            
            % filter and copy data for each PRN in the dataset
            f_ind = 1;
            for i = 1:length(gps_data.GPS_PRN)
                ind = obj.getInd(gps_data, i);
                if ~isempty(ind)
                    % copy filtered data
                    data_out.GPS_PRN(f_ind) = gps_data.GPS_PRN(i);
                    data_out.PRN_data{f_ind} = struct('epoch',[],'raw_SNR',[],'pseudorange',[],'doppler',[],'phase',[]);
                    data_out.PRN_data{f_ind}.epoch = gps_data.PRN_data{i}.epoch(ind);
                    data_out.PRN_data{f_ind}.raw_SNR = gps_data.PRN_data{i}.raw_SNR(ind);
                    data_out.PRN_data{f_ind}.pseudorange = gps_data.PRN_data{i}.pseudorange(ind);
                    data_out.PRN_data{f_ind}.doppler = gps_data.PRN_data{i}.doppler(ind);
                    data_out.PRN_data{f_ind}.phase = gps_data.PRN_data{i}.phase(ind);
                    f_ind = f_ind + 1;
                end
            end % for
            
        end % function
    end % methods
    
    %% Static Methods
    methods (Static=true)
        function b = hasData(gps_data)
            % Determines if the GPS dataset struct has usable data to filter.
            %
            % Outputs: true=contains data, false=no data
            %
            % Inputs:
            %   obj: this instance
            %   gps_data: (single) struct of GPS data defined by makeGpsData
            if (~isstruct(gps_data)) || isempty(gps_data.PRN_data) || ...
                    isempty(gps_data.GPS_PRN) || ...
                    ((length(gps_data.GPS_PRN)) == 1 && gps_data.GPS_PRN(1) == -1)
                b = false;
            else
                b = true;
            end
        end
    end
    
end % classdef