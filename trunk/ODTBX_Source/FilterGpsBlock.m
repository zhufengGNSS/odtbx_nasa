classdef FilterGpsBlock < FilterGpsPrn
    % A FilterGpsData class for filtering GPS measurement data by GPS block
    % type.
    %
    % This filter preserves all data for a PRN in the dataset by comparing
    % it to a GPS SV block type classification.  This class can either
    % include or exclude data compared to the block type.  The GPS SV block
    % type information is provided as a set of data from makeGpsTXID, which
    % links GPS PRN number, the GPS block type, and a user-supplied GPS SV
    % ID (the GPS SV ID is not used by this class).
    %
    % Usage:
    % Construct an object of this class with a cell array of GPS TX ID
    % data, the SV block type of interest, and whether to treat that bock
    % type to be "included" or "excluded".  PRNs in the dataset should be
    % unique.  If the given block value to filter on does not match a block
    % value in the GPs TX ID data then this class will behave as "exclude
    % none" or "include none" as appropriate.  Any PRN in the GPS dataset
    % to be filtered that is not listed in the given GPS TX ID data will be
    % treated as an unknown/unspecified PRN, i.e. removed when "include"
    % is used and preserved when "exclude" is used.
    %
    % Example:
    % Assuming we have a set of GPS TX ID data such as:
    % % PRN 12, block 1(=II/IIA), user ID 441
    % txid{1} = makeGpsTXID(441, 12, 1);
    % % PRN 13, block 1(=II/IIA), user ID 442
    % txid{2} = makeGpsTXID(442, 13, 1);
    % % PRN 14, block 4(=IIF), user ID 443
    % txid{3} = makeGpsTXID(443, 14, 4);
    %
    % ex: To include data only from any block II/IIA SV:
    %     "f = FilterGpsBlock('include',txid,1)".
    % ex: To exclude data any block IIF SV (and include any others):
    %     "f = FilterGpsBlock('exclude',txid,4)".
    % Then call the filter method of this object,
    % "filtered_data = f.filter(input_data)".
    %
    % Note the GPS Satellite Block type: is defined in makeGpsTXID and in
    % gpsmeas (reprinted here for convenience):
    %   1=II/IIA (default)
    %   2=IIR
    %   3=IIR-M
    %   4=IIF
    %   5=III
    %  (These numerical values should be consistent with the gpsmeas GPSBlock
    %   option.)
    %
    % See also: makeGpsData, FilterGpsPrn, FilterGpsData, makeGpsTXID, gpsmeas
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
        TXIDs = []; % the list of GPS TX_ID data, 1xG cell array
        
        block_filter = 1; % the GPS block enum value in TXIDs
    end
    
    %% Methods
    methods (Access=public)
        
        function obj = FilterGpsBlock(flag, TXIDs, blocknum)
            % Constructor that takes an include/exclude flag, a set of TXIDs and a block type.
            %
            % flag is either 'include' or 'exclude', no other value is allowed,
            %   defaults to 'include' if empty.
            %
            % The TXIDs are a 1xG cell array of one or more structs that
            %   describe the GPS transmitters in the gps_meas data, the
            %    GPS_PRN fields must be unique.  See makeGpsTXID.
            %
            % The blocknum is the value of the block type that should be
            % included or excluded, see notes above for the enumeration values.
            %
            % Throws an error if any argument is invalid or if all three
            % arguments are not supplied.  Calling with no arguments creates an
            % object that includes no GPS data.
            
            % CONSTRUCTOR PRE-INIT (no refs to obj here!)
            PRNs = [];
            if nargin == 0
                superargs = {};
            elseif nargin > 0
                if nargin ~= 3
                    error('FilterGpsBlock constructor given incorrect number of arguments (expects 3)');
                end
                
                if isempty(TXIDs) || size(TXIDs,1) ~= 1
                    error('FilterGpsBlock constructor given incorrectly sized argument for TXIDs');
                end
                
                if ~iscell(TXIDs)
                    error('FilterGpsBlock constructor given TXIDs that is not a cell array');
                end
                
                % create the PRNs to use from the TXIDs
                PRNvals = zeros(1,length(TXIDs));
                blockvals = zeros(1,length(TXIDs));
                
                for i = 1:length(TXIDs)
                    PRNvals(i) = TXIDs{i}.GPS_PRN;
                    blockvals(i) = TXIDs{i}.block_type;
                end
                superargs{1} = flag;
                superargs{2} = PRNvals(blockvals == blocknum);
            end
            
            % OBJECT INIT
            % Call the superclass constructor with the
            % empty cell array (no arguments) if nargin == 0
            % otherwise cell array is not empty
            obj = obj@FilterGpsPrn(superargs{:}); % superclass constructor
            
            % set local properties
            if nargin > 0
                obj.TXIDs = TXIDs;
                obj.block_filter = blocknum;
            end
            
            % OBJECT POST-INIT
        end
        
    end % methods
    
end % classdef
