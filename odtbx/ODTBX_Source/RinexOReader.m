classdef RinexOReader < handle
    % Reads a RINEX 2.x formatted GPS observation data file.
    %
    % This handle class reads the header and interprets observations of a
    % RINEX v2.x observation file based on the "# / TYPES OF OBSERV" header.
    % However it only interprets certain RINEX GPS observation types
    % (C1,L1,D1,S1).  This class is designed to read very large datafiles
    % by reading an observation file's body in groups, if required.
    %
    % Usage:
    % Construct the class with the RINEX observation file name and path as in:
    %       ror = RinexOReader('navigator.rnx');
    %
    % The class has two methods to tell the client if the file can be read,
    % if more reading is required, or if the file has been fully read.
    % They are:
    %       ror.isReady, and ror.isDone.
    %
    % The readHeader method can be called separately from the readBody.
    % While this is not required it may be useful to determine if the data
    % file can be processed by this class.  Once the header is read one or
    % more readBody methods can be called until the file is completely read,
    % isDone == 1, or until a parse error is encountered (or the file is
    % completely read), isReady == 0.  Each time the readBody methods
    % completes the newly-read measurements are stored in the following
    % properties for reading:
    %   timegps, C1, D1, L1, S1, and bias.
    % These properties are cleared and overwritten with each readBody call.
    % Note that missing data for measurement times in the RINEX file are 
    % reported as zeros.  The number of lines read during each readBody
    % call are controlled by the readlimit property.
    %
    % Assumptions regarding observation data:
    %   The last two digits of each value are LLI, Sig Strength (these are
    %   dropped).
    %
    % Special handling for GPS times (timegps property):
    % Note that GPS time outputs are based on the GPS epoch of Jan 6, 1980
    % and not the standard datenum epoch.  The following correction
    % is used to be consistent with ODTBX's convertTime method:
    % datenum('Jan 6 1980') = 723186 (exact).
    %
    % See also: convertTime
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
    
    % Original MATLAB .m function, read_rnxo.m:
    % Written by Mike Moreau: 2/07/2001
    % Modified: 2005
    % This class:
    % Written by Allen Brown, 2011, for ODTBX
    
    %% properties
    properties (GetAccess=public, SetAccess=public)
        readlimit = 86400; % the max number of body lines to read during readBody()        
    end % properties
    
    properties (GetAccess=public, SetAccess=protected)
        timegps = []; % double, 1xN, measurement time tag, MATLAB datenum (GPS epoch)
        C1=[]; % double 32xN, GPS L1 psuedorange (first index is the PRN number)
        L1=[]; % double 32xN, GPS L1 carrier phase (first index is the PRN number)
        D1=[]; % double 32xN, GPS L1 doppler (first index is the PRN number)
        S1=[]; % double 32xN, GPS L1 signal strength (first index is the PRN number)
        bias=[]; % double 1xN, clock offset (optional)
        numobs=0; % the number of observation types (C1, L1, D1, S1)
        header=cell(1); % cell 1xJ The actual header lines from the file
        fn = ''; % the file path and name to read
    end % properties
    
    properties (Access=protected)
        fid = -1; % file ID
        readState = 0; % 0=ready to read, 1=header read, body pending, 2=body read, more pending, 3=read complete, 4=can't read/error
        % typeind is a mapping of order to measurement type:
        % typeind(1)=C1, typeind(2)=D1, typeind(3)=L1, typeind(4)=S1;
        % order values: 1-4
        typeind=zeros(1,4);
        epoch = [];
        prns = [];
        count = 0; % the number of lines read during the last readBody() call
    end
    
    %% methods
    methods % public
        function obj = RinexOReader(fname)
            % Constructor that opens the file for reading.
            
            obj.fn = fname;
            
            obj.fid = fopen(obj.fn,'r');
            if obj.fid == -1
                error('RinexOReader: Could not open file: %s\n',obj.fn);
            end
        end % function
        
        function delete(obj)
            % Destructor function that closes any open file handles.
            try
                fclose(obj.fid);
                obj.fid = -1;
            catch %#ok<CTCH>
                % intentionally ignored
            end
        end % function
        
        function ready = isReady(obj)
            % Reader is ready to parse the file (1), or not able to parse the file (0).
            ready = (obj.fid ~= -1) && (obj.readState < 4);
        end % function
        
        function done = isDone(obj)
            % Reader is done parsing the file (1), or not yet done parsing the file (0).
            done = (obj.fid == -1) || (obj.readState > 2);
        end % function
        
        function readHeader(obj)
            % Read header text and set up measurements
            
            if obj.readState > 0
                if obj.readState > 3
                    warning('RinexOReader: unable to read header.');
                end
                return;
            end
            
            EOH = 0;
            while (feof(obj.fid) == 0) && (EOH == 0)
                linestr = fgetl(obj.fid);
                obj.header{end+1} = linestr;
                head_desc = linestr(61:end);
                % look for the number of observations
                if (length(head_desc)>=19 && strcmp('# / TYPES OF OBSERV',head_desc(1:19)))
                    obs_line = linestr(1:60);
                end
                % Look for the end of the header
                if (length(head_desc)>=13 && strcmp('END OF HEADER',head_desc(1:13)))
                    EOH = 1;
                end
            end
            % at end of loop, next line read will be first line of file body
            
            % Interpret measurement types contained in header, set the
            % mapping
            obj.numobs = strread(obs_line,'%f',1);
            compressedstr = sscanf(obs_line,'%s');
            obj.typeind = zeros(1,4);
            ind = 2:3;
            for i = 1:obj.numobs
                label = compressedstr(ind);
                switch (label)
                    case 'C1'
                        obj.typeind(1) = i;
                    case 'D1'
                        obj.typeind(2) = i;
                    case 'L1'
                        obj.typeind(3) = i;
                    case 'S1'
                        obj.typeind(4) = i;
                end
                ind = ind+2;
            end
            
            if obj.numobs > 5
                fclose(obj.fid);
                obj.fid = -1;
                error('not set up to read more than 5 observations types')
            end
            
            obj.readState = 1;
            
        end % function
        
        function readBody(obj)
            % Read the body portion of the file until hitting the storage limit.
            %
            % Parses the body of the RINEX observation file into
            % measurements.  The following arrays are re-filled with data:
            % timegps, C1, D1, L1, S1, bias.  Use isReady() and isDone() to
            % determine if this function can be called again to obtain more
            % data.
            
            if obj.fid == -1
                obj.readState = 4;
                warning('RinexOReader: unable to read body.');
                return;
            end
            
            if obj.readState < 1
                obj.readHeader;
            end
            
            if obj.readState ~= 1 && obj.readState ~= 2
                return;
            end
            
            obj.epoch = zeros(8,obj.readlimit);
            obj.prns = zeros(32,obj.readlimit);
            obj.bias = zeros(1,obj.readlimit);
            
            %for i = 1:obj.numobs
            %eval(['obj.',obj.types(i,:),' = zeros(32,86400);']);
            %end
            if obj.typeind(1) ~= 0
                obj.C1 = zeros(32,obj.readlimit);
            else
                obj.C1 = [];
            end
            if obj.typeind(2) ~= 0
                obj.D1 = zeros(32,obj.readlimit);
            else
                obj.D1 = [];
            end
            if obj.typeind(3) ~= 0
                obj.L1 = zeros(32,obj.readlimit);
            else
                obj.L1 = [];
            end
            if obj.typeind(4) ~= 0
                obj.S1 = zeros(32,obj.readlimit);
            else
                obj.S1 = [];
            end
            
            % Read data
            % Start loop through time
            obj.count = 0;
            
            % Read the first line following the end of header line
            linestr = fgetl(obj.fid);
            
            while feof(obj.fid) == 0
                obj.count = obj.count + 1;
                
                % Read first line of record with epoch, numobs, prns (preceded by signal identifier), bias (optional)
                %  98  6 21  2  0 41.0162948  0 12G26G23G 8G17G 5G10G 9G30G21G25G 1G 6-0.148499620
                %  format string: 1X,I2,4(1X,I2),F11.7,2X,I1,I3,12(A1,I2),F12.9
                %
                % (if there are more than 12 SVs, the additional prn ids follow on the next line
                %  format string: 32X,12(A1,I2)
                obj.epoch(:,obj.count) = str2num(linestr(1:32))';
                numsvs = obj.epoch(8,obj.count);
                % If there is data present... more than zero svs
                if numsvs > 0
                    % Remove satellite system identifier (assumed GPS signals), replace with space
                    svs_str = strrep(linestr(33:min(68,length(linestr))),'G',' ');
                    svs = sscanf(svs_str,'%d');
                    % read clock offset if present
                    if length(linestr)>69
                        obj.bias(obj.count) = sscanf(linestr(69:end),'%f');
                    end
                    % More than 12 svs, so read epoch line
                    if numsvs > 12
                        linestr = fgetl(obj.fid);
                        svs_str = strrep(linestr(33:min(68,length(linestr))),'G',' ');
                        svs = [svs;sscanf(svs_str,'%d')];
                    end
                    
                    obj.prns(svs,obj.count) = svs;
                    
                    % Done reading epoch line(s)
                    % read observations for each prn
                    obs_matrix = zeros(numsvs,obj.numobs);
                    for i = 1:numsvs
                        linestr = fgetl(obj.fid);
                        % The number of observation types is known from the header
                        [data,ct] = sscanf(linestr,'%f');
                        if ct ~= obj.numobs
                            obj.readState= 4;
                            fclose(obj.fid);
                            obj.fid = -1;
                            error('RinexOReader: read %d observations, expected %d at record %d, epoch: %s',ct,obj.numobs,obj.count,calendar2datestr(obj.epoch(1:6,obj.count)));
                        end
                        
                        % Remove the LLI flag and signal strength flag by dropping 4th and 5th decimal places
                        obs_matrix(i,:) = fix(data'.*1000)./1000;
                    end
                    % Write out the observations read at this epoch
                    for i = 1:obj.numobs
                        %eval(['obj.',obj.types(i,:),'(svs,obj.count) = obs_matrix(:,i);']);
                        switch (obj.typeind(i))
                            case 1 % C1
                                obj.C1(svs,obj.count) = obs_matrix(:,i);
                            case 2 % D1
                                obj.D1(svs,obj.count) = obs_matrix(:,i);
                            case 3 % L1
                                obj.L1(svs,obj.count) = obs_matrix(:,i);
                            case 4 % S1
                                obj.S1(svs,obj.count) = obs_matrix(:,i);
                        end
                    end
                    
                end % if numsvs > 0
                
                % bail when hit limit
                if obj.count == obj.readlimit
                    obj.readState = 2; % we will need to continue
                    break;
                end
                
                linestr = fgetl(obj.fid);
            end % while
            
            if feof(obj.fid)
                obj.readState = 3; % done
                fclose(obj.fid);
                obj.fid = -1;
            end
            
            % Find invalid times in epoch block
            ind = find(~([obj.epoch(6,:)<0]|[obj.epoch(6,:)>=60] | ...
                [obj.epoch(5,:)<0]|[obj.epoch(5,:)>=60] | ...
                [obj.epoch(4,:)<0]|[obj.epoch(4,:)>=24] | ...
                [obj.epoch(3,:)<1]|[obj.epoch(3,:)>31] | ...
                [obj.epoch(2,:)<1]|[obj.epoch(2,:)>12] ));
            
            obj.epoch = obj.epoch(:,ind);
            obj.prns = obj.prns(:,ind);
            obj.bias = obj.bias(ind);
            
            %for i = 1:obj.numobs
            %    eval([obj.types(i,:),' = ',obj.types(i,:),'(:,ind);']);
            %end
            obj.C1 = obj.C1(:,ind);
            obj.D1 = obj.D1(:,ind);
            obj.L1 = obj.L1(:,ind);
            obj.S1 = obj.S1(:,ind);
            
            % Compute time in gps format, but do not adjust for leap seconds
            %timegps = utc2gps(obj.epoch(1:6,:)',0)';
            % add the two-digit century prefix:
            obj.epoch(1,:) = (obj.epoch(1,:) > 79)*1900 + obj.epoch(1,:);
            obj.epoch(1,:) = (obj.epoch(1,:) < 79)*2000 + obj.epoch(1,:);
            % convert to the GPS time datenum convention:
            obj.timegps = (datenum(obj.epoch(1:6,:)') - 723186)';
            
            if obj.readState == 2
                fprintf(1,'RinexOReader read %d records, read incomplete.\n',obj.count);
            else
                fprintf(1,'RinexOReader read %d records, read complete.\n',obj.count);
            end
            
        end % function
        
        function set.readlimit(obj,value)
            % Sets (resets) the readlimit property and re-allocates the measurement arrays.
            %
            % Use this method to change the readlimit, which controls how
            % much of a file can be read on a single readBody() call.
            % Lower values reduce memory footprint while reading large
            % files but require more frequent readBody() calls.  A lower
            % limit of 10 lines is enforced.
            %
            % Note: array size changes from this set call take effect on
            % the next readBody() call.
            
            if ~isnumeric(value)
                error('RinexOReader set readlimit given invalid argument');
            end
            if value < 10
                error('RinexOReader set readlimit given unrealistically low value');
            end
            obj.readlimit = value;         
            
        end % function
        
    end % methods
    
end % classdef