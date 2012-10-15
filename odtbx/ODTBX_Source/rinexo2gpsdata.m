function [varargout] = rinexo2gpsdata(rinexfile, RX_ID, filter, memszlimit, fdir)
%RINEXO2GPSDATA Read RINEX 2.x formatted GPS observation data into a makeGpsData struct.
%
% This function reads in the given RINEX observation file and stores the
% results into one or more GPS measurement data structs (see makeGpsData).
% The RINEX observation file header is stored in the
% RX_meta.obs_metadata field of the GPS measurement data structs.
%
% This function can handle extremely large files that might exceed MATLAB
% memory by instead storing the processed GPS observation data as a series
% of .mat files instead of returning a MATLAB struct.  The default
% behavior of this function returns a single gps_data struct when no
% memszlimit argument is supplied.  However, when this argument is supplied
% it is used to split up the observation data into a series of .mat files
% (at least one file).  The full file pathnames are returned as a cell
% array.  An optional directory argument can be used to specify the file
% locations.
%
% The GPS measurements can optionally be filtered using a FilterGpsData
% instance.  If multiple filters are desired then use the FilterGpsMulti
% class which is optimized for performance.  This filter argument can
% be left empty.  This function will read very large observation files a
% chunk at a time and apply the filter as it processes.
%
% Note: This reader will only interpret certain RINEX GPS observation types
% (C1,L1,D1,S1)
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   rinexfile       char    1xN     String specifiying the input file
%   RX_ID           double  1       User-defined unique receiver identifier
%                                   (required).  This should be unique to
%                                   the entire receive system (antenna to
%                                   measurement output) and may be used by
%                                   other tools to ensure the correct data
%                                   is being used together.
%   filter          object  1       (Optional) An object of FilterGpsData
%                                   type (several sub-classes are
%                                   available, including FilterGpsMulti
%                                   which can apply multiple FilterGpsData
%                                   filters).
%   memszlimit      double  1       (Optional) The maximum number of megabytes
%                                   that a GPS measurement struct should
%                                   occupy in memory.  If supplied, this
%                                   function returns a cell array of the
%                                   full file pathnames written during
%                                   processing.  (Note, a terrestrial
%                                   receiver might see ~52MB/day,
%                                   consistently 16 SVs at 1 sec sampling.)
%                                   This limit is applied after filtering
%                                   and is a guideline rather than a hard
%                                   limit.  Must be greater than zero.
%   fdir            char    1xN     (Optional, with memszlimit), the
%                                   directory location to save the files.
%
%   OUTPUTS
%   varargout, one item of either:
%   gps_data        struct  1       A struct with the GPS measurements, see
%                                   makeGpsData for description
%   OR
%   filenames       cell    1xN     A cell array of full file pathnames for
%                                   each .mat file written by this
%                                   function.  Each .mat file contains a
%                                   gps_dat struct with a subset of the
%                                   data from the RINEX observation file.
%
% See also: makeGpsData, RinexOReader
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

gps_data = makeGpsData();
fn = cell(0); % filename cell array output

%% Arg checks
dofilter = 0;
if exist('filter','var') && ~isempty(filter)
    if ~isa(filter,'FilterGpsData')
        error('gps_phys_params: Given invalid filter type, aborting.');
    else
        dofilter = 1;
    end
end

dofiles = 0; % 0= don't process as files, 1=do process as files
if exist('memszlimit','var') && ~isempty(memszlimit)
    if memszlimit <= 0
        error('rinexo2gpsdata: memszlimit illegal value (%d).', memszlimit);
    else
        dofiles = 1;
        
        if exist('fdir','var') && ~isempty(fdir)
            % test the user-supplied dir:
            if ~exist(fdir,'dir')
                error('rinexo2gpsdata: given unknown directory (%s).', fdir);
            end
            if fdir(end) ~= filesep
                % add a trailing filesep
                fdir = [fdir filesep];
            end
        else
            % use the temp dir
            fdir = tempdir;
        end
    end
end

%% open file & read header
rr = RinexOReader(rinexfile);
if ~rr.isReady || rr.isDone
    error('Failed to open file: %s\n',rinexfile);
end

% Set a read limit on the file reader to be somewhat consistent with the
% memszlimit, but don't allocate too much.
if dofiles
    rl = ceil((memszlimit *1024 * 1024)/(5*8)); % roughly consistent
    if rl > 86400
        rr.readlimit = 86400;
    else
        rr.readlimit = rl;
    end
end

% parse the header
rr.readHeader;
if ~rr.isReady || rr.isDone
    error('Failed to parse header of file: %s\n',rinexfile);
end

% populate the metadata:
gps_data.RX_meta.meas_file = rinexfile;
gps_data.RX_meta.RX_ID = RX_ID;
gps_data.RX_meta.obs_metadata = rr.header;

% earliest and latest epochs seen during parsing:
tearly = Inf;
tlate = 0;

%% parse the body
while (rr.isReady && ~rr.isDone)
    rr.readBody;
    if ~isempty(rr.timegps)
        timegps = rr.timegps;
        C1 = rr.C1;
        D1 = rr.D1;
        L1 = rr.L1;
        S1 = rr.S1;
    else
        error('Failed to parse data from file: %s\n',rinexfile);
    end
    
    % find which time indices have at least some data for each PRN
    testvec = (C1 ~= 0) | (L1 ~= 0) | (S1 ~= 0) | (D1 ~= 0);
    
    % fill in the gps_data PRN data, find out where to start filling in new
    % data when we have some
    prnlen = length(gps_data.GPS_PRN);
    if prnlen == 1 && gps_data.GPS_PRN(1) == -1
        % not really populated, so start in index 1
        prnind = 1;
    else
        % start after the last index
        prnind = prnlen+1;
    end    
    for i = 1:32
        tind = find(testvec(i,:));
        if ~isempty(tind)
            if all(gps_data.GPS_PRN ~= i)
                % fill in new data for this PRN
                gps_data.GPS_PRN(prnind) = i;
                gps_data.PRN_data{prnind}.epoch =       timegps(tind);
                gps_data.PRN_data{prnind}.raw_SNR =     S1(i,tind);
                gps_data.PRN_data{prnind}.pseudorange = C1(i,tind);
                gps_data.PRN_data{prnind}.doppler =     D1(i,tind);
                gps_data.PRN_data{prnind}.phase =       L1(i,tind);
                
                % remember our earliest and latest epochs
                emax = max(gps_data.PRN_data{prnind}.epoch);
                emin = min(gps_data.PRN_data{prnind}.epoch);
                if emax > tlate
                    tlate = emax;
                end
                if emin < tearly 
                    tearly = emin;
                end
                
                % bump the index
                prnind = prnind + 1;
            else
                % add to any current data
                sameind = find(gps_data.GPS_PRN == i);
                gps_data.PRN_data{sameind}.epoch =       [gps_data.PRN_data{sameind}.epoch       timegps(tind)];
                gps_data.PRN_data{sameind}.raw_SNR =     [gps_data.PRN_data{sameind}.raw_SNR     S1(i,tind)];
                gps_data.PRN_data{sameind}.pseudorange = [gps_data.PRN_data{sameind}.pseudorange C1(i,tind)];
                gps_data.PRN_data{sameind}.doppler =     [gps_data.PRN_data{sameind}.doppler     D1(i,tind)];
                gps_data.PRN_data{sameind}.phase =       [gps_data.PRN_data{sameind}.phase       L1(i,tind)];
                
                % remember our earliest and latest epochs
                emax = max(gps_data.PRN_data{sameind}.epoch);
                emin = min(gps_data.PRN_data{sameind}.epoch);
                if emax > tlate
                    tlate = emax;
                end
                if emin < tearly 
                    tearly = emin;
                end
            end
        end
    end % for
    
    %% filter the data
    if dofilter
        gps_data = filter.filter(gps_data);
    end
    
    %% reorder the PRN data
    if gps_data.GPS_PRN > 1
        
%         % Debugging:
%         fprintf(1,'Sorting GPS data required:\n');
%         fprintf(1,'Starting PRNs: %s\n',num2str(gps_data.GPS_PRN));
%         for i = 1:length(gps_data.GPS_PRN)
%             fprintf(1,'PRN %d,\tepoch length: %d\n',gps_data.GPS_PRN(i), length(gps_data.PRN_data{i}.epoch));
%         end
        
        [prns, newind] = sort(gps_data.GPS_PRN);
        if any(newind ~= 1:length(newind))
            
            tmp = cell(1,length(prns));
            ids = zeros(1,length(prns));
            
            % apply the new sort indices to the PRN data
            for i = 1:length(newind)
                tmp{i} = gps_data.PRN_data{newind(i)};
                ids(i) = gps_data.GPS_PRN(newind(i));
            end
            gps_data.GPS_PRN = ids;
            gps_data.PRN_data = tmp;
        end
        
%         % Debugging:
%         fprintf(1,'Sorting GPS data finished:\n');
%         fprintf(1,'Ending PRNs: %s\n',num2str(gps_data.GPS_PRN));
%         for i = 1:length(gps_data.GPS_PRN)
%             fprintf(1,'PRN %d,\tepoch length: %d\n',gps_data.GPS_PRN(i), length(gps_data.PRN_data{i}.epoch));
%         end
    end
    
    %% parse as files
    if dofiles
        while 1
            % Check the ratio compared to our limit, if under, keep
            % parsing.  If over, split out by time up to the size limit.
            
            szratio = gpsstructsize(gps_data)/(memszlimit *1024 * 1024); % bytes
            if szratio < 1.0
                break; % keep reading
            else
                % create a time window from the earliest epoch to a good
                % guess based on time and the memory ratio
                ta = tearly;
                tb = ((tlate-tearly)/szratio)+tearly;
                
                % filter by this time window
                f = FilterGpsTimeWindow(ta, tb);
                gps_dat = f.filter(gps_data);
                
                % save them, 
                fname = createfilename(gps_dat, ta, tb, RX_ID);
                filepath = [fdir fname];
                save(filepath, 'gps_dat');
                fn{end+1} = filepath; %#ok<AGROW>
                
                % retain the rest for further reading
                tearly = tb;
                f = FilterGpsTimeWindow(tb, tlate);
                gps_data = f.filter(gps_data);
                
                % for user status
                fprintf(1,'rinexo2gpsdata: saving PRNs: [%s]\n\t\t to file: %s\n',num2str(gps_dat.GPS_PRN),filepath);
                
            end
        end % while
    end
    
end % while

%% outputs
if dofiles
    % save any remaining data to a file
    gps_dat = gps_data;
    clear gps_data;
    fname = createfilename(gps_dat, tearly, tlate, RX_ID);
    filepath = [fdir fname];
    save(filepath, 'gps_dat');
    fn{end+1} = filepath;

    % for user status
    fprintf(1,'rinexo2gpsdata: saving PRNs: [%s]\n\t\t to file: %s\n',num2str(gps_dat.GPS_PRN),filepath);

    % output the filenames
    varargout{1} = fn;
else
    % output the struct
    varargout{1} = gps_data;
end

end % function

function sz = gpsstructsize(gps_data)
% Helper function to determine the size of the gps data struct, in
% bytes.

% we'll ignore metadata and assume its small

sz = length(gps_data.GPS_PRN) * 8; % count of bytes
if length(gps_data.GPS_PRN) == 1 && gps_data.GPS_PRN(1) == -1
    return;
end
for i = 1:length(gps_data.GPS_PRN)
    sz = sz + gpsprnsize(gps_data.PRN_data{i}.epoch);
end
end % function

function sz = gpsprnsize(epochdata)
% Helper function to determine the size of a single PRN's data in a gps
% data struct, in bytes, given the epoch dataset.
sz = length(epochdata) * 5 * 8;
end % function

function fn = createfilename(gps_dat, ta, tb, RX_ID)
% Helper function to create a standardized filename based on times and PRNs

% the PRN strings
minprn = min(gps_dat.GPS_PRN);
maxprn = max(gps_dat.GPS_PRN);
% time date/time stamps, in ISO 8601 to make good filenames
% (the bias is from convertTime, for the GPS epoch bias)
sa = datestr(ta+723186, 30); 
sb = datestr(tb+723186, 30); 
if minprn == maxprn
    fn = sprintf('%d_%s-%s_%d.mat', RX_ID, sa, sb, minprn);
else
    fn = sprintf('%d_%s-%s_%d-%d.mat', RX_ID, sa, sb, minprn, maxprn);
end
end % function

