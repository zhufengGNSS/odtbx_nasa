function [phys_param, fil_gps_meas] = gps_phys_params(gps_meas, RX_state_source, RX_ID, epoch, tsc, xsc, qatt, Rotation2ECI, ant_bod, gps_alm_file, TX_IDs, filter)
%
% Creates GPS physical parameter data sets from GPS measurements.
%
% gps_phys_params calculates the physical parameters of the relationship
% between GPS spacecraft transmitter antennas and a specified receiver
% "box" and antenna that would be useful for constructing a link budget.
% Note, the term "receiver box" is used because this function is designed
% for a terrestrial receiver or a receiver onboard an orbiting spacecraft.  
% The gps_meas contains the measurement times for one or more GPS spacecraft
% measurements from a receiver.  The supplied receiver box time history of 
% position, velocity (tsc, xsc), and optionally attitude (in qatt), are used
% to interpolate the receiver box state at the measurement times.  The
% relationship between the specified frame (in qatt) and the receiver
% antenna can be set via ant_bod if desired.  The xsc and qatt frame is
% defaulted to ECI if Rotation2ECI is not supplied.  Otherwise Rotation2ECI
% converts from the xsc-qatt frame to ECI.  The GPS constellation
% information is specified in the Yuma almanac file.
%
% This tool performs some cross-checking on the data to ensure the data is
% compatible.  The given TX_IDs allow a GPS PRN will be processed,
% otherwise the tool assumes that PRN is unknown.  Data from gps_meas that 
% doesn't match a GPS_PRN in TX_IDs are not processed.  Additionally the
% tool checks to ensure the given gps_meas data is for the intended
% receiver, specified in RX_ID.
%
% The GPS measurements can optionally be filtered before processing using
% a FilterGpsData instance.  If multiple filters are desired then 
% use the FilterGpsMulti class which is optimized for performance.  The
% filtered GPS measurements are returned as an optional output argument.  
% This filter argument can be left empty.
%
% Note, with an orbiting receiver the sampling of tsc, xsc, and qatt is 
% important.  Good range rate results depend more directly on the sampling
% of the orbit.  Few samples on the orbit will result in worse position or
% velocity interpolation error.  Interpolation of qatt assumes
% a constant angular velocity between any two points.  This assumption also
% is based on well-sampled qatt data where the attitude change between sample
% times is not large.
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   gps_meas        struct  1       A struct to hold raw GPS measurements
%       gps_meas fields:
%   .RX_meta        struct  1       A struct to hold the metadata
%                                   associated with this gps_meas, see
%                                   below for its fields
%   .GPS_PRN        double  1xM     The PRN numbers for the data in this
%                                   struct.  Each GPS_PRN matches once cell
%                                   in the PRN_data.  -1 means "no data"
%   .PRN_data       cell    1xM     Measurement data for a PRN
%       PRN_data fields:
%   .epoch          double  1xT     measurement time, using the ODTBX 
%                                   convention of datenum with convertTime 
%                                   GPS epoch
%   .raw_SNR        double  1xT     raw signal strength measurement
%   .pseudorange    double  1xT     raw pseudorange measurement
%   .doppler        double  1xT     raw doppler measurement
%   .phase          double  1xT     raw phase measurement
%
%   RX_state_source string  1       Identifier for the source of the receiver
%                                   box state, attitude, and antenna attitude.
%                                   This is user-defined and is only used to
%                                   associate the given receiver state data with a
%                                   particular receiver system configuration.
%                                   May be left empty, [].
%   RX_ID           double  1       User-defined unique receiver identifier
%   epoch           double  1       UTC epoch of the receiver box state times, as
%                                   a MATLAB serial date, see ConvertTime
%   tsc             double  1xN	    receiver box state times, UTC (sec from epoch in
%                                   options).  These are NOT the measurement times,
%                                   these define the receiver box trajectory.
%   xsc             double  6xN     receiver box state [position;velocity] (km), with
%                                   tsc.  (Assumed ECI unless Rotation2ECI
%                                   is used.)
%   qatt            double  4xN,[]  (optional) receiver box body quaternion history where the
%                                   quaternion specifies the user receiver box
%                                   attitude with respect to the same coordinate
%                                   frame that the position and velocity are given
%                                   in.  The quaternion is specified as
%                                   [sin(theta/2)*e_vec;cos(theta/2)] where e_vec is the
%                                   unit vector representing the axis of rotation
%                                   and theta is the angle of rotation about that
%                                   axis.
%                                   If the quaternion is not supplied, then the
%                                   receiver box attitude is aligned with the
%                                   coordinate frame that the position and velocity
%                                   are given in, i.e., the default is [0 0 0 1]'.
%                                   Additionally, the orientation matrix
%                                   of the antenna wrt to the receiver box body can
%                                   be specified in ant_bod.
%   Rotation2ECI    @functionname   A function that allows for transformation 
%                                   of the states in xsc (and qatt), from 
%                                   an arbitrary coodinate system to ECI at a
%                                   given time. Returns a 6x6 matrix.
%                                   Input must be a rotation function handle
%                                   that uses time as an input.  The time
%                                   input is sepcified like epoch.  See the
%                                   embedded IdentRot() function for
%                                   details.
%   ant_bod         double  3x3,[]  unit orthogonal matrix representing the
%                                   rotation of the antenna with respect to the
%                                   spacecreaft body.  The Z-axis of the antenna
%                                   frame is the boresite.  If not specified, then
%                                   an identity rotation is assumed.
%   gps_alm_file    string  1       File and path containing the GPS almanac data
%   TX_IDs          double  1xG     cell array of one or more structs that 
%                                   describe the GPS transmitters in the 
%                                   gps_meas data, the GPS_PRN fields must be
%                                   unique.  GPS_PRNs in the gps_meas are not used
%                                   if they are not set in the TX_IDs.  See
%                                   makeGpsTXID.
%   filter          object  1       (Optional) An object of FilterGpsData
%                                   type (several sub-classes are
%                                   available, including FilterGpsMulti
%                                   which can apply multiple FilterGpsData
%                                   filters).
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   phys_param      cell    1xZ,[]  Cell array of GPS physical parameter structs
%                                   each containing data for one GPS PRN/SV from
%                                   the given gps_meas.  Sized Z for the
%                                   remaining measurements after filtering.
%                                   May be empty if no data could be
%                                   processed.
%           each phys_param fields:
%       .epoch      double	1xZ	    time, UTC epoch of the measurement, as
%                                   a MATLAB serial date, see ConvertTime
%       .TX_az      double	1xZ     transmitter antenna azimuth angle (rad)
%       .TX_el      double	1xZ     transmitter antenna boresight elevation
%                                   angle (rad)
%       .RX_az      double	1xZ     receiver antenna azimuth angle (rad)
%       .RX_el      double	1xZ     receiver antenna boresight elevation angle
%                                   (rad)
%       .range      double	1xZ     absolute value of the transmitter-receiver
%                                   range (km)
%       .range_rate	double	1xZ the transmitter-receiver relative
%                                   velocity projected along the LOS (km/s)
%       .GPS_yaw	double	1xZ     the GPS yaw model yaw angle of the
%                                   satellite (deg)
%       .meta       struct	1       metadata for this phys_param
%           meta fields:
%       .RX_meta	struct	1       Measurement metadata from the gps_meas 
%                                   input, see makeGpsData
%       .RX_state_source	string  1   see inputs
%       .TX_ID      struct	1	GPS SV/PRN identifier
%       .TX_state_source	string  1   The GPS almanac file used for
%                                   processing.
%       .gen_date double	1       MATLAB datenum date/time of when this
%                                   data was processed
%
%   fil_gps_meas    struct  1       A struct of the filtered GPS
%                                   measurements from the input argument
%                                   gps_meas.
%
%
% See also: makeGpsData, makeGpsTXID, getgpsmeas, gpsmeas, convertTime,
% FilterGpsData, FilterGpsMulti
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

% Constants
d2r = pi/180; % radians to degrees conversion

%% check inputs
% existence of args
if isempty(gps_meas) || isempty(RX_state_source) || ...
        isempty(RX_ID) || isempty(epoch) || isempty(tsc) || ...
        isempty(xsc) || ...
        isempty(gps_alm_file) || isempty(TX_IDs)
    error('gps_phys_params: One or more required arguments were empty, aborting.');
end

if size(tsc,1) ~= 1
    error('gps_phys_params: Invalid receiver box state time size, aborting.');
end

if size(xsc,1) ~= 6 
    error('gps_phys_params: Invalid receiver box state size, aborting.');
end

if size(tsc,2) ~= size(xsc,2)
    error('gps_phys_params: receiver box state and time sizes do not match, aborting.');
end

if isempty(Rotation2ECI)
    %For coordinate frames other than ECI. This value allows
    %the user to have x in any coordinate frame as long as the
    %rotation to ECI goes with it. Input must be a pointer to a
    %rotation function with time as the input.
    Rotation2ECI = @IdentRot;
end

if ~isempty(qatt)
    % empty qatt is handled below, otherwise cross-check the dimensions
    if size(tsc,2) ~= size(qatt,2)
        error('gps_phys_params: receiver box attitude and time sizes do not match, aborting.');
    end    
end

% check gps_meas struct
if ~isstruct(gps_meas)
    error('gps_phys_params: Invalid GPS measurement struct, aborting.');
end
if ~isfield(gps_meas,'RX_meta')
    error('gps_phys_params: Invalid GPS measurement struct, missing metadata, aborting.');
end
if ~isfield(gps_meas.RX_meta,'RX_ID')
    error('gps_phys_params: Invalid GPS measurement struct, missing RX ID in metadata, aborting.');
end

% cross-check given RX_ID & gps_meas.RX_meta.RX_ID
if RX_ID ~= gps_meas.RX_meta.RX_ID
    error('gps_phys_params: Given receiver ID does not match receiver ID GPS data, aborting.');
end

%% prep
phys_param = cell(1,0); % default to empty cell array, alloc as we go

%% apply filters
if ~exist('filter','var') || isempty(filter)
    fil_gps_meas = gps_meas;
else
    if ~isa(filter,'FilterGpsData')
        error('gps_phys_params: Given invalid filter type, aborting.');
    end
    
    fil_gps_meas = filter.filter(gps_meas);
end

%% post-filtering processing of data
% check fil_gps_meas after filtering
if length(fil_gps_meas.GPS_PRN) < 1
    % warn here because this could occur with upstream filtering
    warning('No data in the GPS measurement data, (overfiltered?).'); %#ok<WNTAG>
    return;
elseif any(diff(sort(fil_gps_meas.GPS_PRN)) == 0)
    % this really shouldn't happen - the must be an error in how fil_gps_meas
    % was constructed.
    error('Multiple identical GPS PRNs in a (badly constructed) GPS measurement data struct.');
end

% How many sets of PRN data in fil_gps_meas
prnlen = length(fil_gps_meas.GPS_PRN);

ppind = 0; % index for phys_param

% PRNs contained in the transmitter IDs
txidprns = zeros(1,length(TX_IDs));
for i = 1:length(txidprns)
    txidprns(i) = TX_IDs{i}.GPS_PRN;
end

% receiver box time:
tsc_utc = tsc/86400 + epoch; % UTC in datenum (matlab serial date)
minsctime = min(tsc_utc); % max s/c state time in utc, datenum
maxsctime = max(tsc_utc); % max s/c state time in utc, datenum

xsc_t = xsc';

%% main processing loop over the GPS PRN data
for prnind = 1:prnlen
    
    fprintf(1,'gps_phys_params: Processing PRN %d...\n',...
        fil_gps_meas.GPS_PRN(prnind));
    
    %% process each gps prn (make a new phys_param for each prn)
    
    if isempty(fil_gps_meas.PRN_data{prnind}.epoch)
        fprintf(1,'gps_phys_params: missing measurements for PRN %d, skipping!\n',...
                fil_gps_meas.GPS_PRN(prnind));
        continue;  % with next prnind
    end
    
    % the physical parameters for this prn, defaulted data
    phys_param_tmp = makePpData(fil_gps_meas);
    phys_param_tmp.epoch = [];
    phys_param_tmp.TX_az = [];
    phys_param_tmp.TX_el = [];
    phys_param_tmp.RX_az = [];
    phys_param_tmp.RX_el = [];
    phys_param_tmp.range = [];
    phys_param_tmp.range_rate = [];
    phys_param_tmp.GPS_yaw = [];
    
    % Find the matching fil_gps_meas.GPS_PRN(prnind) in TX_IDs.GPS_PRN,
    % and use that transmitter ID.
    txidind = find(txidprns == fil_gps_meas.GPS_PRN(prnind));
    if length(txidind) < 1
        warning('No matching GPS PRN was provided in TX_IDs that matches the GPS measurement data (PRN=%d), skipping.',fil_gps_meas.GPS_PRN(prnind)); %#ok<WNTAG>
        continue;  % with next prnind
    elseif length(txidind) > 1
        error('Multiple matching GPS PRNs were provided in TX_IDs that matches the GPS measurement data (PRN=%d), input error, aborting.',fil_gps_meas.GPS_PRN(prnind));
    end
    
    % the metadata
    phys_param_tmp.meta.RX_meta = fil_gps_meas.RX_meta;
    phys_param_tmp.meta.RX_state_source = RX_state_source;
    phys_param_tmp.meta.TX_ID = TX_IDs{txidind};
    phys_param_tmp.meta.TX_state_source = gps_alm_file;
    phys_param_tmp.meta.gen_date = now;
    
    % set the options & params for the inputs to getgpsmeas
    params.num_ant = 1; % Number of receiver antennas
    params.xmit_pattern_dim = 2; % always assume 2D pattern
    params.rec_pattern_dim = 2; % always assume 2D pattern
    params.GPS_SIZE = 1; % the number of GPS satellites
    params.PRN = fil_gps_meas.GPS_PRN(prnind);
    params.doH = 0; % measurement partial not required
    measOptions = odtbxOptions('measurement');
    measOptions = setOdtbxOptions(measOptions, ...
        'epoch', epoch, ... % Time associated with start of simulation, UTC datenum format
        'PrecnNutnExpire', 0.1, ... % one day period until recomputing precession/nutation
        'Rotation2ECI', Rotation2ECI, ...
        'YumaFile', gps_alm_file, ...
        'AntennaOrientation', ant_bod);
    
    %% Interpolate the given receiver box states to the measurement times
    tmeas_gps = fil_gps_meas.PRN_data{prnind}.epoch; % GPS times in convertTime's datenum
    tmeas_utc = convertTime('UTC','GPS',tmeas_gps); % UTC in datenum
    
    % cross check gps meas times vs s/c state times: early check
    if (min(tmeas_utc) < minsctime)
        fprintf(1,'gps_phys_params: some measurement times are earlier than the receiver box data for PRN %d.\n',...
            fil_gps_meas.GPS_PRN(prnind));
        
        % truncate some indices that can't be supported by the state
        fprintf(1,'gps_phys_params: truncating measurements for PRN %d.\n',...
            fil_gps_meas.GPS_PRN(prnind));
        usableind = find(tmeas_utc >= minsctime);
        
        if isempty(usableind)
            fprintf(1,'gps_phys_params: no supported measurements for PRN %d, skipping!\n',...
                fil_gps_meas.GPS_PRN(prnind));
            continue;  % with next prnind
        else
            tmeas_gps = tmeas_gps(usableind);
            tmeas_utc = tmeas_utc(usableind);
        end
    end
    
    % cross check gps meas times vs s/c state times: late check
    if (max(tmeas_utc) > maxsctime)
        fprintf(1,'gps_phys_params: some measurement times are later than the receiver box data for PRN %d.\n',...
            fil_gps_meas.GPS_PRN(prnind));
        
        % truncate some indices that can't be supported by the state
        fprintf(1,'gps_phys_params: truncating measurements for PRN %d.\n',...
            fil_gps_meas.GPS_PRN(prnind));
        usableind = find(tmeas_utc <= maxsctime);
        
        if isempty(usableind)
            fprintf(1,'gps_phys_params: no supported measurements for PRN %d, skipping!\n',...
                fil_gps_meas.GPS_PRN(prnind));
            continue;  % with next prnind
        else
            tmeas_gps = tmeas_gps(usableind);
            tmeas_utc = tmeas_utc(usableind);
        end
    end
    
    lmeas = length(tmeas_utc);
    xmeas = zeros(6,lmeas); % receiver box states at measurement times
    fprintf(1,'gps_phys_params: Interpolating states for PRN %d...\n',...
        fil_gps_meas.GPS_PRN(prnind));

    % lets find the indices of the receiver box time nearest to each
    % measurement time
    q_tmeas_ind = interp1(tsc_utc,1:length(tsc_utc),tmeas_utc,'nearest');

    % adjust the indices so that they are always the same time or
    % earlier than our measurement time
    q_tmeas_ind = q_tmeas_ind - ( tsc_utc(q_tmeas_ind) > tmeas_utc );

    % decrement any index that is at the end of the receiver box time
    % array (for the for the following for loop)
    q_tmeas_ind = q_tmeas_ind - ( q_tmeas_ind == length(tsc_utc) );
    
    if exist('qatt','var') && ~isempty(qatt)
        qattmeas = zeros(4,lmeas); % receiver box attitude at measurement times
    else
        qattmeas = [];
    end
    
    % Convert time to the same time unit of the derivative states (sec)
    % for the interpolator
    tmeas_utc_s = (tmeas_utc - tsc_utc(1))*86400;
    tsc_utc_s   = (tsc_utc - tsc_utc(1))*86400;
    
    % Interpolate the receiver box state to the measurement times.
    % Interpolate the receiver box attitude to the measurement times,
    % but it might not be specified in the args.
    for n=1:lmeas
        
        inda = q_tmeas_ind(n); % first index
        ta = tsc_utc(inda); % first time
        
        % properly decremented above so this isn't a problem:
        tb = tsc_utc(inda+1); % second time
        tmeas = tmeas_utc(n); % measurement time
        tratio = (tmeas-ta)/(tb-ta); % ratio of time from ta to tb, (0.0-1.0)
                
        % Use the Hermite Interpolator.  The JAT version is faster than the
        % MATLAB version, but otherwise their performance should be the 
        % same.
        xmeas(:,n) = jat.matlabInterface.JATIntegrators.HermiteInterpolator(tsc_utc_s,tmeas_utc_s(n),6,xsc_t);
                
        if ~isempty(qattmeas)
            % Interpolate between the quaternions.  We'll assume a constant
            % angular velocity between the two. While this isn't always a
            % good assumption, it should work well enough for
            % regularly-sampled states.
            qattmeas(:,n) = slerp(qatt(:,inda),qatt(:,inda+1),tratio);
        end
    end
    
    %% calculate the physical parameters
    out = getgpsmeas((tmeas_utc-epoch)*86400,xmeas,measOptions,qattmeas,params);
    
    if any(out.health)
        
        if fil_gps_meas.GPS_PRN(prnind) ~= out.prn
            error('gps_phys_params: GPS PRN from getgpsmeas didn''t match provided PRN (PRN=%d), aborting.',fil_gps_meas.GPS_PRN(prnind));
        end
        
        % store the outputs, if there are any healthy outputs
        phys_param_tmp.epoch = tmeas_utc; % re-use the same times
        phys_param_tmp.TX_az = mod(out.TX_az',360)*d2r; % The transmitter azimuth angle (rad)
        phys_param_tmp.TX_el = mod(out.TX_el',360)*d2r; % The transmitter elevation angle (rad)
        phys_param_tmp.RX_az = mod(out.RX_az',360)*d2r; % The receiver azimuth angle (rad)
        phys_param_tmp.RX_el = mod(out.RX_el',360)*d2r; % The receiver elevation angle (rad)
        phys_param_tmp.range = out.range; % range (km)
        phys_param_tmp.range_rate = out.rrate; % range (km/s)
        phys_param_tmp.GPS_yaw = out.GPS_yaw'; % GPS SV yaw angle
    else
        fprintf(1,'gps_phys_params: Almanac ''%s'' does not support PRN %d, skipping!\n',...
                gps_alm_file, fil_gps_meas.GPS_PRN(prnind));
        continue; % with next prnind
    end
    
    if isempty(phys_param_tmp.epoch)
        % final check, no usable data?
        continue; % with next prnind
    end
    
    % store in cell array
    ppind = ppind + 1;
    phys_param{ppind} = phys_param_tmp;
    
end % for prnind

end % function

% Generic template for Rotation2ECI(t).
% Returns a 6x6 transformation matrix that transforms a position and
% velocity in a given frame to ECI at a specified time, t.  Time is
% specified in MATLAB datenum format, in UTC.
function I=IdentRot(t) %#ok<INUSD>
    % This template function always returns an identity transformation.
    I = eye(6);
end