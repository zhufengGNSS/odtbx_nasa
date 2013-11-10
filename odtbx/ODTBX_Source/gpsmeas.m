function [y,H,R,AntLB,dtsv] = gpsmeas(t,x,options,qatt,sv)

% GPSMEAS  Makes GPS based measurements using dynamic C/No tracking threshold
%
% [y,H,R,AntLB] = GPSMEAS(t,x,options,qatt) creates GPS measurements based on the 
% information in OPTIONS. The measurement types that can be returned are
% range and range rate or doppler.  However, range rate and doppler are
% mutually exclusive outputs.
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   epoch               datenum                Time associated with start of simulation
%   useRange            {true(default), false} Return range measurement
%   useRangeRate        {true, false(default)} Return range rate measurement
%   useDoppler          {true, false(default)} Return instantaneous doppler measurement
%   rSigma                                     Measurement covariance
%   GPSBand             {'L1'(default), 'L2'}  L-band carrier
%   YumaFile            filename               GPS almanac file
%   Rotation2ECI        @functionname          A function that allows for transformation of the state
%                                              vector, x, from an arbitrary coodinate system to ECI at a
%                                              given time. Returns a 6x6 matrix.  Input must be a pointer to a
%                                              rotation function that uses time as an input. See the
%                                              embedded IdentRot() function for details.
%   AntennaPointing     1 x num_ant            Specify attitude profile for each antenna
%                                              (1) zenith pointing or (-1) nadir pointing wrt geocentric LVLH
%                                              (2) parallel or (-2) antiparallel to the Earth-Sun vector
%                                              (3) ecliptic north or south (-3)
%                                              (4) fixed with respect to apogee zenith, or (-4) nadir vector apogee 
%                                                selected as the point of highest altitude (ephemeris must include apogee)
%                                              (5) body fore and (-5) aft directions relative to geocentric LVLH
%                                              (6) body port and (-6) starboard directions relative to geocentric LVLH
%                                              NOTE: The AntennaPointing
%                                              parameter is used only if
%                                              the receiver antenna pattern
%                                              specified by AntennaPattern
%                                              is 1-D.  If the pattern is
%                                              2-D, then the attitude of the
%                                              antennas is determined from
%                                              the spacecraft body attitude
%                                              quaternion specified in the
%                                              input variable x as well as
%                                              the AntennaOrientation
%                                              parameter in this options
%                                              structure.
%   AntennaPattern      filenames              Specify antenna pattern for each antenna, existing antennas are:
%                                              sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
%                                              omni.txt                  - zero dB gain,  180 degree half beamwidth
%                                              trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
%                                              ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
%                                              ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
%
%                                              The format for antenna
%                                              patttern files is as
%                                              follows:
%
%                                              % Header lines
%                                              % .
%                                              % .
%                                              % .
%                                              % Header lines
%                                              % angle [deg] gain [db]
%                                              0.0  gain1
%                                              .    .
%                                              .    .
%                                              angN gainN
%
%                                              The comment character for
%                                              antenna pattern files is the
%                                              % character and there can be
%                                              as many header lines as
%                                              desired. The antenna pattern
%                                              is specified as two columns
%                                              of numbers where the first
%                                              column contains the angle
%                                              off the antenna boresight in
%                                              units of degrees and the
%                                              second column contains the
%                                              gain value in units of dB.
%                                              The angular resolution at
%                                              which these data are
%                                              specified is up to the user.
%                                              An example antenna pattern
%                                              file named
%                                              sensysmeas_ant.txt is
%                                              provided in the
%                                              ODTBX_Data\antennae directory.
%
%   AntennaMask         radians                Global receiving antenna mask angle.  The lesser of the 
%                                              maximum defined angle of the pattern and this number used
%   AtmosphereMask      kilometers             Mask altitude above Earth's surface (km)
%                                              Troposphere mask radius ~50 km
%                                              Ionosphere mask radius ~(500-1000 km)
%   NoiseTemp           Kelvin                 System noise temp [K], space pointing antenna = 290
%                                              System noise temp [K], earth pointing antenna = 300
%   AtmAttenuation      dB                     Attenuation due to atmosphere (should be negative) [dB]
%   TransPowerLevel     {1, 2(default), 3}     Transmitter Power Level  (1-minimum, 2-typical, 3-max)
%   TransPowerOffset    dB                     Global transmitter power offset
%   GPSBlock            {1(default), 2, 3, 4}  GPS Satellite Block :
%                                               1=II/IIA, (elevation only)
%                                               2=IIR, (elevation only)
%                                               3=IIR-M, (elevation only)
%                                               4=IIF, (unsupported)
%                                               5=III (elevation only)
%                                               6=Ficticious 2D antenna pattern (az/el) based on IIA
%   TransAntMask        Radians                The actual mask used is the lesser of this mask and the limit of the defined pattern
%                                              mask = 70 deg includes entire defined pattern
%                                              mask = 42 deg includes only main and first side lobes
%                                              mask = 26 deg includes only main lobe
%   ReceiverNoise       dB                     Noise figure of receiver/LNA
%   RecConversionLoss   dB                     Receiver implementation, A/D conversion losses [dB]
%                                              Novatel: L = -4.0,  Plessey: L = -1.5		
%   SystemLoss          dB                     System losses, in front of LNA
%   LNAGain             dB                     LNA gain (trimble pre-amp spec = 42-48)
%   CableLoss           dB                     Cable losses (after LNA)
%   RecAcqThresh        dB-Hz                  Receiver acquisition threshold
%   RecTrackThresh      dB-Hz                  Receiver tracking threshold. Can only utilize a tracking
%                                              threshold that is different from the acquisition threshold 
%                                              when it is run for all times at once (i.e. batch mode).  
%   DynamicTrackRange   dB                     Dynamic tracking range of receiver, or maximum difference  
%                                              in power levels tracked simultaneously. If the difference
%                                              in snrs between two satellites is more that link_budget.DynamicTrackRange,
%                                              the weaker of the two will not be considered visible. 
%   PrecnNutnExpire     days                   For how long, in days, a
%                                              computed value for
%                                              precession and nutation of
%                                              the Earth can be used before
%                                              new values need to be
%                                              recomputed.  Reusing values
%                                              improves performance of
%                                              ECI-to-ECEF conversions
%                                              which are done inside the
%                                              GPS computations.  The
%                                              default value is 1 day.
%   AntennaOrientation  3 x 3 x num_ant        unit orthogonal matrices
%                                              representing the rotation of
%                                              each antenna with respect to
%                                              the body.  The Z-axis of the 
%                                              antenna frame is the 
%                                              boresite.  This option
%                                              parameter is used only if
%                                              the corresponding antenna
%                                              pattern is 2-D.  If it's 1-D
%                                              then the antenna boresite is
%                                              computed based on the
%                                              antenna pointing specified
%                                              in AntennaPointing
%                                              parameter. Default is
%                                              identity matrix.
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each gps satellite are grouped together in
% rows. Thus, requesting range and range rate, the output y will look like:
%   y = [   range_sv1(t1)       range_sv1(t2)      ...  range_sv1(tn);
%           rangeRate_sv1(t1)   rangeRate_sv1(t2)  ...  rangeRate_sv1(tn);
%           range_sv2(t1)       range_sv2(t2)      ...  range_sv2(t2);
%           rangeRate_sv2(t1)   rangeRate_sv2(t2)  ...  rangeRate_sv2(t2);
%           range_sv3(t1)       range_sv3(t2)      ...  range_sv3(t2);
%           rangeRate_sv3(t1)   rangeRate_sv3(t2)  ...  rangeRate_sv3(t2)
%           ...                 ...                ...  ...
%           range_sv32(t1)      range_sv32(t2)     ...  range_sv32(t2);
%           rangeRate_sv32(t1)  rangeRate_sv32(t2) ...  rangeRate_sv32(t2)]
%
%   NaN's are output in the y elements where measurements cannot be
%   supported due to geometric or signal constraints.
%
%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      t            (1xN)	measurement times (secs from epoch in options)
%      x            (6xN)   spacecraft state [position;velocity] (km)
%      options      (1x1)   data structure (see above description)
%      qatt         (4xN)   spacecraft body quaternion history where the
%                           quaternion specifies the user spacecraft
%                           attitude with respect to the same coordinate
%                           frame that the position and velocity are given
%                           in.  The quaternion is specified as
%                           [sin(theta/2)*e_vec;cos(theta/2)] where e_vec is the
%                           unit vector representing the axis of rotation
%                           and theta is the angle of rotation about that
%                           axis.  The quaternion information is used IF a 
%                           2-D antenna is specified for the user satellite.  
%                           If the quaternion is not supplied and the 2-D 
%                           antenna is specified, then it will assume that 
%                           spacecraft attitude is aligned with the 
%                           coordinate frame that the position and velocity 
%                           are given in, i.e., the default is [0 0 0 1]'.
%                           Additionally, the orientation matrix
%                           of each antenna wrt to the spacecraft body can 
%                           be specified through the options strucutre.  If
%                           the antennae orientation is not specified, the
%                           default is identity, i.e., aligned with the body
%                           axes.
%
%   OUTPUTS
%      y            (MxN)   measurements (km, and km/s or Hz)
%      H            (Mx6xN) measurement partials matrix
%      R            (MxMxN) measurement covariance
%      AntLB        {num_antx1}  Cell array containing link budget for each
%                           antenna.  Note that "masked" refers to an applied
%                           bias due to visibility or blockage.  The fields
%                           are:
%           Halpha_r    (MxN)   The receiver elevation angle (rad)
%           Hvis_beta   (MxN)   Logical array where both transmit and
%                               receive antenna elevation angles are
%                               within antenna mask angle limits
%           Hvis_CN0    (MxN)   Logical array where CN0 is above the
%                               acquisition/tracking threshold
%           HCN0        (MxN)   Masked Signal carrier to noise ratio
%           HAd         (MxN)   Masked attenuation from R^2 losses (dB)
%           HAr         (MxN)   Receive antenna gain (dB)
%           HAP         (MxN)   Masked budget gain before receiver antenna (dB)
%           HRP         (MxN)   Masked budget gain before receiver amplifiers
%                               and conversion (dB)
%           HAt         (MxN)   Masked transmit antenna gain (dB)
%
%   keyword: measurement
%   See also odtbxOptions, getgpsmeas, gpslinkbudget
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Kevin Berry         09/07/2007      Original
%   Russell Carpenter   07/18/2008      Various bug fixes & improvements
%                             (for other changes, see the svn repository)
%   Sun Hur-Diaz        02/05/2011      Corrected Doppler computation
%   Sun Hur-Diaz        02/10/2011      Implemented 2-D receive antenna usage
%   Kevin Berry         02/28/2011      Corrected Jacobian calculation
%   Sun Hur-Diaz        03/03/2011      Refactored physical paramater
%                                       calculations
%   Allen Brown         03/03/2011      Refactored gps link budget
%   Kevin Berry         03/08/2011      Modified doppler and added range
%                                       rate for consistency with other 
%                                       ODTBX measurement models
%   Phillip Anderson    10/16/2013      Abstracted out link budget code

d2r             = pi/180;

%% Get values from options
useRange        = getOdtbxOptions(options, 'useRange', true );
useRangeRate    = getOdtbxOptions(options, 'useRangeRate', false );
useDoppler      = getOdtbxOptions(options, 'useDoppler', false );
link_budget      = getOdtbxOptions(options, 'linkbudget', []);

%% Link Budget values required from gpsmeas 
% (set them if something isn't already there)
link_budget = linkbudget_default(link_budget, 'GPSBand', 'L1');
link_budget = linkbudget_default(link_budget, 'AntennaPattern', {'sensysmeas_ant.txt','sensysmeas_ant.txt'});
    %  Specify antenna pattern for each antenna, existing antennas are:
    %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
    %     omni.txt                  - zero dB gain,  180 degree half beamwidth
    %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
    %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
    %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
num_ant = length(link_budget.AntennaPattern); %hasn't been tested for >4 antennas
link_budget = linkbudget_default(link_budget, 'RXAntennaMask', 180*d2r);
link_budget = linkbudget_default(link_budget, 'AtmosphereMask', 50); % km
    %  Troposphere mask radius ~50 km
    %  Ionosphere mask radius ~(500-1000 km)
link_budget = linkbudget_default(link_budget, 'NoiseTemp', 300); % K
    % System noise temp [K], space pointing antenna = 290
    % System noise temp [K], earth pointing antenna = 300
link_budget = linkbudget_default(link_budget, 'AtmAttenuation', 0.0); % dB
    % attenuation due to atmosphere (should be negative) [dB]
link_budget = linkbudget_default(link_budget, 'TransPowerLevel', 2); % 1-minimum, 2-typical, 3-max
link_budget = linkbudget_default(link_budget, 'TransPowerOffset', 0.0); % dB, global offset
link_budget = linkbudget_default(link_budget, 'GPSBlock', 1 ); 
    % 1 - II/IIA
    % 2 - IIR
    % 3 - IIR-M
    % 4 - IIF
    % 5 - III  - currently use IIR-M gain pattern and -158.5 min power
    % 6 - Fictitious 3D antenna pattern based on IIA
link_budget = linkbudget_default(link_budget, 'TXAntennaMask', pi );  % in rad
    %  The actual mask used is the lesser of this mask and the limit of the defined pattern
    %  Note:  mask = 70 deg includes entire defined pattern
    %         mask = 42 deg includes only main and first side lobes
    %         mask = 26 deg includes only main lobe
link_budget = linkbudget_default(link_budget, 'ReceiverNoise', -3 );  % dB, Noise figure of receiver/LNA
link_budget = linkbudget_default(link_budget, 'RecConversionLoss', -1.5 );  % dB
    % Receiver implementation, A/D conversion losses [dB]
    %   Novatel: L = -4.0 	
    %   Plessey: L = -1.5		
link_budget = linkbudget_default(link_budget, 'SystemLoss', 0 ); % dB, System losses, in front of LNA
link_budget = linkbudget_default(link_budget, 'LNAGain', 40 ); % dB, LNA gain (trimble pre-amp spec = 42-48)
link_budget = linkbudget_default(link_budget, 'CableLoss', -2 ); % dB, Cable losses (after LNA)
link_budget = linkbudget_default(link_budget, 'RecAcqThresh', 32 ); % dB-Hz, Receiver acquisition threshold
link_budget = linkbudget_default(link_budget, 'RecTrackThresh', link_budget.RecAcqThresh ); % dB-Hz, Receiver tracking threshold
link_budget = linkbudget_default(link_budget, 'DynamicTrackRange', 15 ); % dB
    % Dynamic tracking range of receiver, or maximum difference  
    % in power levels tracked simultaneously. If the difference
    % in snrs between two satellites is more that link_budget.DynamicTrackRange,
    % the weaker of the two will not be considered visible. 

% Reassign the options structure with any changed/default link budget values
options = setOdtbxOptions(options, 'linkbudget', link_budget);
    
% Error messages associated with user inputs
if( useDoppler && useRangeRate)
    error('gpsmeas is not designed to handle both RangeRate and Doppler at the same time')
end
if(length(t) == 1 && link_budget.RecAcqThresh ~= link_budget.RecTrackThresh)
    warning('WarnTests:CNOLIM',['gpsmeas can only utilize a tracking '...
        'limit that is different from the acquisition limit when it is '...
        'run for all times at once (i.e. batch mode).'])
end

% Count number of measurement types
nMTypes         = sum([useRange useRangeRate useDoppler]);


%  Select L-band carrier for simulation
%  Models civil signal power on selected carrier (C/A on L1, L2C on L2)
GPSFreq.L1      = 1575.42e6;    % Hz
GPSFreq.L2      = 1227.6e6;     % Hz
GPSFreq.L5      = 1176.45e6;    % Hz
freq            = GPSFreq.(link_budget.GPSBand);
% link_budget.Frequency = freq;
options = setOdtbxOptions(options, 'frequencyTransmit', freq);
% r_mask          = EARTH_RADIUS + AtmMask;	% Atmosphere mask radius (km)
 

%Set up block type for antenna pattern selection, indexed by PRN
block=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

if nargin < 5
    GPS_SIZE = 32;
else
    params.PRN = sv;
    GPS_SIZE = 1;
    link_budget.GPSBlock = block(sv);
end
                    
%% Input Evaluation
% This Section includes GPS satellite gain pattern file names, reference transmitter
%  power values, and error checking of parameters specified in definitions file.
% Power settings all referenced to meet required minimum power at EOE
% NOTE: GPS antenna pattern data beyond 30 degrees considered LM
%       proprietary for blocks other than II/IIA
% NOTE: Many of the antenna patterns below are 1D (elevation only).

% GPS SV (transmitter) antenna file definitions:
%  The format for 1D (elevation) antenna patttern files is as follows:
% 
%  % Header lines
%  % .
%  % .
%  % .
%  % Header lines
%  % angle [deg] gain [db]
%  0.0  gain1
%  .    .
%  .    .
%  angN gainN
% 
%  The comment character for antenna pattern files is the
%  "%" character and there can be as many header lines as
%  desired. The antenna pattern is specified as two columns
%  of numbers where the first column contains the angle
%  off the antenna boresight in units of degrees and the
%  second column contains the gain value in units of dB.
%  The angular resolution at which these data are
%  specified is up to the user, but the range is usually from 0 deg to 
%  90 deg or less.
%
%  The format for 2D antenna patttern files is as follows:
% 
%  % Header lines
%  % .
%  % .
%  % .
%  % Header lines
%  % One row with a placeholder scalar and then the applicable azimuth 
%  % angles [deg]
%  % Then on subsequent rows: the elevation angle [deg]
%  % deg or less and the gain value [db] for each applicable azimuth angle
%
%  -50  0.0     az2     ...     azM
%  0.0  gain1,1 gain1,2 ...     gain1,M
%  .    .
%  .    .
%  elN  gainN,1 gainN,2 ...     gainN,M
% 
%  The comment character for antenna pattern files is the
%  "%" character and there can be as many header lines as
%  desired. The antenna pattern is specified as a table
%  of numbers where the first column contains the angle
%  off the antenna boresight (elevation) in units of degrees and the
%  subsequent columns contains the gain value in units of dB.
%  The angular resolution at which these data are
%  specified is up to the user, but the elevation range is usually from 0 
%  deg to 90 deg or less.  The azimuth angles should range from 0 to < 360
%  providing complete azimuth coverage.

switch freq
    case GPSFreq.L1
        switch link_budget.GPSBlock
            case 1  % II/IIA
                xmit_pattern = 'GPSIIA_L1MEAN.txt';	 % GPS SV gain pattern file
                xmit_power = 13.4;
            case 2  % IIR
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                xmit_pattern = 'GPSIIR_ALL_L1MEAN_30deg.txt';	 % GPS SV gain pattern file
                %xmit_pattern = 'GPSIIR_ALL_L1MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 13.5;
            case 3  % IIR-M
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                xmit_pattern = 'GPSIIRM_P1_L1MEAN_30deg.txt';	 % GPS SV gain pattern file
                %xmit_pattern = 'GPSIIRM_P1_L1MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 12.8;
            case 4  % IIF
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                error('Unsupported Option: No IIF antenna patterns defined');
            case 5  % III  - currently use IIR-M gain pattern and -158.5 min power
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                %xmit_pattern = 'GPSIIRM_P1_L1MEAN.txt';	 % GPS SV gain pattern file
                xmit_pattern = 'GPSIIRM_P1_L1MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 12.8;
            case 6  % Ficticious 3D antenna pattern based on IIA
                xmit_pattern = 'GPSIIA_L1_3Dantenna.txt';	 % GPS SV gain pattern file
                xmit_power = 13.4;
        end
    case GPSFreq.L2
        switch link_budget.GPSBlock
            case 1  % II/IIA
                error('Unsupported Option: No IIA L2 antenna pattern defined');	 % GPS SV gain pattern file
            case 2  % IIR
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                xmit_pattern = 'GPSIIR_ALL_L2MEAN_30deg.txt';	 % GPS SV gain pattern file
                %xmit_pattern = 'GPSIIR_ALL_L2MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 6.9;
            case 3  % IIR-M
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                xmit_pattern = 'GPSIIRM_P1_L2MEAN_30deg.txt';	 % GPS SV gain pattern file
                %xmit_pattern = 'GPSIIRM_P1_L2MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 9.4;
            case 4  % IIF
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                error('Unsupported Option: No IIF antenna patterns defined')
            case 5  % III  - currently use IIR-M gain pattern and -158.5 min power
                warning('ODTBX:gpsmeas:gpsblock', 'Only GPS Block II/IIA is publicly available.');
                xmit_pattern = 'GPSIIRM_P1_L2MEAN_30deg.txt';	 % GPS SV gain pattern file
                xmit_power = 10.9;
            case 6  % Ficticious 3D antenna pattern based on IIA
                error('Unsupported Option: No 2D antenna patterns defined for L2 frequency.');
        end
end

% Apply global offset to transmitter power
xmit_power = xmit_power + link_budget.TransPowerOffset;

switch link_budget.TransPowerLevel    % Set power into transmitter
    case 1
        P_sv = xmit_power;
    case 2
        P_sv = xmit_power + 1.5;
    case 3
        P_sv = xmit_power + 3.0;
end

nn = length(t);
% Set receiver antenna loop number
loop = max([1,num_ant]);

%% Determine the dimension of the antenna pattern data
TXpattern = load(xmit_pattern);  % [x,2]
if size(TXpattern,2) > 2
    xmit_pattern_dim = 2;
else
    xmit_pattern_dim = 1;
end
rec_pattern_dim = ones(loop,1);
RXpattern = cell(loop,1);
for ANT = 1:loop
    RXpattern{ANT} = load(link_budget.AntennaPattern{ANT});  % [x,2]
    if size(RXpattern{ANT},2) > 2
        rec_pattern_dim(ANT) = 2;
    end
end

%% Compute the gps physical relationships between the transmitters and receivers
params.xmit_pattern_dim = xmit_pattern_dim;
params.rec_pattern_dim = rec_pattern_dim;
params.GPS_SIZE = GPS_SIZE;
params.num_ant = num_ant;
if nargout > 1
    params.doH = 1;
else
    params.doH = 0;
end

if nargin < 4
    qatt = [];
end

% Calculate the physical parameters
out = getgpsmeas(t,x,options,qatt,params);

%% Generate link budget structures from calculated data

% set link budget params in structs
TX_link.P_sv = P_sv;

% Transmitter and Receiver antenna patterns
RX_link.pattern = RXpattern;
TX_link.pattern = TXpattern;

%% Calculate Link Budget
[AntLB, HVIS] = calc_linkbudgets(out, options, RX_link, TX_link);

%% Pull information from measurements helpful in calculating outputs
% the physical parameter results:
epoch = out.epoch;       % epoch of first time [1x1]
% TX_az = out.TX_az*d2r;   % The transmitter azimuth angle (rad) [nn x GPS_SIZE]
% TX_link.alpha = out.TX_el*d2r; % The transmitter elevation angle (rad) [nn x GPS_SIZE]
% RX_link.alpha = out.RX_el*d2r; % The receiver elevation angle (rad)

% Note, the receiver angles from "out" are below, in the ANT loop
Hrange = out.range;      % [GPS_SIZE x nn]
Hrrate = out.rrate;      % [GPS_SIZE x nn]
% rgps_mag = out.rgps_mag; % [nn x GPS_SIZE]
% health = out.health;     % the health indicator, [nn x GPS_SIZE]
% dtsv = out.dtsv;          %individual satellite clock bias which is added to the epoch time to reflect the GPS time of measurement (using AFO and AF1)

% Auxiliary parameters useful for H calculation
if params.doH == 1
    eciRotation = out.eciRotation;  % Used in measurement Jacobian
    los_3d = out.los_3d;
    gps_vel_tot = out.gps_vel_tot;
    sat_vel_tot = out.sat_vel_tot;
    Rotation2ECI = out.Rotation2ECI;
end

%% Create Outputs
y = NaN(size(Hrange,1)*(nMTypes),size(Hrange,2));
H = NaN(size(Hrange,1)*(nMTypes)*(nargout>1),6,size(Hrange,2));
nancheck=[NaN 0];
ind_y = 0;
ind_H = 0;
for n=1:size(Hrange,1)

    if( useRange )
        ind_y = ind_y + 1;
        y(ind_y,:) = Hrange(n,:)+nancheck(HVIS(n,:)+1);
    end
    if( useRangeRate )
        ind_y = ind_y + 1;
        y(ind_y,:) = Hrrate(n,:)+nancheck(HVIS(n,:)+1);
    end
    if( useDoppler )
        ind_y = ind_y + 1;
        y(ind_y,:) = -(freq/c)*Hrrate(n,:)+nancheck(HVIS(n,:)+1);
    end

    if nargout > 1,
        dtsv = zeros(3,nn);
        dr = zeros(3,nn);
        dv = zeros(3,nn);
        for nt=1:nn
            R2ECI    = Rotation2ECI(t(nt)+epoch);
            M        = R2ECI(1:3,1:3)' * eciRotation(1:3,1:3,nt)';
            dr(:,nt) = M * los_3d(:,nt,n); %pos diff vector
            dv(:,nt) = M * (gps_vel_tot(:,nt,n) - sat_vel_tot(:,nt)); %vel diff vector
        end
 
        udr = unit(dr);
        rho = sqrt(sum(dr.^2));
        Hvv = -shiftdim(udr,-1);
        Hrr = -shiftdim(udr,-1);
        Hvr = -shiftdim(dv.*repmat(1./rho,3,1) ...
            - dr.*repmat(dot(dr,dv)./rho.^3,3,1),-1);
        Hrv = zeros(size(Hvr));
        if( useRange )
            ind_H = ind_H +1;
            H(ind_H,:,:) = [Hrr Hrv];
        end
        if( useRangeRate )
            ind_H = ind_H +1;
            H(ind_H,:,:) = [Hvr Hvv];
        end
        if( useDoppler )
            ind_H = ind_H +1;
            H(ind_H,:,:) = -(freq/c)*[Hvr Hvv];
        end
    end
end

%% Set the measurement covariance output
if nargout > 2,
    %get sigma out of the options
    sigmaDefault = ones(1,size(y,1))*1e-3;
    sigma        = getOdtbxOptions(options, 'rSigma', sigmaDefault );
    if( length(sigma)~=length(sigmaDefault) )
        disp('WARNING: In gpsmeas, length(rSigma) does not match total number of measurements');
        disp('   Setting rSigma = default (1e-3) for all measurements');
        sigma = sigmaDefault;
    end
    sigma           = diag(sigma);
    R = repmat(sigma.^2,[1,1,length(t)]);
end