% Demo for GPS link budget and antenna pattern analysis tools.
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

echo on

% ----------------------------------
% GPS Enhancements Tools demo script
% ----------------------------------
%
% This script demonstrates the use of the GPS Enhancements Toolset tools 
% in ODTBX to investigate GPS link budget parameters and antenna patterns.  
% These tools are: gps_phys_param, gps_est_cno, and gps_gain.
%
% This demo will walk a user through some sample uses of these three tools
% and the accompanying data filters to take GPS measurement data and
% analyze transmitter and receiver antenna patterns.
%
% Note, the measurement data in this demo is produced using the ODTBX 
% gpsmeas function whereas the real input data would be from a GPS receiver.
%
pause % Hit RET to begin

% First, we will set up some constants for later...
r2d=180/pi; % rad to deg
d2r=pi/180; % deg to rad

% Next, we specify the GPS almanac file that we wish to use.  The almanac 
% should cover the time period for the GPS measurement data we wish to analyze.
yuma_file = 'Yuma1134_example.txt';

% --------------------------------
% Artificial Data Generation Begin
% --------------------------------
% For this demo, since we don't have direct GPS receiver measurements,
% we will create a set of GPS measurements by using a support script and
% ODTBX tools.  We are modeling a receiver onboard an orbiting spacecraft
% with a zenith-pointing antenna for seven hours.
%
% (Please ignore this part of this demo since normal analyst use is to read
% in GPS measurement data from the GPS receiver and all of this 
% processing would not be necessary.)  This may take some time to run.
%
pause % Hit RET to start artificial data generation

echo off;
gps_antenna_tools_support;
echo on;

% Now we have finished creating the GPS measurement data.  Usually this is
% read in from a RINEX Observation file using the rinexo2gpsdata function.
% That function usually creates some metadata for traceability, so we will
% assign that metadata manually right now:
%
% Each receiver/antenna combination should have its own identifier.  This
% is very helpful when dealing wih large datasets and where receivers
% and/or antennas may change.  Let us define this receiver/antenna as #1138.
gps_meas.RX_meta.RX_ID              = 1138;

% The rinexo2gpsdata function records the filename the data came from, we
% don't have one so we'll choose:
gps_meas.RX_meta.meas_file          = '_synthetic_no_file';

% And usually RINEX Observation file header lines are also recorded, we
% don't have those either so we'll just use:
gps_meas.RX_meta.obs_metadata{1}    = 'Synthetic measurement data from GPS Enhancements Tools demo.';

% --------------------------------
% Artificial Data Generation End
% --------------------------------
% Now we have the GPS measurements in the format that the rinexo2gpsdata
% function creates.  Remember that the steps so far have been unique to this demo.
%
pause % Hit RET to see a description of this GPS measurement data

% Before proceeding further, let us look at the measurement data.  The GPS
% measurements are in the structure gps_meas.  It contains the receiver 
% metadata, the list of GPS PRNs that were observed, and a cell array 
% holding a set of data for each PRN.  The GPS_PRN value is the PRN for
% the PRN_data of the same index.  We have measurements for 27 satellites
% over this seven hour period.
gps_meas

% The receiver metadata we set above, but if this was read in from a file
% we could look at it directly:
gps_meas.RX_meta

pause % Hit RET to continue

% Finally, each PRN has the GPS measurements from the receiver.  For PRN
% #1, they are:
gps_meas.PRN_data{1}

pause % Hit RET to see some plots of the GPS measurement data
plot_gps_prange(gps_meas, [31 30 1]); % Plot pseudorange
plot_gps_cno(gps_meas, [31 30 1]);    % Plot C/No

pause % Hit RET to start using the GPS Enhancement Tools

% --------------------------------------------------------------
% GPS Enhancements Tools Step 1: Pre-Process Physical Parameters
% --------------------------------------------------------------
% Since the receiver measurements contain time tags, one of the first
% things that must be done is to create the physical relationship between
% the transmitter and receiver for each measurement.  This physical 
% relationship includes the relative range and range rate between the 
% transmitter and receiver as well as the azimuth and elevation angles for
% both the receiver and transmitter antennas.  In addition, the modeled yaw
% angle of the GPS SV transmitter may be of interest and is included.  Once
% calculated these physical measurements shouldn't change.  They are
% usually stored for subsequent use because they are the foundation
% for link budget analysis.
%
pause % Hit RET to continue

% GPS PRN assignments can change but it is the unique vehicle and
% transmitter that we want to analyze.  We have to keep track of the
% current PRN against a particular vehicle.  Let's say we're only
% interested in eight of the GPS satellites we've observed.  Here is the
% identifying data for our eight satellites.  The ID is our own unique tag
% for managing the data.
TX_ID{1} = makeGpsTXID(345, 1, 1); % ID=345, PRN=1, block II/IIA
TX_ID{2} = makeGpsTXID(346, 2, 1); % ID=346, PRN=2, block II/IIA
TX_ID{3} = makeGpsTXID(347, 3, 1); % ID=347, PRN=3, block II/IIA
TX_ID{4} = makeGpsTXID(348, 4, 1); % ID=348, PRN=4, block II/IIA
TX_ID{5} = makeGpsTXID(349, 5, 2); % ID=349, PRN=5, block IIR
TX_ID{6} = makeGpsTXID(350, 6, 3); % ID=350, PRN=6, block IIR-M
TX_ID{7} = makeGpsTXID(351, 7, 4); % ID=351, PRN=7, block IIF
TX_ID{8} = makeGpsTXID(352, 8, 5); % ID=352, PRN=8, block III
%
pause % Hit RET to continue

% Currently we have measurement data for 27 satellites but we are
% interested in only eight satellites, so let's reduce the measurement dataset.
% We can filter by GPS block type.  Let's create a filter to keep only 
% Block II/IIA:
block_filter = FilterGpsBlock('include',TX_ID,1);

% In fact, let's focus only on one of our satellites, PRN #1 for this demo.
% We will create another filter to further reduce the measurement dataset
% and speed up processing.  (Note that other filter types are available,
% including filtering based on a window of time.)
prn_filter = FilterGpsPrn('include',1);

% While we could use only one filter to accomplish selecting PRN #1 data, 
% let us demonstrate how to combine multiple filters together so we can use
% them when processing physical parameters:
combinedfilter = FilterGpsMulti(block_filter, prn_filter);
%
pause % Hit RET to continue

% As we bring this data together we can note the receiver state source
% information for future reference.  We will use this script name instead
% of a data file.
RX_state_source = 'gps_antenna_tools_demo.m';

% In order to create this physical relationship data we will need to relate
% the measurements to the position, velocity, and orientation of the
% receiver and the GPS almanac.  This occurs in the gps_phys_params
% function.  
%
% The gps_meas are our GPS measurements.  The 1138 is the
% receiver/antenna identifier that we use to ensure our measurement data
% matches the receiver states (the epoch, t, x, and qatt, all in ECI).
% The dcm argument specifies the antenna-to-(spacecraft) body 
% orientation.  The yuma_file is the GPS almanac and the TX_ID combines our
% knowledge of the GPS satellites with the measurements.  Finally we use
% our filters to reduce the data set down to one satellite.
%
phys_param = gps_phys_params(gps_meas, ...
    RX_state_source, 1138, epoch, t, x, qatt, [], ...
    dcm('ax2',pi/2), ...
    yuma_file, TX_ID, combinedfilter);
%
pause % Hit RET to see the results

% The physical parameter results are a cell array for each GPS PRN.  We
% only have one PRN of interest from our filter and this is what it
% contains:
phys_param{1}
%
pause % Hit RET to see the metadata

% This data can be saved for future use and/or kept in MATLAB for further
% processing.  It contains metatdata about the receiver, transmitter, and
% the receiver state source so that the analyst can re-use it in the future
% and understand its pedigree.  The metadata is:
phys_param{1}.meta
phys_param{1}.meta.RX_meta
phys_param{1}.meta.TX_ID
%
pause % Hit RET to see some plots of these physical parameters
plot_gps_pp_rrdot(phys_param{1});    % Plot range and range rate
plot_gps_pp_azel(phys_param{1}, 1);  % Plot receive azimuth and elevation angles
plot_gps_pp_azel(phys_param{1}, 2);  % Plot transmit azimuth and elevation angles
plot_gps_pp_gpsyaw(phys_param{1});   % Plot GPS satellite yaw angle
%
pause % Hit RET to continue


% --------------------------------------------------------------
% GPS Enhancements Tools Step 2: Pre-Process Link Estimates
% --------------------------------------------------------------
% Now that we have some physical parameters that support GPS link budget
% analysis, we could calculate several link budget variables.  In our next
% step, let us calculate an estimated C/No.  We might use this to perform
% gross adjustments on our link budget parameters by comparing the
% estimated C/No to a measured value.
%
pause % Hit RET to continue

% The link budget is comprised of several items for the transmitter and the
% receiver.  Each antenna has a set of pattern data that specifies a 2D
% antenna gain model based on boresight azimuth and elevation.  (Here we
% will set the data to the antenna patterns used during our measurement
% generation so that they match well.)
%
rx_pat = load(truth.rec_pattern);
tx_pat = load('GPSIIA_L1_3Dantenna.txt');

% The transmitter's power is kept consistent with the PRN (based on block
% type, etc.)
TX_link.P_sv = 14.9; % dB

% The receiver has several values of interest. (These values are being set 
% to coincide with the measurement generation but an analyst would use 
% values obtained from a variety of sources and estimates.)
RX_link.Nf = truth.Nf;      % Noise figure of receiver/LNA [dB]
RX_link.L = truth.L;        % Receiver A/D conversion losses [dB]
RX_link.freq = 1575.42e6;   % GPS L1 frequency, [Hz]
RX_link.Ts = truth.Ts;      % Receiver system noise temp [K]
RX_link.As = truth.As;      % Receiver system losses, in front of LNA [dB]
RX_link.Ae = truth.Ae;      % attenuation due to atmosphere [dB]

% We put these together to estimate the C/No with the physical parameters.
% Note, additional filters could be used here to focus the dataset for this
% analysis.  This tool can be run many times with different link budget
% parameters using the same physical parameters.
CN0_est = gps_est_cno(phys_param{1}, RX_link, TX_link, rx_pat, tx_pat);
%
pause % Hit RET to plot results
plot_est_cno(CN0_est, 1);   % Plot estimated C/No
%
pause % Hit RET to continue


% --------------------------------------------------------------
% GPS Enhancements Tools Step 3: Process GPS Gain Data
% --------------------------------------------------------------

% So far we have the physical parameters that support these measurements
% and we have an idea of the values for our link budget parameters.  
%
% Now we will turn the link budget equation around and compute the gain 
% values for one of our antennas using the GPS measurements.  For example,
% if we knew one antenna well, either the transmit or receive antenna, but 
% the other antenna was poorly characterized, we can use the raw GPS signal
% strength measurement from our receiver to estimate the gain for the 
% unknown antenna as a function of azimuth and elevation angle.
%
% Let us use the gps_gain tool to demonstrate reconstructing antenna gain
% values both ways.  We will use the raw GPS signal strength from our 
% receiver measurements, the corresponding physical parameters for those 
% measurements, the link budget parameters, and the "known" antenna gain 
% pattern.
% Case 1) Use the "known" receiver antenna gain pattern to calculate the
%         "unknown" transmitter antenna gain measurements.
% Case 2) Alternately, use the "known" transmitter antenna gain pattern to 
%         calculate the "unknown" receiver antenna gain measurements.
%
pause % Hit RET to compute

% Call gps_gain for the transmitter gain values:
tx_gain = gps_gain(gps_meas, phys_param{1}, RX_link, TX_link, rx_pat, 1);

% Call gps_gain for the receiver gain values:
rx_gain = gps_gain(gps_meas, phys_param{1}, RX_link, TX_link, tx_pat, 2);

% Now that the calculations are done, let us see how we did.
%
% Given that we generated our own measurements and used the same link 
% budget and antenna pattern data for our tools, we should have a perfect
% match between the calculated gain values and the actual antenna patterns
% that we used.  Let us co-plot the two actual antenna patterns with our
% calculated gains.
%
pause % Hit RET to plot

% plot the calculated receive gain values
plot_gps_gain_2d(1, rx_gain, 1);

% plot the true receive rx_pattern
plot_gps_gain_3d(1, rx_gain, 1, rx_pat);

% plot the calculated transmit gain values
plot_gps_gain_2d(2, tx_gain, 1);

% plot the true transmit rx_pattern
plot_gps_gain_3d(2, tx_gain, 1, tx_pat);


% ---------------------------------------
% GPS Enhancements Tools demo script end
% ---------------------------------------
% This concludes this demo.
echo off;