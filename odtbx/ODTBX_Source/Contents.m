% ORBIT_DETERMINATION_TOOLBOX
%
% Functions
%   alm2xyz                  - Compute R and V from GPS alamanc data
%   alm_pos_az_el            - Calculates GPS SV antenna azimuth and elevation for a given receiver state.
%   antennaCalc              - 
%   antennaVis               - 
%   Attitude                 - Class to store spacecraft attitude information
%   AutoDX                   - Compute the step size dX that minimizes the error in dF/dX,
%   Body                     - Class representing a simple solar system body (i.e. planet, asteroid, comet, etc.)
%   bplane                   - Function to calculate B-Plane parameters for a hyperbolic orbit
%   calc_linkbudgets         - Calculates link budget and visibility between satellite targets
%   calc_obliquity           - Compute Earth obliquity of the ecliptic
%   calc_JulianCenturies2000 - Compute the time in Julian Centuries since the ephemeris epoch of 2000
%   cameradyn                - Models dynamics for camera calibration parameters
%   Camera                   - Class to model an onboard imaging device used for optical navigation
%   CCDCamera                - Extends Camera to include CCD pixel mapping and distortion effects.
%   chi2cdf_odtbx            - ODTBX version of chi2cdf
%   chi2inv_odtbx            - ODTBX version of chi2inv
%   chi2pdf_odtbx            - ODTBX version of chi2pdf
%   chkestopts               - Check some of the estimator options for consistency
%   ChPartModel_basic        - Calculate the signal delay (in km) caused by a signal traveling 
%   clkdyn2                  - A 2-state clock error dynamics model.
%   clkdyn3                  - A 3-state clock error dynamics model.
%   clkprop2                 - A 2-state clock error propagator.
%   clkprop3                 - A 3-state clock error propagagor.
%   compValsTol              - Compares values using a tolerance.
%   comp_bs_3d               - Computes unit vectors describing antenna boresite 
%   copyGpsMetaData          - Copies the metadata from a GPS data structure into a blank GPS data structure.
%   corrmat                  - Correlation Matrix.
%   covmake                  - Covariance maker.
%   covprop                  - Discrete propagation of covariance matrix.
%   covsmpl                  - Covariance sampling function.
%   dataCache                - Creates and manipulates a data cache.
%   dcm                      - Direction Cosine Matrix.
%   dcm2q                    - Converts the direction cosine matrix to quaternion
%   ddhermite                - interpolates based on modified divided difference.
%   ddormeas                 - Makes delta-differenced oneway range measurements.
%   EarthOrbitPlot           - Creates a 3-D plot of an Earth orbit trajectory
%   editmeas                 - Measurement editing - scalar or vector
%   errellip                 - 2-D and 3-D Error Ellipse plotter.
%   estbat                   - Batch Estimator.
%   estbatCon                - Batch Estimator.
%   estpar                   - Particle Filter Estimator.
%   estseq                   - Sequential Estimator.
%   estspf                   - Sigma-Point Filter. 
%   estsrif                  - Square-Root Information Filter (SRIF)
%   estval                   - Evaluation of estimator outputs (plotting).
%   EOPLib                   - This class facilitates the CIO based celestial to terrestrial transformation
%   IERSLib                  - This class facilitates the CIO based celestial to terrestrial transformation.
%   fgprop                   - Lagrangian two-body propagation.
%   fixzoom                  - Fixes the zoom function so that zero stays in view when vertically zooming.
%   FilterEstCnoThresh       - Filters GPS data structs against an estimated C/No threshold.
%   FilterGpsBlock           - A FilterGpsData class for filtering GPS measurement data by GPS block
%   FilterGpsCnoThresh       - A FilterGpsData class for filtering GPS measurement data via a C/No
%   FilterGpsData            - An abstract class for filtering GPS measurement data from makeGpsData.
%   FilterGpsExplicitTimes   - A FilterGpsData class for filtering GPS measurement data by selecting
%   FilterGpsMulti           - A FilterGpsData class that uses multiple FilterGpsData filters on GPS measurement data.
%   FilterGpsPrn             - A FilterGpsData class for filtering GPS measurement data by PRN
%   FilterGpsTimeWindow      - A FilterGpsData class for filtering GPS measurement data by selecting
%   FilterPpData             - An abstract class for filtering GPS physical parameter data from gps_phys_params.
%   FilterPpEl               - A FilterPpData class for filtering GPS physical parameter data on antenna
%   FilterPpExplicitTimes    - A FilterPpData class for filtering GPS physical parameter by time.
%   gainpenalty_mask         - Check and apply additional mask angle penalty.
%   gasdyn                   - Simulates comet outgassing forces
%   get_linkbudget           - Calculate a link budget between the receiver x1 with num_ant antennas and
%   getgpsmeas               - Computes physical parameters required for GPS based measurements
%   getBodyShadow            - Determine whether any of the specified celestial bodies is shadowing the spacecraft.
%   getOdtbxOptions          - Get ODTBX Options structure parameter.
%   getPulsarData            - Retrieves pulsar data for use by XNAVMEAS
%   gpslinkbudget            - Calculates the GPS transmitter-receiver link budget.
%   gpsmeas                  - Makes GPS based measurements using dynamic C/No tracking threshold
%   gpsmeasAbs               - function [y,H,R,AntLB] = gpsmeasAbs(varargin)
%   gpspseudomeas            - Creates GPS pseudo-measurements (one-way measurements) biased by onboard clock errors.
%   gps_est_cno              - Calculates the estimated carrier to noise ratio from physical parameters and link budget data.
%   gps_gain                 - Calculates antenna gains from GPS physical parameters, CN0, and link budget data.
%   gps_phys_params          - Creates GPS physical parameter data sets from GPS measurements.
%   gsCoverage               - Checks gsmeas for measurement availability and creates a possible tracking schedule.
%   gsmeas                   - Makes ground station based measurements.
%   gsmeasCon                - Makes ground station based measurements.
%   HermiteInterpolator      - Interpolates based on modified divided difference.
%   impulsiveBurn            - Impulsive Maneuver (for RESTARTRECORD).
%   initialize_cov           - Calculates an initial covariance matrix.
%   initialize_nbody         - Initializes and validates inputs for NBODYPM
%   integ                    - Integration of state, and possibly STM and process noise.
%   integev                  - Integration of state, and possibly STM, etc. w/event handling.
%   IonoModel_basic          - Calculate the signal delay (in km) caused by a signal traveling 
%   kalmup                   - Kalman update with sigma-edit.
%   kep2cart                 - Convert Keplerian orbital elements to cartesian states.
%   kepanom                  - Solves for one angle of anomaly given another using Kepler's Equation
%   kepel                    - Compute Keplerian orbital elements from Cartesian states.
%   kepprop2b                - Propagates Keplerian elements using 2-body dymanics
%   lightTimeCorrection      - Calculates light time corrections including optional environmental delays
%   lincov_kf                - Linear covariance analysis for Kalman Filter.
%   linkbudget               - Calculates the transmitter-receiver link budget.
%   linkbudget_default       - Sets a required linkbudget value to a default if
%   lnrmeas                  - Makes Lunar Relay based measurements.
%   LOSDoppler               - Line of sight instantaneous doppler between two objects
%   LOSRange                 - Line of sight range measurement between two objects
%   LOSRangeCon              - Line of sight range measurement between two objects
%   LOSRangeRate             - Line of sight range rate measurement between two objects 
%   LOSRangeRateCon          - Line of sight range rate measurement between two objects 
%   LOSUnit                  - Line of sight unit vector between two objects (provides angle information to the estimator)
%   makeGpsData              - Creates a blank GPS data structure for holding raw GPS measurements.
%   makeGpsTXID              - Creates a blank data structure for describing a GPS transmitter.
%   makePpData               - Creates a blank data structure for holding GPS physical parameters.
%   nbodypm                  - N-Body Point Mass dynamics equations.
%   observ                   - Observability grammian, and possibly mapped innovations.
%   odtbxOptions             - Returns a valid ODTBX options structure of the input type 
%   ominusc                  - Innovations, meas. partials, and possibly innovation covariance.
%   ominuscCon               - Innovations, meas. partials, and possibly innovation covariance.
%   opnavmeas                - Simulates measurements from optical landmark tracking.
%   opt_sched                - Calculates the sensitivities of a measurement schedule.
%   plot_est_cno             - Plots the estimated Carrier to Noise ratio from gps_est_cno vs time.
%   plot_gps_cno             - Plots the raw GPS signal to noise ratio vs time.
%   plot_gps_gain_2d         - Plots the gain from gps_gain vs azimuth and elevation angles (2D plots).
%   plot_gps_gain_3d         - Plots the gain from gps_gain and/or antenna gain pattern vs az and el angles (3D plots).
%   plot_gps_pp_azel         - Plots antenna az and el angle physical parameters from gps_phys_params vs time.
%   plot_gps_pp_gpsyaw       - Plots the modeled GPS yaw angle physical parameter from gps_phys_params vs time.
%   plot_gps_pp_rrdot        - Plots the range and range rate physical parameters from gps_phys_params vs time.
%   plot_gps_prange          - Plots the raw GPS pseudorange vs time.
%   plot_ominusc             - Plots the observed measurement minus the computed
%   polydyn                  - Force model for exterior gravitation of a polyhedron.
%   q2dcm                    - Converts from a quaternion [vector;scalar] to a direction cosine matrix
%   read_rnxn                - Read rinex formated GPS broadcast ephemeris data
%   read_spephem             - Read in SP Ephemeris formatted ASCII text file
%   read_stkephem            - Read in STK Ephemeris (formatted ASCII text file)
%   read_yuma                - Read in yuma format GPS almanac file
%   relativityLightDelay     - Returns light time delay error from relativistic effects
%   rinexo2gpsdata           - Read RINEX 2.x formatted GPS observation data into a makeGpsData struct.
%   RinexOReader             - Reads a RINEX 2.x formatted GPS observation data file.
%   rnxnEval                 - Evaluates ICD-200 GPS Ephemeris at user specified time.  Uses ephemerides
%   rotate_results           - 3D rotation for estimator output
%   rotransf                 - Rotational Transformation of Position, Velocity, & Acceleration.
%   rrdotang                 - Calculate range, range rate, doppler, angles between two objects.
%   rrdotangCon              - Calculate range, range rate, doppler, angles between two objects.
%   rrdotlt                  - Calls rrdotang with the light time correction
%   scrunch                  - "Vectorize" a square symmetric matrix.
%   sensmos                  - Sensitivity Mosaic Plot: A Checkerboard plot of the sensitivity matrix.
%   setOdtbxOptions          - Create/alter OPTIONS structure.
%   sigextract               - Extract standard deviations from "scrunched" covariances.
%   slerp                    - Interpolates between two quaterntions (uniform angular velocity, fixed axis).
%   solve_consider           - solve_consider Summary of this class goes here
%   srpAccel                 - Compute solar radiation pressure acceleration
%   staEarthTide             - Ground station position adjustment for solid earth tides
%   stationData              - Provides database information about selected ground stations
%   sun_vec                  - Computes unit vector from origin of ECI frame to sun
%   tdrssmeas                - Makes TDRSS based measurements.
%   tdrssmeasCon             - Makes TDRSS based measurements.
%   tdrssmeas_basic          - Makes TDRSS based measurements.
%   testMeasPartials         - compares the numerically estimated H matrix to the H matrix provided by a measurement model.
%   TropoModel_basic         - Calculate the signal delay (in km) caused by a signal traveling 
%   ttdelay                  - Compute the station time tag delay error
%   unit                     - Efficiently compute unit vectors and norms for vector arrays.
%   unscrunch                - Recreate a matrix that has been vectorized with SCRUNCH.
%   validateOdtbxOptions     - Checks to see if the given structure is valid.
%   varpiles                 - Variance Sandpiles.
%   visibility_constraints   - 
%   xmat                     - Cross-product matrix for arrays of vectors.
%   xnavmeas                 - Calculates spacecraft to distant RSO range difference & range difference rate measurements. 
