function options = odtbxOptions(type)

% ODTBXOPTIONS  Returns a valid ODTBX options structure of the input type 
%   with empty [] default values for all fields. With no output specified, 
%   the list of valid options for the specified type is displayed.
%
%   OPTIONS = odtbxOptions(TYPE) returns a valid JAT options structure
%   for the input type. Valid types are listed below. If no type is
%   entered, a structure is returned that has all valid options. If no
%   return value is requested, a list of all valid options for that input
%   type is displayed.
%
%   Valid input types are:
%       'estimator'
%       'forceModels'
%       'measurement'
%
%   New options can easily be entered by adding the new option name to the
%   appropriate list. To add a new type of options, the name must be added
%   to the 'types' list, a new options structure must be created, and this
%   new structure must be added to the allOptions cell array.
%
%   The setOdtbxOptions and getOdtbxOptions functions only require an 
%   adequate number of input characters to uniquely identify an option. 
%   However, the code will run more efficiently if the whole option name
%   is input exactly.
%
%   Note, although the optionsType field exists in the ODTBX options
%   structure, it is for internal use only.
%
% ESTIMATOR PROPERTIES
% --------------------
% The estimator options structure contains the properties used by the OD
% Toolbox estimators. This structure was managed using estset and estget in
% Release 1.
%
% OdeSolver - Matlab ODE solver for ODEFUN. [ {ode113} | ode45 | ...]
%   Set this property to whichever Matlab or JAT ODE solver is to be used
%   to integrate ODEFUN between measurement data updates.
%
% OdeSolvOpts - Matlab ODE solver options. [ {empty matrix} ]
%   This property is a container for the ODE solver options structure that
%   may be changed with ODESET.
%
% DatVectorized - Vectorized DATFUN, h(t,x).  [ 1 | {0} ]
%   Set this property 'on' (1) if DATFUN is coded so that h(t,[x1 x2 ...])
%   returns [h(t,x1) h(t,x2) ...].  Note that this property only gets
%   checked if a numerical Jacobian is computed.
%
% DatJTolerance - DATFUN numerical Jacobian threshold. [ {scalar} | vector]
%   Set this property to provide a threshold of significance for the
%   function whose gradient NUMJAC will approximate, i.e. the exact value
%   of a component Y(i) with abs(Y(i)) < DatJTolerance(i) is not important.
%   All components of DatJTolerance must be positive. Default in ominusc() 
%   is 1e-6.
%
% DatJPattern - DATFUN Jacobian sparsity pattern. [ sparse matrix ]
%   Set this property to a sparse matrix S with S(i,j) = 1 if component i
%   of h(t,x) depends on component j of x, and 0 otherwise.
%
% UpdateIterations - Number of update iterations. [ scalar {1} | vector ]
%   Set this property to specify the number of times to iterate the
%   measurement update. Value set must be > 0. If the measurements are 
%   uncorrelated with each other, and the vector's length is the same as 
%   the number of measurements, then a separate iteration for each 
%   measurement may be applied. estseq() & estsrif() defaults to {1} while
%   estbat() defaults to {3}.
%
% MonteCarloCases - Number of Monte Carlo cases to run. [ integer {1} ]
%   Set this property to an integer number of Monte Carlo cases to run.
%   Value must be > 0.
%
% MonteCarloSeed - Random number generator state(s). [ scalar {1} | vector]
%   Set this property to an integer J to reset the random number generator 
%   to its Jth state. See RANDN or RAND for details. If set to a vector of 
%   integers, the random number generator will be reset according to each
%   succeeding element of the matrix at the beginning of each case. This 
%   is especially helpful for re-examining particular subsets of a large 
%   number of Monte Carlo runs, if the random number states at the 
%   beginnings of those cases are known (e.g. from the ESTSOL.case output
%   solution field). If only a single integer is specified but there is 
%   more than one Monte Carlo case, then a list seeds (one per case) will 
%   be created by beginning with the integer provided and incrementing 
%   from it. chkestopts() defaults to {NaN}.
%
% ValidationCase - Output to a file and compare with ODEAS. [ integer {0} ]
%   This will allow the user to output the measurements and trajectory from
%	ODTBX to a file and compare against ODEAS for the validation test case
%   if set to 1. The validation test cases are numbered and based on the 
%   case a different output function is used to compare the values.
%
% SchmidtKalman - Flag for using Schmidt Kalman filter in estseq [{0} | 1 ]
%   Set this to 1 to use the Schmidt Kalman filter instead of the Kalman
%   filter in ESTSEQ.  
%
% UseProcNoise - Flag for using process noise in estbat. [ {1} | 0 ]
%   Set this flag to 0 to ignore process noise. Set this flag to 1 to 
%   include process noise.
%
% EditRatio - Threshold for editing measurements [scalar {9} | vector ]
%   Value that is compared to the quadratic form dy'*inv(Pdy)*dy when
%   deciding whether to accept or reject measurements. For example, for a 
%   "3-sigma edit" requires setting EditRatio equal to 9. If there is more 
%   than one measurement type, a vector of values for the EditRatio may be 
%   supplied where the length of the vector is equal to the number of 
%   measurement types. In this case separate EditRatio values can be 
%   specified for each measurement type. If only one value is provided, it
%   will be used for all measurement types. To disable measurement editing, 
%   set EditFlag = 2.  The default value is 9 corresponding to a 
%   "3-sigma edit."
%
% EditFlag - Flag to set measurement editing. [ scalar {1} | vector ]
%   Setting EditFlag equal to 2 will force a measurement to be processed
%   regardless of EditRatio. Setting EditFlag equal to 1 will cause a 
%   measurement to be accepted if it passes the EditRatio test. In this 
%   case, if the EditRatio test is not passed then the EditFlag is returned
%   as 0 to indicate that the measurement was rejected. EditFlags can be 
%   specified on a per-measurement-type basis in a vector, just like the 
%   EditRatio, or a single value can be provided that will then be used 
%   for all measurement types.
% 
% EditVector - Flag to vectorize measurement editing [ {0} | 1 ]
%   Indicates whether the measurement editing will be performed as (1) a 
%   vector or (0) as a scalar.  The default is a scalar editing.
%
% UpdateVectorized - Flag to vectorize measurement processing. [ {1} | 0 ]
%   Flag determines whether multiple measurements are processed 
%   simultaneously in a vector or individually as scalars. The default
%   behavior is to process multiple measurements in a vector. Valid values 
%   for this flag are:
%
%      0 - Process measurements invidually as scalars.
%      1 - Process measurements together in a vector (default).
%
%   Note: This option is only available in estseq and estsrif.
%
% Refint - Determines the intermediate time steps [ >= 0 {3} | <0 ]+
%    Determines the saving of data between measurement points.
%    >= 0 - Number of intermediate points between measurements. This is a 
%       fixed value for the entire time span. (default = 3)
%    < 0  - Let the variable step integrator, integ, determine the 
%       intermediate steps based on the true reference trajectory 
%       propagation
%
%    Note: Currently this option is implemented in estseq and estsrif only.
%
% UseSmoother - Flag for smoother in sequential estimators. [ 1 | {0} ] 
%   Set this flag to 0 to not use smoother. Set this flag to 1 to use
%   smoother.  Default is for smoother to be off.
%   
%   Note: Currently this option is implemented in estsrif only, but could 
%      also be used in estseq.
%
% FORCEMODEL PROPERTIES
% ---------------------
% The forceModels options structure contains the properties used by the JAT
% force models for orbit integration. These models can be used with either
% Matlab or JAT integrators.
%
% epoch - Simulation start time (UTC) in DATENUM format. [ scalar >=0 {0} ]
%   Set this property to specify the start time of the simulation. It is
%   used to determine planetary and ECI positions for force models.
%
% mass - Spacecraft mass in kg. [ scalar >= 0 {1000} ]
%   Set this property to specify the mass of the spacecraft.
%
% dragArea - Atmospheric drag area in m^2 [ scalar >= 0 {10} ]
%   Set this property to specify the cross-sectional area used in the
%   computation of atmospheric drag.
%
% srpArea - Solar radiation pressure area in m^2. [ scalar >= 0 {10} ]
%   Set this property to specify the cross-sectional area used in the
%   computation of solar radiation pressure forces.
%
% cD - Spacecraft coefficient of drag. [ scalar >= 0 {2.2} ]
%   Set this property to specify the drag coefficient used in the
%   computation of atmospheric drag.
%
% cR - Spacecraft coefficient of reflectivity. [ 0 <= scalar <= 1 {0.7} ]
%   Set this property to specify the reflectivity coefficient used in the
%   computation of solar radiation pressure forces.
%
% earthGravityModel - Earth gravitational model. 
%   [ {'2Body'} | 'JGM2' | 'JGM3' | 'EGM96' | 'GEMT1' | 'WGS84' | ...
%   'WGS84_EGM96' ]
%   Set this property to specify which of the Earth gravity models to use
%   for orbit propagation.
%
% gravDegree - Degree of gravity model. [ scalar > 1 {20} ]
%   Set this flag to specify the degree of the Earth gravitational model.
%
% gravOrder - Order of gravity model. [ scalar > 1 {20} ]
%   Set this flag to specify the order of the Earth gravitational model.
%
% useSolarGravity - Include solar gravity forces. [ true | {false} ]
%   Set this property to 'true' to include solar gravitational forces in
%   orbit propagation.
%
% useLunarGravity - Include lunar gravity forces. [ true | {false} ]
%   Set this property to 'true' to include lunar gravitational forces in
%   orbit propagation.
%
% useSolarRadiationPressure - Include SRP. [ true | {false} ]
%   Set this property to 'true' to include solar radiation pressure forces
%   in orbit propagation.
%
% useAtmopsphericDrag - Include atmospheric drag. [ true | {false} ]
%   Set this property to 'true' to include atmospheric drag forces in
%   orbit propagation.
%
% atmosphereModel - Earth atmosphere model. [ 'NRL' | {'HP'} ]
%   Set this property to specify the Earth atmosphere model to use for drag
%   computations. 
%     'NRL' is the NRL-MISESE-00 model   - valid from ground to space
%     'HP'  is the Harris-Priester model - valid for altitudes between 
%                                          100 km and 2000 km.
%
% nParameterForHPModel - HP model n parameter. [ 2 <= scalar <= 6 {2} ]
%   Set this property to specify the inclination parameter used in the
%   Harris-Priester atmosphere model. The acceptable range of values is 2-6
%   where 2 corresponds to a low-inclination orbit and 6 to a
%   high-inclination orbit.
%
% f107Average - 81-day averaged 10.7cm solar flux. [ scalar >= 0 {150} ]
%   Set this property to specify the 81-day averaged 10.7cm solar flux used in
%   the Harris-Priester and NRL-MISES atmosphere models.
%
% f107Daily - Daily 10.7cm solar flux. [ scalar >= 0 {150} ]
%   Set this property to specify the daily 10.7cm solar flux used in the
%   NRL-MISES atmosphere model.
%
% ap - Daily magnetic flux index. [ scalar >= 0 {4} ]
%   Set this property to specify the daily magnetic flux index used in the
%   NRL-MISES atmosphere model.
%
% CentralBody - Central body. 
%   [ 'Sun' | 'Jupiter Barycenter' | {'Earth'}]
%   Set this property to specify the central body used in the n-body
%   dynamics model. If given a corresponding spice file and GM, primitive 
%   bodies can be used.
% 
% PointMasses - Point masses. 
%   ['Jupiter Barycenter' | 'Earth Barycenter' | {'Sun'} ]
%   Set this property to specify the point masses used in the n-body
%   dynamics model. If given a corresponding spice file and GM, primitive
%   bodies can be used.
% 
% SpiceFiles - Include Spice files. [ 'which('wld18782.15')' | {[]} ]
%   Set this property to specify the spice files to be included in
%   consideration for point masses and/or the central body. Requires a
%   correlating GM value.
% 
% GM - GM for an associated Spice ID [ '2004660' scalar >= 0 | {[]} ]
%   Set this property to specify the GM for the correlating spice ID. The
%   user is required to provide the spice file in SpiceFiles.
% 
% GM_CB - GM value of the central body [ scalar >= 0 | {3.9860e5} ]
%   Set this property to specify the GM value of the correlating central  
%   body as given in CentralBody. 
% 
% GM_PM - GM value(s) of the point mass(es) [ scalar >= 0 | {1.3271e11} ]
%   Set this property to specify the GM value(s) of the correlating point
%   mass(es) as given in PointMasses.
% 
% MEASUREMENT PROPERTIES
% ----------------------
% The measurement options structure contains the properties used by the 
% OD Toolbox measurement models: GSMEAS, LOSRANGE, LOSRANGERATE, LOSDOPPLER
%
% epoch - Simulation start time (UTC) in DATENUM format. [ scalar >= 0 {0}]
%   Set this property to specify the start time of the simulation. It is
%   used to determine ECI positions.
%
% rangeType - Range computation type.  [ 1way | {2way} ]
%   Set this property to '1way' to compute a 1-way range measurement as
%   described in the DATSIM math spec. Setting the property to '2way'
%   computes the 2-way range value.
%
% useRange - Compute range measurements. [ {true} | false ]
%   Set this property to 'true' to calculate range measurements.
%
% useRangeRate - Compute range rate measurements. [ {true} | false ]
%   Set this property to 'true' to calculate range rate measurements.
%
% useDoppler - Compute Doppler measurements. [ true | {false} ]
%   Set this property to 'true' to calculate Doppler measurements.
%
% useUnit - Compute Unit Vector measurements. [ true | {false} ]
%   Set this property to 'true' to calculate unit vectors corresponding to
%   angle based measurements.
%
% useAngles - Compute angle measurements. [ true | {false} ]
%   Set this property to 'true' to calculate two angles corresponding to
%   azimuth and elevation relative to the coordinate system of the state.
%
% useLightTime - Light time correction in measurements. [ true | {false} ]
%   Set this property to 'true' to include the light time delay in
%   measurement calculations. This must be set to true in order to include
%   Ionospheric and Tropospheric effects.
%
% useGPSIonosphere - Ionospheric effects in GPS measurements. [ true | {false} ]
%   Set this property to 'true' to include ionospheric effects in GPS
%   measurement calculations. 
%
% useIonosphere - Ionospheric effects in measurements. [ true | {false} ]
%   Set this property to 'true' to include ionospheric effects in ground
%   based measurement calculations. 
%
% useTroposphere - Tropospheric effects in measurements. [ true | {false} ]
%   Set this property to 'true' to include tropospheric effects in ground
%   based measurement calculations. 
% 
% useChargedParticle - Charged particle effects in measurements. [ true | {false} ]
%   Set this property to 'true' to include charged particle effects in
%   ground based interplanetary measurement calculations.
%
% frequencyTransmit - Transmitter frequency in Hz. [ scalar >0 {1.57542e9}]
%   Set this property to specify the transmitter frequency used for Doppler
%   and ionospheric effects calculations.
%
% gsList - Container for JAT groundStationList object. [JAT groundStationList java object ]
%   Use this property to pass in the JAT groundStationList java object that
%   contains information about all groundstations. This object can be
%   created using CREATEGROUNDSTATIONLIST.
%
% gsID - List of groundstations to use. [ cell array of strings ]
%   Set this property to a cell array of groundstation identifier strings 
%   that specify which groundstations to use in measurement calculations.
%
% gsECEF - ECEF locations of ground stations in kilometers. [m x 3 matrix]
%   Set this property to a matrix of ground station locations to speed up
%   gsmeas, or to allow the use of locations not in the NDOSL.
%
% gsElevationConstraint - Minimum elevation in degrees. [ scalar >= 0 {10}]
%   Set this property to specify the minimum elevation necessary for ground
%   station measurements.
%
% rSigma - Measurement covariance [ scalar >= 0 {1e-3} ]
%   Set this property equal to the diagonal of the measurement noise
%   matrix. Used in computing R in gsmeas.
%
% EarthAtmMaskRadius - Effective Earth Occultation Radius [scalar > 0 {6478.12} km]
%   Set this property to specify the effective radius of Earth Occultation
%   used for measurement processing. 
%
% tdrss - Structure of TDRS data to be used in the tdrssmeas function
%   See tdrssmeas for acceptable formats.
%
% relay - Structure of relay satellite data to be used in the lnrmeas function
%   See lnrmeas for acceptable formats.
%
% ddor - Structure of DDOR data to be used in the ddormeas function
%   See ddormeas for acceptable formats.
%
% xnav - Structure of XNAV data to be used in the xnavmeas function
%   See xnavmeas for acceptable formats.
%
% Schedule - Ground station tracking times [ ID, ti, tf ]
%   Tracking times to be used in the gsmeas or xnavmeas function. See 
%   gsmeas for the data format. In the case of xnavmeas, the pulsar 
%   replaces the ground station. This option can be used to specify 
%   measurements at specific times. For example, if pulsar 1 is measured 
%   at times 0 and 3e5 and pulsar 2 is measuremed at times 1e5 and 4e5, 
%   then the Schedule can be specified with exclusive intervals as follows:
%   Schedule = [1 0        1;
%               1 3e5  3e5+1;
%               2 1e5  1e5+1;
%               2 4e5  4e5+1];
%   with measurement times specified as tspan = [0 1e5 3e5 4e5];
%
%   NOTE: tdrssmeas() expects only [ 'ti' 'tf' ] as a start time string and
%      stop time string. gsmeas() expects the [ID, ti, tf] where the start 
%      and stop times are in seconds from epoch.
%
% PrecnNutnExpire - For how long, in days, a computed value for precession 
%   and nutation of the Earth can be used before new values need to be 
%   recomputed  {1 days}.  
%   Reusing values improves performance of ECI-to-ECEF
%   conversions which are done inside the GPS computations.  
%   The default value is 1 day.
%
% Camera - Class to model an onboard imaging device in OPNAVMEAS
%
% Attitude- Class to store spacecraft attitude information
%
% CalibrationFlag- Flag to indicate if camera calibration parameter partial
%   derivatives should be returned in OPNAVMEAS
%
%   keyword: options
%   See also GETODTBXOPTIONS, SETODTBXOPTIONS, VALIDATEODTBXOPTIONS, 
%   CREATEJATWORLD, CREATEGROUNDSTATIONLIST, GSMEAS, GPSMEAS, LOSRANGE,
%   LOSRANGERATE, LOSDOPPLER, ODESET, ODEGET, INITIALIZE_NBODY
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
%   Derek Surka         07/20/2007      Original JATOptions
%   Derek Surka         08/08/2007      Added atmosphere model parameters
%   Derek Surka         09/07/2007      Renamed and revised
%   Derek Surka         09/18/2007      Added help
%   Kevin Berry         06/05/2008      Seperated Iono from GPS Iono
%   Kevin Berry         09/08/2008      Added useChargeParticle option
%   Kevin Berry         09/28/2008      Added tdrss input
%   Sun Hur-Diaz        01/12/2009      Added Schmidt-Kalman filter option
%   Allen Brown         02/24/2009      Updated comments
%   Allen Brown         02/25-26/2009   Removed UpdateVectorized,
%                                       documented Schedule,
%                                       removed: EditRatio, EditFlag, 
%                                       MonteCarloSeed, MeasurementSeed,
%                                       CovarianceSeed, ValidationFilename,
%                                       MeasOpts, ForceOpts, useGPS,
%                                       useClockError, useMeasurementNoise,
%                                       useAntenna, noise, OutputInterval,
%                                       displayWarnings, dopplerType
%   Kevin Berry         06/15/2009      Added EarthAtmMaskRadius
%   Kevin Berry         06/25/2009      Added time scale comments
%   Kevin Berry         07/01/2009      Added Angle based measurement type
%   Brent Wm. Barbee    08/24/2009      Added documentation for the "use
%                                       process noise" flag (UseProcNoise)
%                                       that I added previously.
%   Kevin Berry         09/03/2009      Changed Angle measurement to Unit
%   Kevin Berry         09/28/2009      Added ddor data structure
%   Kevin Berry         02/02/2010      Added gsECEF parameter for gsmeas
%   Kevin Berry         04/06/2010      Added xnav data structure
%   Sun Hur-Diaz        04/14/2010      Added refint
%   Rosemary Huang      07/07/2010      Added UseSmoother
%   Russell Carpenter   02/10/2011      Added useAngles option
%   Benjamin Asher      05/10/2011      Added: CentralBody, PointMasses,
%                                       SpiceFiles, GM, GM_CB, GM_PM
%   K. Getzandanner     08/15/2011      Added: CalibrationFlag

types = {
    'estimator'
    'forceModels'
    'measurement'
    'all'};  % this must be the last entry in types


estimatorOptions = {
    'OdeSolver'
    'OdeSolvOpts'
    'DatVectorized'
    'DatJTolerance'
    'DatJPattern'
    'UpdateIterations'
    'MonteCarloCases'
    'ValidationCase'
    'SchmidtKalman'
    'UseProcNoise'
    'EditRatio'
    'EditFlag'
    'EditVector'
    'MonteCarloSeed'
    'UpdateVectorized'
    'refint'
    'UseSmoother'
    'Particles'
    };

forceModelOptions = {
    'epoch'
    'mass'
    'dragArea'
    'srpArea'
    'cD'
    'cR'
    'earthGravityModel'
    'gravDegree'
    'gravOrder'
    'useSolarGravity'
    'useLunarGravity'
    'useSolarRadiationPressure'
    'useAtmosphericDrag'
    'atmosphereModel'
    'nParameterForHPModel'      % HP inclination parameter 2-6 (low-high incl)
    'f107Average'               % 81-day averaged 10.7 cm solar flux
    'f107Daily'                 % daily 10.7 cm solar flux
    'ap'                        % daily magnetic flux index
    'CentralBody'
    'PointMasses'
    'SpiceFiles'
    'GM'
    'GM_CB'
    'GM_PM'
    'AttitudeSigma'
    'AttitudeTimeConstant'
    'OpticalBiasSigma'
    'OpticalBiasTimeConstant'
    'isCCD'
    'modelid'};

measurementOptions = {
    'epoch'                     % datenum epoch
    'rangeType'                 % '1way' or '2way'
    'useRange'
    'useRangeRate'
    'useDoppler'
    'useUnit'
    'useAngles'
    'useLightTime'
    'useGPSIonosphere'
    'useIonosphere'
    'useTroposphere'
    'useChargedParticle'
    'frequencyTransmit'         % transmitter frequency
    'gsList'                    % JAT GroundStationList object
    'gsID'                      % cell array of Ground Station ID's to use
    'gsECEF'                    % matrix of ECEF ground station locations
    'gsElevationConstraint'     % groundstation elevation constraint
    'rSigma'
    'EarthAtmMaskRadius'
    'tdrss'
    'relay'
    'ddor'
    'xnav'
    'Schedule'
    'PrecnNutnExpire'
    'OpticalSigma'
    'Camera'
    'Attitude'
    'CalibrationFlag'
    'solvefor'
    'dynamicConsider'
    'localConsider'
    'linkbudget'
    };

allOptions = {
    estimatorOptions
    forceModelOptions
    measurementOptions
    };

if( nargin==0 )
    type = 'all';
end

typeIndex = getIndex(type,types);  % can use types(:,1) for just names

if( nargout==0 )
    switch typeIndex
        case 1
            printEstimatorOptionsHelp;
        case 2
            printForceOptionsHelp;
        case 3
            printMeasOptionsHelp;
        case 4
            printEstimatorOptionsHelp;
            printForceOptionsHelp;
            printMeasOptionsHelp;
    end

else

    % Get the list of field names
    if( typeIndex == length(types) ) % requesting all options
        fieldList = {};
        for j = 1:typeIndex-1
            fieldList = union(fieldList,allOptions{j});
        end
    else
        fieldList = sort(allOptions{typeIndex});
    end

    % Create the fields with empty [] values
    for j = 1:length(fieldList)
        options.(fieldList{j}) = [];
    end
    options.optionsType = types{typeIndex}; % internal field used for input checking

end

end

function printEstimatorOptionsHelp()
fprintf('\n');
fprintf('ESTIMATOR OPTIONS\n');
fprintf('--------------------------\n');
fprintf('                OdeSolver: [ @ode45 | {@ode113} | ... ]\n');
fprintf('              OdeSolvOpts: [ ode options struct | {empty matrix} ]\n');
fprintf('            DatVectorized: [ 1 | {0} ]\n');
fprintf('            DatJTolerance: [ scalar | vector ]\n');
fprintf('              DatJPattern: [ sparse matrix ]\n');
fprintf('         UpdateIterations: [ integer>0 {1} | vector of integers>0 ]\n');
fprintf('          MonteCarloCases: [ scalar>0 {1} ]\n');
fprintf('           MonteCarloSeed: [ NaN | scalar | vector ]\n');
fprintf('           ValidationCase: [ integer {0} ]\n');
fprintf('            SchmidtKalman: [ {0} | 1 ]\n');
fprintf('             UseProcNoise: [ {1} | 0 ]\n');
fprintf('                EditRatio: [ {9} | scalar | vector ]\n');
fprintf('                 EditFlag: [ {1} | 2 | scalar | vector ]\n');
fprintf('               EditVector: [ {0} | 1 ]\n');
fprintf('         UpdateVectorized: [ {1} | 0 ]\n');
fprintf('                   refint: [ >= 0 {3} | < 0]\n');
fprintf('              UseSmoother: [ 1 | {0} ]\n');
fprintf('\n');
end


function printForceOptionsHelp()
fprintf('\n');
fprintf('FORCEMODELS OPTIONS\n');
fprintf('--------------------------\n');
fprintf('                    epoch: [ scalar>=0 {0} ]\n');
fprintf('                     mass: [ scalar>=0 {1000} ]\n');
fprintf('                 dragArea: [ scalar>=0 {10} ]\n');
fprintf('                  srpArea: [ scalar>=0 {10} ]\n');
fprintf('                       cD: [ scalar>=0 {2.2} ]\n');
fprintf('                       cR: [ 0 <= scalar <= 1 {0.7} ]\n');
fprintf('        earthGravityModel: [ {''2Body''} | ''JGM2'' | ''JGM3'' | ''EGM96'' \n');
fprintf('                           | ''GEMT1'' | ''WGS84'' | ''WGS84_EGM96'' ]\n');
fprintf('               gravDegree: [ integer>1 {20} ]\n');
fprintf('                gravOrder: [ integer>1 {20} ]\n');
fprintf('          useSolarGravity: [ ''true'' | {''false''} ]\n');
fprintf('          useLunarGravity: [ ''true'' | {''false''} ]\n');
fprintf('useSolarRadiationPressure: [ ''true'' | {''false''} ]\n');
fprintf('       useAtmosphericDrag: [ ''true'' | {''false''} ]\n');
fprintf('          atmosphereModel: [ {''HP''} | ''NRL'' ]\n');
fprintf('     nParameterForHPModel: [ 2 <= integer <= 6 {2} ]\n');
fprintf('              f107Average: [ scalar>=0 {150} ]\n');
fprintf('                f107Daily: [ scalar>=0 {150} ]\n');
fprintf('                       ap: [ scalar>=0 {4} ]\n');
end

function printMeasOptionsHelp()
fprintf('\n');
fprintf('MEASUREMENT OPTIONS\n');
fprintf('--------------------------\n');
fprintf('                    epoch: [ scalar>=0 {0} ]\n');
fprintf('                rangeType: [      ''1way'' | {''2way''}  ]\n');
fprintf('                 useRange: [    {''true''} | ''false''   ]\n');
fprintf('             useRangeRate: [    {''true''} | ''false''   ]\n');
fprintf('               useDoppler: [      ''true'' | {''false''} ]\n');
fprintf('                  useUnit: [      ''true'' | {''false''} ]\n');
fprintf('             useLightTime: [      ''true'' | {''false''} ]\n');
fprintf('         useGPSIonosphere: [      ''true'' | {''false''} ]\n');
fprintf('            useIonosphere: [      ''true'' | {''false''} ]\n');
fprintf('           useTroposphere: [      ''true'' | {''false''} ]\n');
fprintf('       useChargedParticle: [      ''true'' | {''false''} ]\n');
fprintf('        frequencyTransmit: [ scalar>0 {1.57542e9} ]\n');
fprintf('                   gsList: [ JAT groundStationList java object ]\n');
fprintf('                     gsID: [ cell array of strings ]\n');
fprintf('                   gsECEF: [ m x 3 matrix of ECEF locations]\n');
fprintf('    gsElevationConstraint: [ scalar>=0 {10} ]\n');
fprintf('                   rSigma: [ scalar>=0 {1e-3} ]\n');
fprintf('       EarthAtmMaskRadius: [ scalar>0 {6478.12} ]\n');
fprintf('                    tdrss: [ stucture for tdrssmeas.m ]\n');
fprintf('                    relay: [ stucture for lnrmeas.m ]\n');
fprintf('                     ddor: [ stucture for ddormeas.m ]\n');
fprintf('                     xnav: [ stucture for xnavmeas.m ]\n');
fprintf('                 Schedule: [ ground station tracking times for gsmeas.m ]\n');
fprintf('          PrecnNutnExpire: [ days {1} ]\n');
fprintf('             OpticalSigma: [ scalar>=0 {6e-10} ]\n');
fprintf('                   Camera: [ Camera Object ]\n');
fprintf('                 Attitude: [ Attitude Object ]\n');
fprintf('               linkbudget: [ structure for gpsmeas.m ]\n');
fprintf('\n');
end
