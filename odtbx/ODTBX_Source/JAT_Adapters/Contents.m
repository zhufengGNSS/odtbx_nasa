% JAT_ADAPTERS
%
% Functions
%
%   convertTime             - Converts the input time from one timescale to another
%   createGroundStationList - Returns a JAT NASA_GroundStationList object
%   createJATWorld          - Creates a JATModel java object for propagating earth orbits.
%   ecef2LLA                - Returns the latitude, longitude and altitude of an ECEF position.
%   ephemDE405              - Computes the position and velocity of the input planet at the given input time.
%   getGroundStationInfo    - Returns the value of the specified ground station parameter.
%   getIERSTimes            - Returns the applicable times for the JAT IERS Parameters.
%   getIndex                - Returns the index within strArray of the entry that starts with str.
%   getJatRK8Options        - Get adaptor options for jatRK8.
%   jatChargedParticleModel - Returns the light time delay error due to charged particles
%   JATConstant             - Returns the value from JAT of the input constant from the specified source.
%   jatDCM                  - Returns the requested Direction Cosine Matrix from JAT.
%   jatForces               - Returns derivatives of JAT Force models (state in meters).
%   jatForces_km            - Returns derivatives of JAT Force models (state in km).
%   jatGPSIonoModel         - Returns the ionosphere range error for GPS code
%   jatHCW                  - Returns time derivatives of Hills-Clohessy-Wiltshire equations.
%   jatIonoDelayModel       - Ionospheric delay distance between satellites and ground stations
%   jatRK8                  - The wrapper for the JAT Runge-Kutta 8th order integrator.
%   jatStaAzEl              - Returns Azimuth and Elevation of satellites relative to ground stations
%   jatTropoModel           - Returns tropospheric refraction errors for range and elevation measurements.
%   jatWorldPropagatorRK4   - JAT 4th order Runge-Kutta integrator (Earth orbits)
%   jatWorldPropagatorRK8   - JAT 8th order Runge-Kutta integrator (Earth orbits)
%   jatWorldPropagatorRKF78 - JAT Runge-Kutta-Fehlberg variable step integrator (Earth orbits)
%   jD2MatlabTime           - Returns the datenum of the input Julian Date
%   jD2MJD                  - Converts Julian Date to Modified Julian Date format.
%   LLA2ecef                - Returns the ECEF position from latitude, longitude and altitude
%   matlabTime2JD           - Returns the Julian Date of the input datenum
%   matlabTime2MJD          - Returns the Modified Julian Date of the input time
%   mJD2JD                  - Converts Modified Julian Date to Julian Date.
%   mJD2MatlabTime          - Returns the datenum of the input Modified Julian Date
%   setJatRK8Options        - Set adaptor options for jatRK8.
%   jatForcesCon            - Returns derivatives of JAT Force models (state in meters).
