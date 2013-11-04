function     [gps_pos, gps_vel, vel_tot,dtsv] = rnxnEval(time,ephs, sv)
% Evaluates ICD-200 GPS Ephemeris at user specified time.  Uses ephemerides
% read in using read_rnxn function.  Several ephemerides may or may not be
% present depending on the RINEX file that is read with read_rnxn.  If
% multiple ephemerides are present for a PRN, this code will use the most
% recent relative to the TOE.  This function may be used to calculate GPS
% state for a single PRN or for the entire constellation.
%
% INPUTS:       time = total GPS Seconds [numberTimes x 1] double from start of GPS
%                               time (Jan 6,1980). See ODTBX convertTime() m-file.
%
%                    ephs = ICD-200 ephemerides, [numEphemerides x 34]
%                               Depending on RINEX File read in, there
%                               could be several ephemerides per PRN, as
%                               one elapses another becomes usable. 
%                    column breakdown (see read_rnxn.m file)
%                               Column==ID====Description====================
%                               1 short int prn;
%                               2 logical vflg;                                   0=No valid data
%                               3 long TofXmission;  Time of subframe 1 transmission, sec of week
%                               4 short int s1hlth;                        Subframe 1 health code
%                               5 short int codeL2;                               Code on L2 flag
%                               6 short int wkn;         GPS week at time of subframe 1 reception
%                               7 short int L2Pdata;                               L2 P data flag
%                               8 short int ura;                             Satellite's URA code
%                               9 short int iodc;                            Issue of data, clock
%                              10  double tgd;                               Group delay parameter
%                              11  short int tocwk;                  GPS week corresponding to toc
%                              12  double toc;          Reference time of clock data parameter set
%                              13  double af0;             Clock correction polynomial coefficient
%                              14  double af1;             Clock correction polynomial coefficient
%                              15  double af2;             Clock correction polynomial coefficient
%                              16  short int iode;                        Issue of data, ephemeris
%                              17  double crs;          Sine harmonic correction to orbital radius
%                              18  double deltan;            Mean motion delta from computed value
%                              19  double m0;                                  Mean anomaly at TOE
%                              20  double cuc;        Cosine harmonic correction to orbital radius
%                              21  double ecc;                                        Eccentricity
%                              22  double cus;          Sine harmonic corr to argument of latitude
%                              23  double sqrta;                     Square root of semimajor axis
%                              24  short int toewk;                  GPS week corresponding to toe
%                              25  double toe;                Reference time of ephemeris data set
%                              26  short int fti;                                     Fit interval
%                              27  double cic;                 Cosine harmonic corr to inclination
%                              28  double om0;                              Right ascension at TOE
%                              29  double cis;                   Sine harmonic corr to inclination
%                              30  double in0;                                  Inclination at TOE
%                              31  double crc;        Cosine harmonic correction to orbital radius
%                              32  double olc;                          Argument of perigee at TOE
%                              33  double omd;                             Rate of right ascension
%                              34  double idot;                                Rate of inclination
%
% OUTPUTS:   
%                    NOTE:  Data data is assigned by PRN in either the
%                    column number or 3-d dimension of a 3-d array
%                    (~,~,PRN)
%       
%                    gps_pos = WGS-84 ECEF Position [km]
%                                    [3 x numberTimes x 32]  OR [3 x numberTimes x 1]
%                    gps_vel  = WGS-84 ECEF Velocity [km/s]
%                                    [3 x numberTimes x 32]  OR [3 x numberTimes x 1]
%                    vel_tot   =  UNUSED, currently the same as gps_vel
%                    dtsv       = Accumulated Clock Bias for PRN

%% Pre-allocate Arrays and Define Constants
MU = 3.986005e14;              %  [m^3/s^2]
OMEGA_EARTH = 7.2921151467e-5;  % WGS-84 Earth's rotation rate [rad/s] from GPS ICD 1991

numTimes = length(time);
week = time(1)/(86400*7);
gpsSecs = time - floor(week)*86400*7;  % seconds into current GPS week

if nargin == 3
    gps_pos = NaN(3,numTimes,1);
    gps_vel  = NaN(3, numTimes,1);
    vel_tot = NaN(3, numTimes,1);
    dtsv = NaN(1,numTimes);
    prnList = sv;
elseif nargin ==2
    gps_pos = NaN(3,numTimes,32);
    gps_vel  = NaN(3, numTimes,32);
    vel_tot = NaN(3, numTimes,32);
    dtsv = NaN(32,numTimes);
    prnList = [1:32];
else
    warning('Not enough arguments provided to rnxnEval()');
end

%% Loop thru each PRN and Time to get GPS State(s)
for prn = prnList  % cycle thru all or just 1 PRN
    % Find the ephemerides that pertain to the PRN
    rows = find(ephs(:,1) == prn);
    
    if ~isempty(rows)
        ephs_k = ephs(rows,:);
    else
        continue;  % Ephemeris Unavailable for PRN
    end
    
    deltaToe = 0;  % Time elapsed from TOE in ephemeris
    
    for t_i = 1:numTimes % cycle thru all times
        % If just staring, find and use most recent ephemeris for PRN
        % If time is over 2 hours away from TOE, grab more recent ephemeris
        if t_i == 1 || deltaToe > 7200
            eph_k = getMostRecentEph(time(t_i), ephs_k);
        end
        
        % At each time step calculate time from TOE
        deltaToe = abs(eph_k(:,24)*7*86400 + eph_k(:,25) -  time(t_i));
        
        gpsWeek    = eph_k(6);         %GPS week
        toc     =   eph_k(12);   %Reference time of SV clock data
        af0     =   eph_k(13);   %Clock correction polynomial coefficient 1
        af1     =   eph_k(14);   %Clock correction polynomial coefficient 1
        af2     =   eph_k(15);   %Clock correction polynomial coefficient 1
        iode    =   eph_k(16);   %Issue of data ephemeris
        Crs     =   eph_k(17);   %Radius correction sin component
        delta_n =   eph_k(18);   %Mean motion
        M0     =   eph_k(19);   %Mean Anomaly at toe
        Cuc     =   eph_k(20);   %Argument of latitude correction cos component
        ecc     =   eph_k(21);   %Eccentricity
        Cus     =   eph_k(22);   %Argument of lattitude correction sin component
        sqrt_a  =   eph_k(23);   %Semi-Magor axis root
        week    =   eph_k(24);   %TOE week number
        toe     =   eph_k(25);   %TOE GPS seconds
        Cic     =   eph_k(27);   %Inclination correction cos component
        OMEGA_0 =   eph_k(28);   %Right ascention of the ascending node
        Cis     =   eph_k(29);   %Inclination correction sin component
        inc0   =   eph_k(30);   %Inclination at toe
        Crc     =   eph_k(31);   %Radius correction cos component
        omega   =   eph_k(32);   %Arguement of perigee
        OMEGA_dot = eph_k(33);   %Right ascention rate
        inc_dot =   eph_k(34);   %Inclination rate
        
        a         = (sqrt_a)^2;    % compute semimajor axis (meters)
        n_0       = sqrt(MU/a^3);  % computed mean motion (rad/s)
        n         = n_0 + delta_n; % corrected mean motion (rad/s)
        
        % Determine the elapsed time into the current gps week
        t_epoch = gpsSecs(t_i);
        
        % Elapsed time since the ephemeris has become valid
        t_k = t_epoch - toe;
        
        % Compute the Mean anomaly at time t
        M = M0 + (n*t_k);
        M = rem(M,2*pi);         
        
        %Solve Keppler's equation for the eccentric anomaly
        E = M;                      % Typically a good guess for the first iteration
        E_old = E+1;                % Make E_old large enough that the loop begins
        while sum(abs(E - E_old) >= 0.0000000005) > 0
            E_old = E;
            E = M + (ecc*sin(E_old));
        end
        
        % True Anomaly
        v = atan2(((sqrt(1 - ecc.^2))*sin(E)),(cos(E) - ecc));
        
        % Arg of Latitude
        phi = v + omega;
        
        % Generate the updates based upon the broadcast corrections
        delta_u = Cus*sin(2*phi) + Cuc*cos(2*phi);
        delta_r = Crs*sin(2*phi) + Crc*cos(2*phi);
        delta_i = Cis*sin(2*phi) + Cic*cos(2*phi);
        
        % Compute the updated orbital parameters
        u       = phi +delta_u;                     % Corrected argument of latitude
        r       = a*(1 - ecc*cos(E)) + delta_r;    % Corrected radius
        inc     = inc0 + delta_i + (inc_dot*t_k);    % Corrected inclination
        
        % Corrected longitude of ascending node
        OMEGA = OMEGA_0 + ((OMEGA_dot-OMEGA_EARTH)*t_k) - (OMEGA_EARTH.*toe);
        
        % Compute the satellite position
        x = r*((cos(OMEGA)*cos(u)) - (sin(OMEGA)*sin(u)*cos(inc)));
        y = r*((sin(OMEGA)*cos(u)) + (cos(OMEGA)*sin(u)*cos(inc)));
        z = r*sin(inc)*sin(u);
        pos = [x; y; z;];
        
                
        % Compute Sat Velocity
        p = a.*(1-ecc.^2);			   % semi-latus rectum
        h = sqrt(MU*p);				   % angular momentum
        vel(1) = (h/r)*((x*ecc*sin(v))/p - (cos(OMEGA)*sin(u) + sin(OMEGA)*cos(u).*cos(inc)));
        vel(2) = (h/r)*((y*ecc*sin(v))/p - (sin(OMEGA)*sin(u) - cos(OMEGA)*cos(u).*cos(inc)));
        vel(3) = (h/r)*((z*ecc*sin(v))/p + sin(inc)*cos(u));
        
        % Compute EARTH_RATE cross R
        o_cross_r = cross([0 0 OMEGA_EARTH],pos);
        
        % Remove the earth rate from the satellite velocity obtain a proper ECEF
        % velocity
        vecef = vel - o_cross_r;
        
        if length(prnList) == 1
            gps_pos(:,t_i,1) = pos/1000;
            gps_vel(:,t_i,1) = vecef'/1000;
            vel_tot(:,t_i,1) = vel'/1000; 
            dtsv(1,t_i) = af0 + af1*t_k + af2 * t_k * t_k;
        else
            gps_pos(:,t_i,prn) = pos/1000;
            gps_vel(:,t_i,prn) = vecef'/1000;
            vel_tot(:,t_i,prn) = vel'/1000;
            dtsv(prn,t_i) = af0 + af1*(t_epoch - toc) + af2 * (t_epoch - toc)^2;
        end
    end % end of times loop
end  % end of prnLIst Loop

end  % end of rnxnEval function

function eph_k = getMostRecentEph(gpsTotSecs,ephs_k)
% function grabs the appropriate ephemeris for the prn
% if the time is within a +/- 2 hours of the TOE it grabs that ephemeris
%
% INPUTS:
%               gpsTotSecs = total GPS seconds since 1980
%               ephs_k        = rinex ephemerides for satellite k (see above function)
%
% OUTPUTS: eph_k = [1x34] array of ephemeris info
%

%     gpsTime = (gpsWeek * 7 * 86400 + gpsSec);
tDiff = abs(ephs_k(:,24)*7*86400 + ephs_k(:,25) -  gpsTotSecs);
[~, rowToe] = min(tDiff);

% Store and return ephemeris for PRN, one with most recent TOE
eph_k = ephs_k(rowToe,:);

end
