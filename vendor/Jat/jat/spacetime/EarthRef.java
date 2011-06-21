/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
 * administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 *
 * This file is part of JAT. JAT is free software; you can
 * redistribute it and/or modify it under the terms of the
 * NASA Open Source Agreement, version 1.3 or later.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * NASA Open Source Agreement for more details.
 *
 * You should have received a copy of the NASA Open Source Agreement
 * along with this program; if not, write to the NASA Goddard
 * Space Flight Center at opensource@gsfc.nasa.gov.
 *
 */

package jat.spacetime;
import jat.cm.Constants;
import jat.eph.DE405;
import jat.eph.DE405_Body;
import jat.math.MathUtils;
import jat.matvec.data.Matrix;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.util.FileUtil;

/**
 * <P>
 * The EarthRef Class provides methods of converting between various types of Earth reference systems.
 * Reference: Satellite Orbits by Montenbruck and Gill. This is basically their C++ code converted to Java.
 * The interface to this class is an input Terrestrial time in MJD format, which is available from either
 * the CalDate or GPSTimeFormat classes.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class EarthRef extends BodyCenteredInertialRef implements BodyRef {

	private static final long serialVersionUID = 1L;
    
    public boolean use_moon = false;
    public boolean use_sun = false;
    //* TODO Watch this
    private boolean debug_polar = false;
    private boolean debug_geons = false;
    
    //private double GMST_REF, MJD_UT1_0;
    //private boolean gmst_initialized = false;
        
    /** Earth's rotation rate in rad/s.
     */
    //public final static double omega_e = 7.2921157746E-05;  // earth rotation rate
    public final static double omega_e = 7.292115E-05;  // IERS 1996 conventions
    //public final static double omega_e = 7.2921151467E-05;  // WGS-84
    
    //public double omega_e_dynamic = omega_e;
    
    /** Equatorial radius of earth in m from WGS-84
     */
    //public final static double R_Earth = 6378.137e3;      // Radius Earth [m]; WGS-84
    public final static double R_Earth = 6378.1363e3;      // Radius Earth [m]; STK JGM3
    
    /** Flattening factor of earth from WGS-84
     */
    //public final static double f_Earth = 1.0/298.257223563; // Flattening; WGS-84
    //public final static double f_Earth = 0.003353; // STK HPOP
    public final static double f_Earth = 0.00335281; // STK HPOP - 2
    
    /** Earth gravity constant in m^3/s^2 from JGM3
     */
    public final static double GM_Earth    = 398600.4415e+9;    // [m^3/s^2]; JGM3
    
    /** Two PI.
     */
    private static final double pi2 = 2.0*MathUtils.PI;
    
    /** Psi angle
     */
    private double dpsi;  // Celestial pole offsets, computed or gotten from IERS.
    /** Epsilon nutation angle.
     */
    private double deps;  // Celestial pole offsets, computed or gotten from IERS.
    
    private double Om; // Ascending node of Moon see Astro Almenac
    
    /** X Pole Coordinate in radians. */
    private double x_pole  = 0.0;          // Pole coordinate [rad]. From IERS Bulletin A.
    
    /** Y Pole Coordinate in radians. */
    private double y_pole  = 0.0;          // Pole coordinate [rad]. From IERS Bulletin A.
    
    /** Flag indicating poles are non-zero and should be calculated */
    private boolean use_poles = false;
    
    private Matrix T;
    
    private Matrix E;
    
    protected DE405 jpl_ephem;
    
    /** Geocentric position of the sun [m] in the J2000 (DE405) inertial frame 
     */
    private VectorN r_sun;
    
    /** Geocentric position of the moon [m] in the J2000 (DE405) inertial frame 
     */
    private VectorN r_moon; 

    /* The UT1 time, for when calculating the state of the 
     * Earth (e.g. translating to Greenwich Meridian System) */ 
    private double t_mjd_ut1;

    /* The TT time, for when calculating the state of the 
     * Earth (e.g. translating to True of Date) */ 
    private double t_mjd_tt;
    
    /** Construct an EarthRef object at a given time
     * @param mjd_utc Universal coordinated time in modified julian date
     */
    public EarthRef(double mjd_utc){
        this(new Time(mjd_utc));
    }
    
    /** Construct an EarthRef object at a given time
     * @param t0 JAT time
     */
    public EarthRef(Time t0){
      this(t0.mjd_ut1(), t0.mjd_tt(), false, false);
    }
    
    /** Construct an EarthRef object at a given time
     * @param MJD_UT1 UT1 time in MJD format.
     * @param MJD_TT TT time in MJD format.
     */
    public EarthRef(double MJD_UT1, double MJD_TT){
      this(MJD_UT1, MJD_TT, false, false);
    }
    
    /** Construct an EarthRef object at a given time
     * @param MJD_UT1 UT1 time in MJD format.
     * @param MJD_TT TT time in MJD format.
     * @param boolean whether Earth-to-Moon information should be precalculated
     * @param boolean whether Earth-to-Sun information should be precalculated
     */
    public EarthRef(double MJD_UT1, double MJD_TT, boolean use_moon, boolean use_sun){
        super(DE405_Body.EARTH);
      	t_mjd_ut1 = MJD_UT1;
      	t_mjd_tt = MJD_TT;
        this.T = trueOfDate(t_mjd_tt);
        this.E = tod2ecef(MJD_UT1, MJD_TT).times(T);

        this.use_sun = use_sun;
        this.use_moon = use_moon;
        if(this.use_moon || this.use_sun){
            String fs, dir_in;
            fs = FileUtil.file_separator();
            try{
                dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
            }catch(Exception e){
                dir_in = "C:/Code/Jat/jat/eph/DE405data/";
            }
            jpl_ephem = new DE405(dir_in);
 
            if(this.use_sun) compute_JPL_Sun_Vector(MJD_TT);
            if(this.use_moon) compute_JPL_Moon_Vector(MJD_TT);
        }
     }

    
    /** Set the IERS (Earth Rotation) Data
     * @param x x_pole coordinate (rad)
     * @param y y_pole coordinate (rad)
     */
    public void setIERS(double x, double y){
        this.x_pole = x*Constants.arcsec2rad;
        this.y_pole = y*Constants.arcsec2rad;
        use_poles = (x_pole != 0) || (y_pole != 0);
        //this.T = trueOfDate();
        //this.E = eci2ecef();
    }
    
    /**
     * Initialize the JPL DE405 ephemerides and initialize the Moon vector.
     *
     */
    public void initializeMoonEphem(double MJD_TT){
    	this.use_moon = true;
    	if(this.use_moon || this.use_sun){
            String fs, dir_in;
            fs = FileUtil.file_separator();
            try{
                dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
            }catch(Exception e){
                dir_in = "C:/Code/Jat/jat/eph/DE405data/";
            }
            jpl_ephem = new DE405(dir_in);
        }
    	//double JD_TDB = MJD_TT+2400000.5;
    	compute_JPL_Moon_Vector(MJD_TT);
    }
    
    /**
     * Initialize the JPL DE405 ephemerides and initialize the Sun vector.
     *
     */
    public void initializeSunEphem(double MJD_TT){
    	this.use_sun = true;
    	if(this.use_moon || this.use_sun){
            String fs, dir_in;
            fs = FileUtil.file_separator();
            try{
                dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
            }catch(Exception e){
                dir_in = "C:/Code/Jat/jat/eph/DE405data/";
            }
            jpl_ephem = new DE405(dir_in);
        }
    	//double JD_TDB = MJD_TT+2400000.5;
    	compute_JPL_Sun_Vector(MJD_TT);
    }

    /** Returns the (precalculated) ECI to ECEF (ICRF to ITRF) Transformation Matrix.
     * @return ECI to ECEF Transformation Matrix.
     */    
    public Matrix ECI2ECEF(){
        return this.E;
    }

    /** Returns the (precalculated) J2000 to True of Date (ICRF to TOD) Transformation Matrix.
     * @return J2000 to TOD Transformation Matrix.
     */        
    public Matrix TOD(){
        return this.T;
    }
    
    /** Computes the mean obliquity of the ecliptic. Uses Mjd_TT (Terrestrial Time).
     *@return  Mean obliquity of the ecliptic
     */
    public double MeanObliquity(double MJD_TT) {
        double T = (MJD_TT - TimeUtils.MJD_J2000)/36525.0;
        //double T = (this.MJD_TDB - MJD_J2000)/36525.0; //* used prior to IERS1996 convention
//      TODO
        return MathUtils.DEG2RAD *( 23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
        //* Following taken from GTDS documentation
        //return MathUtils.DEG2RAD *( 23.43929111-0.0130047*T-(0.1639e-6)*T*T+(0.5036e-6)*T*T*T);
        //* Debug - updated the epsilon-0 value according to the Supplemental Astro Almanac (old above)
        //return MathUtils.DEG2RAD *( 23.4392911111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
        //return Constants.arcsec2rad *( 84381.448-(46.8150+(0.00059-0.001813*T)*T)*T );
        //* Debug Vallado p 209
//        double T = Time.TTtoTDB(MJD_TT);
//        T = (T-Time.MJD_J2000)/36525.0;
//        return MathUtils.DEG2RAD * (23.439291 - 0.0130042*T-1.64e-7*T*T+5.04e-7*T*T*T);
//        //return Constants.arcsec2rad * (84381.948-(46.8150+(0.00059-0.001813*T)*T)*T );
    }
    
//*** <- The following is a backup of MeanObliquity as verified Summer 2005 -> ***
//    /** Computes the mean obliquity of the ecliptic. Uses Mjd_TT (Terrestrial Time).
//     *@return  Mean obliquity of the ecliptic
//     */
//    public double MeanObliquity(double MJD_TT) {
//        double T = (MJD_TT - TimeUtils.MJD_J2000)/36525.0;
//        //double T = (this.MJD_TDB - MJD_J2000)/36525.0; //* used prior to IERS1996 convention
////      TODO
//        return MathUtils.DEG2RAD *( 23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
//        //* Debug - updated the epsilon-0 value according to the Supplemental Astro Almanac (old above)
//        //return MathUtils.DEG2RAD *( 23.4392911111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
//        //return Constants.arcsec2rad *( 84381.448-(46.8150+(0.00059-0.001813*T)*T)*T );
//        //* Debug Vallado p 209
////        double T = Time.TTtoTDB(MJD_TT);
////        T = (T-Time.MJD_J2000)/36525.0;
////        return MathUtils.DEG2RAD * (23.439291 - 0.0130042*T-1.64e-7*T*T+5.04e-7*T*T*T);
////        //return Constants.arcsec2rad * (84381.948-(46.8150+(0.00059-0.001813*T)*T)*T );
//    }
    
    /** Transformation of equatorial to ecliptical coordinates. Uses Mjd_TT (Terrestrial Time).
     *   @return  Transformation matrix
     */
    public RotationMatrix EclMatrix(double MJD_TT) {
        RotationMatrix out = new RotationMatrix(1, MeanObliquity(MJD_TT));
        return out;
    }
    
    /** Precession transformation of equatorial coordinates. Uses Mjd_TT (Terrestrial Time).
     *  @return  Precession transformation matrix
     */
    public Matrix PrecMatrix(double MJD_TT) {
        double Mjd_1 = TimeUtils.MJD_J2000;  // Epoch given (Modified Julian Date TT)
        double Mjd_2 = MJD_TT;   // Epoch to precess to (Modified Julian Date TT)
        //double Mjd_2 = this.MJD_TDB;   //* used prior to IERS1996 convention
        
        // Constants
        final double Arcs = 1.0/MathUtils.ARCSEC2RAD;
        final double T  = (Mjd_1-TimeUtils.MJD_J2000)/36525.0;
        final double dT = (Mjd_2-Mjd_1)/36525.0;
        
        // Variables
        double zeta,z,theta;
        
        // Precession angles
        zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+
        ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
        z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs; 
        theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-
        ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;
        
        // Precession matrix
        RotationMatrix r1 = new RotationMatrix(3, -z);
        RotationMatrix r2 = new RotationMatrix(2, theta);
        RotationMatrix r3 = new RotationMatrix(3, -zeta);
        
        Matrix out = r1.times(r2.times(r3));
        return out;
    }
    
    /** Computes Nutation in longitude and obliquity using the IAU 1980 nutation theory.
     * Uses Mjd_TT (Terrestrial Time).
     */
    public void NutAngles(double MJD_TT) {
        double Mjd_TT = MJD_TT;
        //double Mjd_TT = this.MJD_TDB; //* used prior to IERS1996 convention
        
        // Constants
        
        final double T  = (Mjd_TT-TimeUtils.MJD_J2000)/36525.0;
        final double T2 = T*T;
        final double T3 = T2*T;
        final double rev = 360.0*3600.0;  // arcsec/revolution
        
        final int  N_coeff = 106;
        final long [][] C = {
            {  0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 },   //   1
            {  0, 0, 0, 0, 2,   20620,    2,   -8950,    5 },   //   2
            { -2, 0, 2, 0, 1,     460,    0,    -240,    0 },   //   3
            {  2, 0,-2, 0, 0,     110,    0,       0,    0 },   //   4
            { -2, 0, 2, 0, 2,     -30,    0,      10,    0 },   //   5
            {  1,-1, 0,-1, 0,     -30,    0,       0,    0 },   //   6
            {  0,-2, 2,-2, 1,     -20,    0,      10,    0 },   //   7
            {  2, 0,-2, 0, 1,      10,    0,       0,    0 },   //   8
            {  0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 },   //   9
            {  0, 1, 0, 0, 0,   14260,  -34,     540,   -1 },   //  10
            {  0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 },   //  11
            {  0,-1, 2,-2, 2,    2170,   -5,    -950,    3 },   //  12
            {  0, 0, 2,-2, 1,    1290,    1,    -700,    0 },   //  13
            {  2, 0, 0,-2, 0,     480,    0,      10,    0 },   //  14
            {  0, 0, 2,-2, 0,    -220,    0,       0,    0 },   //  15
            {  0, 2, 0, 0, 0,     170,   -1,       0,    0 },   //  16
            {  0, 1, 0, 0, 1,    -150,    0,      90,    0 },   //  17
            {  0, 2, 2,-2, 2,    -160,    1,      70,    0 },   //  18
            {  0,-1, 0, 0, 1,    -120,    0,      60,    0 },   //  19
            { -2, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  20
            {  0,-1, 2,-2, 1,     -50,    0,      30,    0 },   //  21
            {  2, 0, 0,-2, 1,      40,    0,     -20,    0 },   //  22
            {  0, 1, 2,-2, 1,      40,    0,     -20,    0 },   //  23
            {  1, 0, 0,-1, 0,     -40,    0,       0,    0 },   //  24
            {  2, 1, 0,-2, 0,      10,    0,       0,    0 },   //  25
            {  0, 0,-2, 2, 1,      10,    0,       0,    0 },   //  26
            {  0, 1,-2, 2, 0,     -10,    0,       0,    0 },   //  27
            {  0, 1, 0, 0, 2,      10,    0,       0,    0 },   //  28
            { -1, 0, 0, 1, 1,      10,    0,       0,    0 },   //  29
            {  0, 1, 2,-2, 0,     -10,    0,       0,    0 },   //  30
            {  0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 },   //  31
            {  1, 0, 0, 0, 0,    7120,    1,     -70,    0 },   //  32
            {  0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 },   //  33
            {  1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 },   //  34
            {  1, 0, 0,-2, 0,   -1580,    0,     -10,    0 },   //  35
            { -1, 0, 2, 0, 2,    1230,    0,    -530,    0 },   //  36
            {  0, 0, 0, 2, 0,     630,    0,     -20,    0 },   //  37
            {  1, 0, 0, 0, 1,     630,    1,    -330,    0 },   //  38
            { -1, 0, 0, 0, 1,    -580,   -1,     320,    0 },   //  39
            { -1, 0, 2, 2, 2,    -590,    0,     260,    0 },   //  40
            {  1, 0, 2, 0, 1,    -510,    0,     270,    0 },   //  41
            {  0, 0, 2, 2, 2,    -380,    0,     160,    0 },   //  42
            {  2, 0, 0, 0, 0,     290,    0,     -10,    0 },   //  43
            {  1, 0, 2,-2, 2,     290,    0,    -120,    0 },   //  44
            {  2, 0, 2, 0, 2,    -310,    0,     130,    0 },   //  45
            {  0, 0, 2, 0, 0,     260,    0,     -10,    0 },   //  46
            { -1, 0, 2, 0, 1,     210,    0,    -100,    0 },   //  47
            { -1, 0, 0, 2, 1,     160,    0,     -80,    0 },   //  48
            {  1, 0, 0,-2, 1,    -130,    0,      70,    0 },   //  49
            { -1, 0, 2, 2, 1,    -100,    0,      50,    0 },   //  50
            {  1, 1, 0,-2, 0,     -70,    0,       0,    0 },   //  51
            {  0, 1, 2, 0, 2,      70,    0,     -30,    0 },   //  52
            {  0,-1, 2, 0, 2,     -70,    0,      30,    0 },   //  53
            {  1, 0, 2, 2, 2,     -80,    0,      30,    0 },   //  54
            {  1, 0, 0, 2, 0,      60,    0,       0,    0 },   //  55
            {  2, 0, 2,-2, 2,      60,    0,     -30,    0 },   //  56
            {  0, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  57
            {  0, 0, 2, 2, 1,     -70,    0,      30,    0 },   //  58
            {  1, 0, 2,-2, 1,      60,    0,     -30,    0 },   //  59
            {  0, 0, 0,-2, 1,     -50,    0,      30,    0 },   //  60
            {  1,-1, 0, 0, 0,      50,    0,       0,    0 },   //  61
            {  2, 0, 2, 0, 1,     -50,    0,      30,    0 },   //  62
            {  0, 1, 0,-2, 0,     -40,    0,       0,    0 },   //  63
            {  1, 0,-2, 0, 0,      40,    0,       0,    0 },   //  64
            {  0, 0, 0, 1, 0,     -40,    0,       0,    0 },   //  65
            {  1, 1, 0, 0, 0,     -30,    0,       0,    0 },   //  66
            {  1, 0, 2, 0, 0,      30,    0,       0,    0 },   //  67
            {  1,-1, 2, 0, 2,     -30,    0,      10,    0 },   //  68
            { -1,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  69
            { -2, 0, 0, 0, 1,     -20,    0,      10,    0 },   //  70
            {  3, 0, 2, 0, 2,     -30,    0,      10,    0 },   //  71
            {  0,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  72
            {  1, 1, 2, 0, 2,      20,    0,     -10,    0 },   //  73
            { -1, 0, 2,-2, 1,     -20,    0,      10,    0 },   //  74
            {  2, 0, 0, 0, 1,      20,    0,     -10,    0 },   //  75
            {  1, 0, 0, 0, 2,     -20,    0,      10,    0 },   //  76
            {  3, 0, 0, 0, 0,      20,    0,       0,    0 },   //  77
            {  0, 0, 2, 1, 2,      20,    0,     -10,    0 },   //  78
            { -1, 0, 0, 0, 2,      10,    0,     -10,    0 },   //  79
            {  1, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  80
            { -2, 0, 2, 2, 2,      10,    0,     -10,    0 },   //  81
            { -1, 0, 2, 4, 2,     -20,    0,      10,    0 },   //  82
            {  2, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  83
            {  1, 1, 2,-2, 2,      10,    0,     -10,    0 },   //  84
            {  1, 0, 2, 2, 1,     -10,    0,      10,    0 },   //  85
            { -2, 0, 2, 4, 2,     -10,    0,      10,    0 },   //  86
            { -1, 0, 4, 0, 2,      10,    0,       0,    0 },   //  87
            {  1,-1, 0,-2, 0,      10,    0,       0,    0 },   //  88
            {  2, 0, 2,-2, 1,      10,    0,     -10,    0 },   //  89
            {  2, 0, 2, 2, 2,     -10,    0,       0,    0 },   //  90
            {  1, 0, 0, 2, 1,     -10,    0,       0,    0 },   //  91
            {  0, 0, 4,-2, 2,      10,    0,       0,    0 },   //  92
            {  3, 0, 2,-2, 2,      10,    0,       0,    0 },   //  93
            {  1, 0, 2,-2, 0,     -10,    0,       0,    0 },   //  94
            {  0, 1, 2, 0, 1,      10,    0,       0,    0 },   //  95
            { -1,-1, 0, 2, 1,      10,    0,       0,    0 },   //  96
            {  0, 0,-2, 0, 1,     -10,    0,       0,    0 },   //  97
            {  0, 0, 2,-1, 2,     -10,    0,       0,    0 },   //  98
            {  0, 1, 0, 2, 0,     -10,    0,       0,    0 },   //  99
            {  1, 0,-2,-2, 0,     -10,    0,       0,    0 },   // 100
            {  0,-1, 2, 0, 1,     -10,    0,       0,    0 },   // 101
            {  1, 1, 0,-2, 1,     -10,    0,       0,    0 },   // 102
            {  1, 0,-2, 2, 0,     -10,    0,       0,    0 },   // 103
            {  2, 0, 0, 2, 0,      10,    0,       0,    0 },   // 104
            {  0, 0, 2, 4, 2,     -10,    0,       0,    0 },   // 105
            {  0, 1, 0, 1, 0,      10,    0,       0,    0 }    // 106
        };
        
        // Variables
        double  l, lp, F, D;//, Om; //Om changed to class variable
        double  arg;
        
        // Mean arguments of luni-solar motion
        //
        //   l   mean anomaly of the Moon
        //   l'  mean anomaly of the Sun
        //   F   mean argument of latitude
        //   D   mean longitude elongation of the Moon from the Sun
        //   Om  mean longitude of the ascending node
        
        l  = MathUtils.Modulo(  485866.733 + (1325.0*rev +  715922.633)*T
        + 31.310*T2 + 0.064*T3, rev );
        lp = MathUtils.Modulo( 1287099.804 + (  99.0*rev + 1292581.224)*T
        -  0.577*T2 - 0.012*T3, rev );
        F  = MathUtils.Modulo(  335778.877 + (1342.0*rev +  295263.137)*T
        - 13.257*T2 + 0.011*T3, rev );
        D  = MathUtils.Modulo( 1072261.307 + (1236.0*rev + 1105601.328)*T
        -  6.891*T2 + 0.019*T3, rev );
//        Om = MathUtils.Modulo(  450160.280 - (   5.0*rev +  482890.539)*T
//         +  7.455*T2 + 0.008*T3, rev );
//      TODO
        //* Debug - changed to correspond to the Supplemental Astro Almanac (old above)
        //* Warning!!! This value below causes error ~8.816m for ISS ~1.264m Sun-Sync
//        Om = MathUtils.Modulo(  486160.28 - (   5.0*rev +  482890.539)*T
//          +  7.455*T2 + 0.008*T3, rev );
        //* Debug2 - IERS tech report 21 1996 conventions (below)
        Om = MathUtils.Modulo(  450160 - 6962890.2665*T
                +  7.4722*T2 + 0.007702*T3 - 0.00005939*T3*T, rev );
       
        // Nutation in longitude and obliquity [rad]
        final double Arcs = MathUtils.ARCSEC2RAD;
        
        deps = 0.0;
        dpsi = 0.0;
        for (int i=0; i<N_coeff; i++) {
            arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om ) * Arcs;
            dpsi += ( C[i][5]+C[i][6]*T ) * Math.sin(arg);
            deps += ( C[i][7]+C[i][8]*T ) * Math.cos(arg);
        }
        // TODO
        //* Delta Delta corrections
//        double ddpsi = -0.0534;
//        double ddeps = -0.00472;
//        this.dpsi = this.dpsi + ddpsi;
//        this.deps = this.deps + ddeps;
        
        this.dpsi = 1.0E-5 * dpsi*Arcs;
        this.deps = 1.0E-5 * deps*Arcs;
        this.Om = Om*Arcs;

        //* Vallado p 223
//        double mjd_TDB = Time.TTtoTDB(MJD_TT);
//        double TDB = (mjd_TDB-Time.MJD_J2000)/36525.0;
//        double TDB2 = TDB*TDB;
//        double TDB3 = TDB2*TDB;
//        double TDB4 = TDB3*TDB;
//        double r = 360;
//        		//* degrees
//        		l  = MathUtils.Modulo( 134.96340251 + (1325.0*r +  198.8675605)*TDB
//        		        + 0.0088553*TDB2 + (1.4343e-5)*TDB3-(6.797e-6)*TDB4, r );			
//                lp = MathUtils.Modulo( 357.52910918 + (  99.0*r + 359.0502911)*TDB
//                        -  0.0001537*TDB2 - 3.8e-8*TDB3 -3.19e-9*TDB4, r );
//                F  = MathUtils.Modulo( 93.27209062 + (1342.0*r +  82.0174577)*TDB
//                        - 0.0035420*TDB2 + 2.88e-7*TDB3+1.16e-9*TDB4, r );
//                D  = MathUtils.Modulo( 297.85019547 + (1236.0*r + 307.1114469)*TDB
//                        -  0.0017696*TDB2 + 1.831e-6*TDB3 - 8.80e-9*TDB4, r );
//                Om = MathUtils.Modulo( 125.04455501 - (   5.0*r +  134.1361851)*TDB
//                        +  0.0020756*TDB2 + 2.139e-6*TDB3 - 1.65e-8*TDB4, r );
//        final double deg2rad = MathUtils.PI/180.0;
//
//        // Nutation in longitude and obliquity [rad]
//        final double Arcs = MathUtils.ARCSEC2RAD;
//        
//        deps = 0.0;
//        dpsi = 0.0;
//        for (int i=0; i<N_coeff; i++) {
//            arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om ) * deg2rad;
//            dpsi += ( C[i][5]+C[i][6]*TDB ) * Math.sin(arg);
//            deps += ( C[i][7]+C[i][8]*TDB ) * Math.cos(arg);
//        }
//        // TODO
//        //* Delta Delta corrections
////        double ddpsi = -0.0534;
////        double ddeps = -0.00472;
////        this.dpsi = this.dpsi + ddpsi;
////        this.deps = this.deps + ddeps;
//        
//        this.dpsi = 1.0E-5 * dpsi*Arcs;
//        this.deps = 1.0E-5 * deps*Arcs;
//        //this.Om = Om*Arcs;

    }
    
    /** Transformation from mean to true equator and equinox. Uses Mjd_TT (Terrestrial Time).
     *   @return  Nutation matrix
     */
    public Matrix NutMatrix(double MJD_TT) {
//        double Mjd_TT = MJD_TT;
        //double Mjd_TT = this.MJD_TDB; //* used prior to IERS1996 convention
        
        // Mean obliquity of the ecliptic
        double eps = MeanObliquity(MJD_TT);
        
        // Nutation in longitude and obliquity
        NutAngles(MJD_TT);
        
        // Transformation from mean to true equator and equinox
        RotationMatrix r1 = new RotationMatrix(1, (-eps-deps));
        RotationMatrix r2 = new RotationMatrix(3, -dpsi);
        RotationMatrix r3 = new RotationMatrix(1, eps);
        
        Matrix out = r1.times(r2.times(r3));

        // Transformation from mean to true equator and equinox (Vallado)
//        RotationMatrix r1 = new RotationMatrix(1, (-eps-deps));
//        RotationMatrix r2 = new RotationMatrix(3, dpsi);
//        RotationMatrix r3 = new RotationMatrix(1, eps);
//        
//        Matrix out = r1.times(r2.times(r3));
//        out = out.transpose();
        return  out;
    }
    
    /** Transformation from mean to true equator and equinox (low precision).
     * Uses Mjd_TT (Terrestrial Time).
     * @return Nutation matrix
     */
    public Matrix NutMatrixSimple(double MJD_TT) {
        double Mjd_TT = MJD_TT;
        
        // Constants
        final double Arcs = 1.0/MathUtils.ARCSEC2RAD;
        final double T  = (Mjd_TT-TimeUtils.MJD_J2000)/36525.0;
        
        // Variables
        double  ls, D, F, N;
        double  eps;
        
        // Mean arguments of luni-solar motion
        ls = pi2*MathUtils.Frac(0.993133+  99.997306*T);   // mean anomaly Sun
        D  = pi2*MathUtils.Frac(0.827362+1236.853087*T);   // diff. longitude Moon-Sun
        F  = pi2*MathUtils.Frac(0.259089+1342.227826*T);   // mean argument of latitude
        N  = pi2*MathUtils.Frac(0.347346-   5.372447*T);   // longit. ascending node
        
        // Nutation angles
        dpsi = ( -17.200*Math.sin(N)   - 1.319*Math.sin(2*(F-D+N)) - 0.227*Math.sin(2*(F+N))
        + 0.206*Math.sin(2*N) + 0.143*Math.sin(ls) ) / Arcs;
        deps = ( + 9.203*Math.cos(N)   + 0.574*Math.cos(2*(F-D+N)) + 0.098*Math.cos(2*(F+N))
        - 0.090*Math.cos(2*N)                 ) / Arcs;
        
        // Mean obliquity of the ecliptic
        eps  = 0.4090928-2.2696E-4*T;
        
        RotationMatrix r1 = new RotationMatrix(1, (-eps-deps));
        RotationMatrix r2 = new RotationMatrix(3, -dpsi);
        RotationMatrix r3 = new RotationMatrix(1, eps);
        
        Matrix out = r1.times(r2.times(r3));
        return  out;
    }
    
    /** Computation of the equation of the equinoxes.
     * Uses Mjd_TT (Terrestrial Time).
     * Notes: The equation of the equinoxes dpsi*Math.cos(eps) is the right ascension of the
     * mean equinox referred to the true equator and equinox and is equal to the
     * difference between apparent and mean sidereal time.
     * @return Equation of the equinoxes
     */
    public double EqnEquinox(double MJD_TT) {
        // Nutation in longitude and obliquity
        //TODO NutAngles is called twice!!! ******************************* optimizeable
        NutAngles(MJD_TT);
        double arcs = Constants.arcsec2rad;
        // Equation of the equinoxes
//      TODO
//        double out =   dpsi * Math.cos( MeanObliquity(MJD_TT) );
//        return out;
        //* Montenbruck
//        return  dpsi * Math.cos( MeanObliquity() ) 
//        			+ (0.002649*Math.sin(Om)-0.000013*Math.cos(Om))*arcs;
        //* Astro Almenac 96 modified
        //* JAT Validated vs STK
        return  dpsi * Math.cos( MeanObliquity(MJD_TT) )
        			+ (0.002649*Math.sin(Om)+0.000063*Math.sin(2*Om))*arcs;
        //* Astro Almenac 96
//        return  dpsi * Math.cos( MeanObliquity(MJD_TT) )
//					+ (0.00264*Math.sin(Om)+0.000063*Math.sin(2*Om))*arcs;
    }
    
    /** Greenwich Mean Sidereal Time. Uses Mjd_UT1.
     *   @return  GMST in [rad]
     */
    public double GMST(double MJD_UT1) {
        //double Mjd_UT1 = MathUtils.round(MJD_UT1*1.0e6)/1.0e6;
    	double Mjd_UT1 = MJD_UT1;
        
        // Constants
        final double Secs = 86400.0;        // Seconds per day
        
        // Variables
        double Mjd_0,UT1,T_0,T,gmst;
        
        // Mean Sidereal Time
        Mjd_0 = Math.floor(Mjd_UT1);
        UT1   = Secs*(Mjd_UT1-Mjd_0);          // [s]
        T_0   = (Mjd_0  - TimeUtils.MJD_J2000)/36525.0;
        T     = (Mjd_UT1- TimeUtils.MJD_J2000)/36525.0;
        
        //* GEONS      
//        gmst  = (24110.54841 + 8640184.812866*T + 0.093104*T*T -6.2e-6*T*T*T); // [s]
//        //GMST_REF = gmst;
//        //MJD_UT1_0 = MJD_UT1;
//        //gmst_initialized = true;
//    	//- 0.3; // [s]
//        double tmp = pi2*(gmst/Secs)+omega_e*UT1;       // [rad], 0..2pi
//        tmp = MathUtils.Modulo(tmp,pi2);
//        if(tmp < 0) tmp = tmp + pi2;
//        return tmp;
    
        //* JAT Validated vs STK
        gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1+ (0.093104-6.2e-6*T)*T*T;  // [s]
        	//- 0.3; // [s]
        //???gmst  = 24110.54841 + 8640184.812866*T + (0.093104-6.2e-6*T)*T*T;  // [s]
        double tmp = pi2*MathUtils.Frac(gmst/Secs);       // [rad], 0..2pi
        //tmp = tmp + omega_e*UT1;
        tmp = MathUtils.Modulo(tmp,pi2);
        if(tmp < 0) tmp = tmp + pi2;
        return tmp;
      //TODO
        //* debug below IERS 1997
//        double r = 1.002737909350795 + (5.9006e-11 - 5.9e-15 * T_0)*T_0;
//        gmst  = 24110.54841 + 8640184.812866*T_0 + (0.093104-6.2e-6*T_0)*T_0*T_0 + r*UT1;
//        //double tmp =  gmst*Constants.omega_e;
//        double tmp;
//        tmp = pi2*MathUtils.Frac(gmst/Secs);
//        //tmp = MathUtils.Modulo(tmp,pi2);
//        return tmp;
        
        //* Vallado p 191
        //gmst  = 24110.54841 + 8640184.812866*T_0 + (0.093104-6.2e-6*T_0)*T_0*T_0; // [s]
        //gmst  = 100.4606184+36000.77005361*T_0+0.00038793*T_0*T_0-2.6e-8*T_0*T_0*T_0; // [deg]
//        double d1 = 67310.54841;
//        double d2 = (876600.0*3600.0+8640184.812866);
//        double d3 = 0.093104;
//        double d4 = -(6.2e-6);
//        gmst = d1 + d2*T + d3* T*T + d4 * T*T*T;
//        //gmst  = 67310.54841+(876600*3600+8640184.812866)*T+0.093104*T*T - (6.2e-6)*T*T*T; // [s]
//        double tmp = pi2*MathUtils.Frac(gmst/Secs);       // [rad], 0..2pi
//        tmp = MathUtils.Modulo(tmp,pi2);
//        //double tmp = MathUtils.Modulo(gmst*MathUtils.DEG2RAD,pi2);
//        //tmp = tmp + this.omega_e*UT1;//EarthRef.omega_e*UT1;
//        //if(tmp < 0) tmp = tmp + pi2;
//        return tmp;
        
    }
    
    /** Greenwich Apparent Sidereal Time.
     *   @return  GAST in [rad]
     */
    public double GAST(double MJD_UT1, double MJD_TT) {
        double out;
        //if(!gmst_initialized) 
        	//GMST(MJD_UT1);
        //if(this.debug_geons)
        	//out = MathUtils.Modulo( GMST_REF + omega_e*(MJD_UT1-MJD_UT1_0)*86400 + EqnEquinox(MJD_TT), pi2 );
        //else
        	//out = MathUtils.Modulo( GMST(MJD_UT1) + EqnEquinox(MJD_TT), pi2 );//- 2.55020752e-5;//6.55020752e-6;
        out = MathUtils.Modulo( GMST(MJD_UT1) + EqnEquinox(MJD_TT), pi2);// - 6.55020752e-6;
        return out;
    }
    
    /** Transformation from true equator and equinox to Earth equator and
     * Greenwich meridian system.
     * @return Greenwich Hour Angle matrix
     */
    public RotationMatrix GHAMatrix(double MJD_UT1, double MJD_TT) {
        RotationMatrix out = new RotationMatrix(3, GAST(MJD_UT1, MJD_TT));
        return  out;
    }
    
    /** Transformation from pseudo Earth-fixed to Earth-fixed coordinates
     * for a given date. Uses Mjd_UTC if poles are a function of UTC.
     * @return Pole matrix
     */
    public Matrix PoleMatrix() {
    	Matrix out;
    	if (use_poles) {
    		RotationMatrix r1 = new RotationMatrix(2, -this.x_pole);
    		RotationMatrix r2 = new RotationMatrix(1, -this.y_pole);
    		out = r1.times(r2);
    	}
    	else {
    		out = new Matrix(3);
    	}
        return  out;
    }
    
    public void computePole(Time t){
    	//return computePole(t.mjd_utc());
    }
    private void computePole(double mjd_utc){
    	
		double a1 = 0.14926633324398;
		double a2 = -0.34117340426258;
		double a3 = -1.8388673096747;
		double a4 =  0.10320829139742;
		double a5 = 2.2954920265308;
		double a6 = 0.030356241356650;
		double a7 = -0.95611083632580;
		double a8 = -1.7543492230047;
		double a9 =  1.3504353199051;
		double a10 = 2.1867543701143;
		double Tp = 52187.0;
    	
    	//* Case 1_1
//		double a1 = 0.13137346279621E+03;
//		double a2 = 0.31421954870224E+03;
//		double a3 = 0.19024293683469E+02;
//		double a4 =  -0.44553695678711E+03;
//		double a5 = -0.22882597982883E+02;
//		double a6 = 0.70715136051178E+02;
//		double a7 = 0.16834603309631E+03;
//		double a8 = 0.67812896966934E+01;
//		double a9 =  -0.23893760681152E+03;
//		double a10 = -0.80865914225578E+01;
//		double Tp = 52187.0;

    	//* Case 1_6
//		double a1 = -0.097873978689492;
//		double a2 = -0.024595601235973;
//		double a3 = -0.65797118746339;
//		double a4 =  0.095555969865679;
//		double a5 =  0.94278239745940;
//		double a6 =  0.13268750786818;
//		double a7 = -0.29741426364758;
//		double a8 = -0.38409384587563;
//		double a9 =  0.62823543875476;
//		double a10 = 0.53300376467789;
//		double Tp = 51013.0;
		double A = 2*Constants.pi/365.25*(mjd_utc-Tp);
		double C = 2*Constants.pi/435*(mjd_utc-Tp);
		double xp = a1 + a2 *Math.cos(A) + a3* Math.sin(A) + a4*Math.cos(C) + a5*Math.sin(C);
		double yp = a6 + a7 *Math.cos(A) + a8* Math.sin(A) + a9*Math.cos(C) + a10*Math.sin(C);
		xp = xp*Constants.arcsec2rad;
		yp = yp*Constants.arcsec2rad;
		this.x_pole = xp;
		this.y_pole = yp;
		use_poles = true;
	}
    
    /** J2000 to TOD Transformation
     * @return J2000 to ECEF transformation matrix
     */
    
    public Matrix trueOfDate(double MJD_TT){
        Matrix P = PrecMatrix(MJD_TT);
        Matrix N = NutMatrix(MJD_TT);
        Matrix out = N.times(P);
        return out;
    }
    
    /** Compute the ECI to ECEF tranformation matrix
     * @param the UT1 time in MJD format
     * @param the TT time in MJD format
     * @return ECI to ECEF transformation matrix
     */
    public Matrix eci2ecef(double MJD_UT1, double MJD_TT){
        return eci2ecef(MJD_UT1, MJD_TT, true);
    }
    
    /** Compute the ECI to ECEF tranformation matrix
     * @param the UT1 time in MJD format
     * @param the TT time in MJD format
     * @param this EarthRef object, when created, precalculated a True Of Date tranformation
     * matrix for the specified time.  If this is false, it will use the precalculated matrix even
     * if it was calculated for a different time
     * @return ECI to ECEF transformation matrix
     */
    public Matrix eci2ecef(double MJD_UT1, double MJD_TT, boolean recalcTOD){
    	Matrix Ecef;
    	if (MJD_TT == t_mjd_tt) {
    		// Don't bother calculating.  We've already precalculated.
    		Ecef = E;
    	}
    	else {
    		Matrix Tod = (recalcTOD ? trueOfDate(MJD_TT) : T);
    		Matrix A = tod2ecef(MJD_UT1, MJD_TT);
    		Ecef = A.times(Tod);
    	}
        return Ecef;
    }
    
    /** ECI to ECEF Transformation
     * @return ECI to ECEF transformation matrix
     */
    public Matrix eci2ecef(Time t){
        return eci2ecef(t.mjd_ut1(), t.mjd_tt(), false);
    }
    
    public Matrix tod2ecef(double mjd_ut1, double mjd_tt) {
        Matrix G = GHAMatrix(mjd_ut1, mjd_tt);
        
        if(debug_geons) {
        	computePole(mjd_ut1);
        }
        
        Matrix out = G;
        if(use_poles && !debug_polar){
        	out = PoleMatrix().times(G);
        }
    	return out;
    }
    
    /**
     * GEONS transformation between inertial and earth fixed
     * @param recf earth fixed position
     * @param vecf earth fixed velocity
     * @param t time
     * @return vector containing position and then velocity in inertial frame
     */
    public VectorN ecf2eci(VectorN recf, VectorN vecf, Time t){
//  	Compute derivative of GHA Matrix (S) and its transpose
    	double omega = Constants.omega_e;//Constants.WE_WGS84;
    	Matrix C = trueOfDate(t.mjd_tt());
        Matrix Rg = GHAMatrix(t.mjd_ut1(), t.mjd_tt());
        Matrix B;
        if(debug_geons) {
        	computePole(t);
        }
        if(!use_poles || debug_polar){
        	B = new Matrix(3);
        }else{
        	B = PoleMatrix();
        }
        //Matrix Pole = new Matrix(3);
        Matrix A = B.times(Rg);
        Matrix E = A.times(C);
    	VectorN omegaE = new VectorN(0,0,omega);

    	VectorN reci = C.transpose().times(Rg.transpose().times(B.transpose().times(recf)));

    	VectorN rpef = B.transpose().times(recf);
    	VectorN vpef = B.transpose().times(vecf);
    	Matrix Rgdot = new Matrix(3,3);
    	double ag = GAST(t.mjd_ut1(),t.mjd_tt());
    	Rgdot.A[0][0] = -omega*Math.sin(ag);
    	Rgdot.A[0][1] = omega*Math.cos(ag);
    	Rgdot.A[1][0] = -omega*Math.cos(ag);
    	Rgdot.A[1][1] = -omega*Math.sin(ag);
    	VectorN v_pole = B.transpose().times(vecf);
    	VectorN omegar = omegaE.crossProduct(recf);
    	VectorN sum = v_pole.plus(omegar);
    	VectorN sidereal = Rg.transpose().times(sum);
    	VectorN PN = C.transpose().times(sidereal);
    	VectorN veci = PN;
    	//VectorN veci = (C.transpose().times(Rgdot.transpose().times(recf))).plus(C.transpose().times(Rg.transpose().times(B.transpose().times(vecf))));
    	//VectorN veci = C.transpose().times(Rg.transpose().times((B.transpose().times(vecf)).plus(omegaE.crossProduct(recf))));
    	//VectorN veci = E.transpose().times(vecf).plus(omegaE.crossProduct(reci));
    	VectorN out = new VectorN(reci,veci);
    	return out;
    }
    /**
     * GEONS transformation between inertial and earth fixed
     * @param recf earth fixed position
     * @param vecf earth fixed velocity
     * @param t time
     * @return vector containing position and then velocity in inertial frame
     */
    public VectorN eci2ecf(VectorN reci, VectorN veci, Time t){
//  	Compute derivative of GHA Matrix (S) and its transpose
    	double omega = Constants.omega_e;//Constants.WE_WGS84;
    	Matrix C = trueOfDate(t.mjd_tt());
        Matrix Rg = GHAMatrix(t.mjd_ut1(), t.mjd_tt());
        Matrix B;
        if(debug_geons)
        	computePole(t);
        if(debug_polar){
        	B = new Matrix(3);
        }else{
        	B = PoleMatrix();
        }
        //Matrix Pole = new Matrix(3);
        Matrix A = B.times(Rg);
        Matrix E = A.times(C);
    	VectorN omegaE = new VectorN(0,0,omega);

    	VectorN recf = E.times(reci);

    	VectorN rpef = Rg.times(C.times(reci));//st*nut*prec*reci;
    	VectorN vecf = B.times((Rg.times(C.times(veci)).minus(omegaE.crossProduct(rpef))));
    	VectorN out = new VectorN(recf,vecf);
    	return out;
    }
    
    /** Updates the Earth model.
     * 
     * @param MJD_UT1.  Universal Time in modified julian date
     * @param MJD_TT.  Terrestrial Dynamical Time in modified julian date
     */
    public void update(double MJD_UT1, double MJD_TT){
        this.T = trueOfDate(MJD_TT);
        this.E = eci2ecef(MJD_UT1, MJD_TT);
        if(this.use_sun) compute_JPL_Sun_Vector(MJD_TT);
        if(this.use_moon) compute_JPL_Moon_Vector(MJD_TT);
        //this.omega_e_dynamic = get_omega_e(MJD_UT1);
    }
    
    /**
     * @see update(double MJD_UT1, double MJD_TT)
     * @param time
     */
    public void update(Time time){
    	update(time.mjd_ut1(),time.mjd_tt());
    }
    
    /** Computes the Sun's geocentric position using a low precision analytical series.
     * @return Solar position vector [m] with respect to the mean equator and equinox of J2000 (EME2000, ICRF)
     */
    public static VectorN sunVector(double MJD_TT) {
        double Mjd_TT = MJD_TT;
        // Constants
        
        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
        double T   = (Mjd_TT-TimeUtils.MJD_J2000)/36525.0;  // Julian cent. since J2000
        
        
        // Mean anomaly, ecliptic longitude and radius
        
        double M = pi2 * MathUtils.Frac( 0.9931267 + 99.9973583*T);                    // [rad]
        double L = pi2 * MathUtils.Frac( 0.7859444 + M/pi2 +(6892.0*Math.sin(M)+72.0*Math.sin(2.0*M))/1296.0e3); // [rad]
        double r = 149.619e9 - 2.499e9*Math.cos(M) - 0.021e9*Math.cos(2*M);             // [m]
        
        // Equatorial position vector
        
        RotationMatrix R = new RotationMatrix(1, -eps);
        
        VectorN temp = new VectorN(r*Math.cos(L),r*Math.sin(L),0.0);
        VectorN r_Sun = R.times(temp);
        
        return r_Sun;
        
    }
    
    /**
     * Get the geocentric position of the sun [m] in the J2000 (DE405) inertial frame
     * @param t - Time object (not used!)
     * @return Vector 3 [m]
     */
    public VectorN get_JPL_Sun_Vector(Time t){
        // We assume (hope) the time passed in is the same
        // as the current time.
    	
    	// TODO: not using the given Time is an interface violation - this is dangerous
    	
        if(this.use_sun)
            return r_sun;
        else
            return new VectorN(3);
    }

    /**
     * Get the geocentric position of the moon [m] in the J2000 (DE405) inertial frame 
     * at the current time
     * @return Vector 3 [m]
     */
    public VectorN get_JPL_Moon_Vector(){
        if(this.use_moon)
            return r_moon;
        else
            return new VectorN(3);
    }
    /**
     * Return the dynamic rotation rate of the Earth at the given Time.
     * @param t Time.
     * @return Omega [rad/s]
     */
    public double get_omega_e(Time t){
//      Variables
        double Mjd_0;
        
        // Mean Sidereal Time
        Mjd_0 = Math.floor(t.mjd_ut1());
        double Tu   = (Mjd_0  - TimeUtils.MJD_J2000)/36525.0;
        return 7292115.8553e-11 + 4.3e-15 * Tu;
    }
    
    /**
     * Return the dynamic rotation rate of the Earth at the given time.
     * @param mjd_ut1 time.
     * @return Omega [rad/s]
     */
    public double get_omega_e(double mjd_ut1){
//      Variables
        double Mjd_0;
        
        // Mean Sidereal Time
        Mjd_0 = Math.floor(mjd_ut1);
        double Tu   = (Mjd_0  - TimeUtils.MJD_J2000)/36525.0;
        return 7292115.8553e-11 + 4.3e-15 * Tu;
    }
    
    /**
     * Compute the JPL Sun Vector. (Converts to TDB).
     * @param MJD_TT Time.
     */
    private void compute_JPL_Sun_Vector(double MJD_TT){
        if(this.use_sun){
            r_sun = new VectorN(jpl_ephem.get_planet_pos(DE405_Body.GEOCENTRIC_SUN, MJD_TT));
        }
    }
    
    /**
     * Compute the JPL Moon Vector. (Converts to TDB).
     * @param MJD_TT Time.
     */
    private void compute_JPL_Moon_Vector(double MJD_TT){
        if(this.use_moon){
            r_moon = new VectorN(jpl_ephem.get_planet_pos(DE405_Body.GEOCENTRIC_MOON, MJD_TT));
        }
    }
    /**
     * Set the flag whether to calculate the Sun's position.
     * @param b
     */
    public void set_use_sun(boolean b){
        this.use_sun = b;
        if(b){
            if(jpl_ephem == null){
                String fs = FileUtil.file_separator();
                String dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
                jpl_ephem = new DE405(dir_in);
            }
        }
    }
    /**
     * Set the flag whether to calculate the Moon's position.
     * @param b
     */
    public void set_use_moon(boolean b){
        this.use_moon = b;
        if(b){
            if(jpl_ephem == null){
                String fs = FileUtil.file_separator();
                String dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
                jpl_ephem = new DE405(dir_in);
            }
        }
    }
    
    /** Computes the Moon's geocentric position using a low precision analytical series.
     * @return Lunar position vector [m] with respect to the mean equator and equinox of J2000 (EME2000, ICRF).
     */
    public VectorN moonVector(double MJD_TT) {
        double Mjd_TT = MJD_TT;
        
        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
        double T   = (Mjd_TT-TimeUtils.MJD_J2000)/36525.0;  // Julian cent. since J2000
        
        
        // Mean elements of lunar orbit
        
        double L_0 =     MathUtils.Frac( 0.606433 + 1336.851344*T );     // Mean longitude [rev]
        // w.r.t. J2000 equinox
        double l   = pi2*MathUtils.Frac( 0.374897 + 1325.552410*T );     // Moon's mean anomaly [rad]
        double lp  = pi2*MathUtils.Frac( 0.993133 +   99.997361*T );     // Sun's mean anomaly [rad]
        double D   = pi2*MathUtils.Frac( 0.827361 + 1236.853086*T );     // Diff. long. Moon-Sun [rad]
        double F   = pi2*MathUtils.Frac( 0.259086 + 1342.227825*T );     // Argument of latitude
        
        
        // Ecliptic longitude (w.r.t. equinox of J2000)
        
        double dL = +22640*Math.sin(l) - 4586*Math.sin(l-2*D) + 2370*Math.sin(2*D) +  769*Math.sin(2*l)
        -668*Math.sin(lp) - 412*Math.sin(2*F) - 212*Math.sin(2*l-2*D) - 206*Math.sin(l+lp-2*D)
        +192*Math.sin(l+2*D) - 165*Math.sin(lp-2*D) - 125*Math.sin(D) - 110*Math.sin(l+lp)
        +148*Math.sin(l-lp) - 55*Math.sin(2*F-2*D);
        
        double L = pi2 * MathUtils.Frac( L_0 + dL/1296.0e3 );  // [rad]
        
        // Ecliptic latitude
        
        double S  = F + (dL+412*Math.sin(2*F)+541*Math.sin(lp)) * MathUtils.ARCSEC2RAD;
        double h  = F-2*D;
        double N  = -526*Math.sin(h) + 44*Math.sin(l+h) - 31*Math.sin(-l+h) - 23*Math.sin(lp+h)
        +11*Math.sin(-lp+h) - 25*Math.sin(-2*l+F) + 21*Math.sin(-l+F);
        
        double B = ( 18520.0*Math.sin(S) + N ) * MathUtils.ARCSEC2RAD;   // [rad]
        
        double cosB = Math.cos(B);
        
        // Distance [m]
        
        double R = 385000e3 - 20905e3*Math.cos(l) - 3699e3*Math.cos(2*D-l) - 2956e3*Math.cos(2*D)
        -570e3*Math.cos(2*l) + 246e3*Math.cos(2*l-2*D) - 205e3*Math.cos(lp-2*D)
        -171e3*Math.cos(l+2*D) - 152e3*Math.cos(l+lp-2*D);
        
        // Equatorial coordinates
        
        VectorN temp = new VectorN( R*Math.cos(L)*cosB, R*Math.sin(L)*cosB, R*Math.sin(B) );
        RotationMatrix Rmat = new RotationMatrix(1, -eps);
        VectorN r_Moon = Rmat.times(temp);
        
        return r_Moon;
        
    }

    /**
     * Returns the transformation between inertial and body coordinates.
     * @param mjd_body Earth referenced time (Universal Time UT1)
     * @param mjd_inertial Dynamical time (Terrestrial Time)
     */
    public Matrix inertial_to_body(Time t) {
        return this.E;//eci2ecef(t);
    }

    /**
     * Returns the (precalculated) eci2ecef transformation.
     * @see jat.spacetime.BodyRef#body_to_inertial(Time)
     */
    public Matrix body_to_inertial(Time t) {
        //Matrix E = eci2ecef(t);
        return this.E.transpose();
    }
    /**
     * Returnts the dynamic Earth rotation rate. 
     * @see jat.spacetime.BodyRef#get_spin_rate()
     */
    public double get_spin_rate(Time t) {
        return get_omega_e(t);
    }
    /**
     * Returns the mean Earth Radius. 
     * @see jat.spacetime.BodyRef#get_mean_radius()
     */
    public double get_mean_radius() {
        return EarthRef.R_Earth;
    }
    /**
     * Returns the Earth's gravitational constant. 
     * @see jat.spacetime.BodyRef#get_grav_const()
     */
    public double get_grav_const() {
        return EarthRef.GM_Earth;
    }

    /**
     * Returns the trueOfDate trasformation.
     * @see jat.spacetime.BodyRef#trueOfDate(jat.spacetime.Time)
     */
    public Matrix trueOfDate(Time t) {
        return trueOfDate(t.mjd_tt());
    }
  
    /** OD Toolbox interface to Precession transformation of equatorial coordinates.
     *  @return  Precession transformation matrix
     */
    public Matrix PrecMatrix() {
    	return PrecMatrix(t_mjd_tt);    
    }
    
    /** OD Toolbox interface to Transformation from mean to true equator and equinox.
     *  @return  Nutation transformation matrix
     */
    public Matrix NutMatrix() {
    	return NutMatrix(t_mjd_tt);    
    }
    
    /** OD Toolbox interface to Transformation from true equator and equinox to 
     * Earth equator and Greenwich meridian system.
     *  @return  Greenwich Hour Angle transformation matrix
     */
    public Matrix GHAMatrix() {
    	return GHAMatrix(t_mjd_ut1,t_mjd_tt);    
    }
    
    
    /** OD Toolbox interface to Transformation of equatorial to ecliptical coordinates
     *  @return  transformation matrix
     */
    public Matrix EclMatrix() {
    	return EclMatrix(t_mjd_tt);    
    }

    /**
     * Test method.  See Vallado example 3-14.
     * @param args
     */
    public static void main(String[] args) {
        //* Follows Example 3-14 in Vallado
//        CalDate date = new CalDate(1991,4,6,7,51,28.386009);
//        double mjd_utc = date.mjd();
//        Time t = new Time(mjd_utc);
//        FitIERS iers = new FitIERS();
//        iers.process();
//        //double[] param = iers.search(mjd_utc);
//        double[] param = {-0.21959,0.30266,0.402521};
//        t.set_UT1_UTC(param[2]);
//        t.update(0);
//        EarthRef eRef = new EarthRef(t);
//        eRef.setIERS(param[0],param[1]);
//        Matrix E = eRef.eci2ecef(t);
//        double[] xd = {5102.5096, 6123.01152, 6378.1363, -4.7432196, 0.7905366, 5.53375619};
//        VectorN x = new VectorN(xd);
//        Matrix M = eRef.PrecMatrix(t.mjd_tt());
//        VectorN rmod = M.times(x.get(0,3));
//        Matrix N = eRef.NutMatrix(t.mjd_tt());
//        VectorN rtod = N.times(rmod);
//        Matrix SD = eRef.GHAMatrix(t.mjd_ut1(),t.mjd_tt());
//        VectorN rpef = SD.times(rtod);
//        Matrix P = eRef.PoleMatrix();
//        VectorN recef = P.times(rpef);
//        //double[] xvallado = {5102.509433,6123.011473,6378.136478,-4.74321966,0.79053639,5.53375617};
//        double[] xvallado = {-1120.598506,7894.483204,6374.079611,-3.18701800,-2.90527125,5.53765280};
//        VectorN vallado = new VectorN(xvallado);
//        System.out.println("error position: "+vallado.get(0,3).minus(recef));
//        
//        VectorN test = eRef.eci2ecf(x.get(0,3),x.get(3,3),t);
//        VectorN error = test.minus(vallado);
//        System.out.println("error: "+error);
//        
//        RotationMatrix E1 = new RotationMatrix(eRef.ECI2ECEF());
//        VectorN test2 = E1.times(x.get(0,3));
//        RotationMatrix E2 = new RotationMatrix(eRef.ECI2ECEF().transpose());
//        VectorN test3 = E2.times(test2.get(0,3));
//        VectorN error2 = test2.minus(vallado.get(0,3));
//        VectorN error3 = test3.minus(x.get(0,3));
//        System.out.println("error2: "+error2);
//        System.out.println("error3: "+error3);
  
    	
    	EarthRef eRef = new EarthRef(new Time(TimeUtils.MJD_J2000));
    	Time t = new Time(TimeUtils.TTtoUTC(TimeUtils.MJD_J2000));    	
    	Matrix E = eRef.eci2ecef(t.mjd_ut1(),t.mjd_tt());
    	System.out.println(""+E.toString());
//        int year = 1992;
//        int month = 8;
//        int day = 20;
//        int hour = 12;
//        int min = 14;  // UT1
//        
//        Time ut1 = new Time(year,month,day,hour,min,0);
//        EarthRef eRef = new EarthRef(ut1);
//        double gmst = eRef.GMST(ut1.mjd_ut1());
//        double gmst_v = MathUtils.Modulo(-232984181.090015915595*Constants.arcsec2rad,2*Constants.pi);
//        System.out.println("gmst: "+gmst);
//        System.out.println("error_gmst: "+(gmst-gmst_v));
//        
//        System.out.println("done");
    }
    
}
