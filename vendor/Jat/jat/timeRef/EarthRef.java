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

package jat.timeRef;
import jat.matvec.data.*;
import jat.math.*;
import jat.spacetime.BodyCenteredInertialRef;
import jat.spacetime.ReferenceFrame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.*;
import jat.util.FileUtil;
import jat.cm.*;
import jat.eph.*;

/**
 * <P>
 * The EarthRef Class provides methods of converting between various types of Earth reference systems.
 * Reference: Satellite Orbits by Montenbruck and Gill. This is basically their C++ code converted to Java.
 * The interface to this class is an input Terrestrial time in MJD format, which is available from either
 * the CalDate or GPSTimeFormat classes.
 * 
 * For general use of Earth based coordinate systems, see jat.spacetime.EarthRef.java
 *
 * @deprecated
 * @see jat.spacetime.EarthRef
 * @author Richard C. Page III 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 * @see jat.spacetime.EarthRef
 */
public class EarthRef extends BodyCenteredInertialRef implements jat.spacetime.BodyRef {

	private static final long serialVersionUID = 1L;
    
    protected double sim_time = 0;
    public int n = 0;
    public int increment = 0;//1800;
    public boolean use_moon = true;
    public boolean use_sun = true;
    //public boolean use_iers = true;
        
    /** Modified Julian Date of the J2000 Epoch.
     */
    public static final double MJD_J2000 = 51544.5;
    
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
    
    /** Modified Julian Date.
     */
    private double MJD_TT;
    
    private double MJD_TDB;
    
    private double MJD_UTC;
    
    private double MJD_UT1;
    
    private double MJD_UTC_START;
    
    /** Difference between UT1 and UTC.
     */
    private double UT1_UTC=0;          // UT1-UTC time difference [s]. From IERS Bulletin A.
    
    /** X Pole Coordinate in radians.
     */
    private double x_pole  = 0.0;          // Pole coordinate [rad]. From IERS Bulletin A.
    
    /** Y Pole Coordinate in radians.
     */
    private double y_pole  = 0.0;          // Pole coordinate [rad]. From IERS Bulletin A.
    
    private Matrix T;
    
    private Matrix E;
    
    /** Geocentric position of the sun [m] in the J2000 (DE405) inertial frame 
     */
    private VectorN r_sun;
    
    /** Geocentric position of the moon [m] in the J2000 (DE405) inertial frame 
     */
    private VectorN r_moon; 
    
    private DE405 jpl_ephem;
    
    /** Construct an EarthRef object using UTC Time in Modified Julian Date format.
     * @param mjd_UTC UTC time in MJD format.
     */
    public EarthRef( double mjd_UTC ){
        super(DE405_Body.EARTH);
        //accept the input
//        CalDate utc = new CalDate(mjd_UTC);
        this.MJD_UTC = mjd_UTC;
        this.MJD_UTC_START = mjd_UTC;
        this.MJD_TT = CalDate.UTC2TT(mjd_UTC);
        this.MJD_TDB = TimeUtils.TTtoTDB(this.MJD_TT);
        this.MJD_UT1 = this.MJD_UTC + this.UT1_UTC/86400.0;
        this.T = trueOfDate();
        this.E = eci2ecef();
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
        if(this.use_sun) compute_JPL_Sun_Vector();
        if(this.use_moon) compute_JPL_Moon_Vector();
    }
    
    /** Construct an EarthRef object using UTC Time in CalDate format.
     * @param date CalDate object containing UTC time.
     */
    public EarthRef(CalDate date){
        super(DE405_Body.EARTH);
        this.MJD_UTC = date.mjd();
        this.MJD_UTC_START = date.mjd();
        this.MJD_TT = CalDate.UTC2TT(this.MJD_UTC);
        this.MJD_UT1 = this.MJD_UTC + this.UT1_UTC/86400.0;
        this.T = trueOfDate();
        this.E = eci2ecef();
        if(this.use_moon || this.use_sun){
            String fs = FileUtil.file_separator();
            String dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
            jpl_ephem = new DE405(dir_in);
        }
        if(this.use_sun) compute_JPL_Sun_Vector();
        if(this.use_moon) compute_JPL_Moon_Vector();
    }

    // For use with matlab
    public EarthRef( double mjd_UTC , boolean usingmatlab){
        super(DE405_Body.EARTH);
        //accept the input
//        CalDate utc = new CalDate(mjd_UTC);
        this.MJD_UTC = mjd_UTC;
        this.MJD_UTC_START = mjd_UTC;
        this.MJD_TT = CalDate.UTC2TT(mjd_UTC);
        this.MJD_TDB = TimeUtils.TTtoTDB(this.MJD_TT);
        this.MJD_UT1 = this.MJD_UTC + this.UT1_UTC/86400.0;
        this.T = trueOfDate();
        this.E = eci2ecef();
        if(this.use_moon || this.use_sun){
            String dir_in;
//            String fs = FileUtil.file_separator();
            dir_in = "C:/Code/Jat/jat/eph/DE405data/";
            jpl_ephem = new DE405(dir_in);
        }
        if(this.use_sun) compute_JPL_Sun_Vector();
        if(this.use_moon) compute_JPL_Moon_Vector();
    }

    
    /** Set the IERS (Earth Rotation) Data
     * @param x x_pole coordinate (rad)
     * @param y y_pole coordinate (rad)
     * @param d UT1 - UTC difference in s.
     */
    public void setIERS(double x, double y, double d){
        this.x_pole = x*Constants.arcsec2rad;
        this.y_pole = y*Constants.arcsec2rad;
        if(sim_time >= n){
            this.UT1_UTC = d;
            n = n+increment;
        }else{
            this.UT1_UTC = 0.0;
        }
        this.MJD_UT1 = this.MJD_UTC + this.UT1_UTC/86400.0;
        //this.T = trueOfDate();
        //this.E = eci2ecef();
    }
    
    /** Returns the Modified Julian Date of the start of the simulation
     * @return modified julian date of the start of the simulation in UTC
     */
    public double MJD_UTC_START(){
        return this.MJD_UTC_START;
    }
        
    /** Returns the ECI to ECEF (ICRF to ITRF) Transformation Matrix.
     * @return ECI to ECEF Transformation Matrix.
     */    
    public Matrix ECI2ECEF(){
        return this.E;
    }

    /** Returns the J2000 to True of Date (ICRF to TOD) Transformation Matrix.
     * @return J2000 to TOD Transformation Matrix.
     */        
    public Matrix TOD(){
        return this.T;
    }
    
    /** Return the current Terrestrial time in MJD format.
     * @return MJD of current Terrestrial time.
     */
    public double mjd_tt(){
        return this.MJD_TT;
    }
    
    public double mjd_tdb(){
        return this.MJD_TDB;
    }
    /** Return the current UTC time in MJD format.
     * @return MJD of current UTC time.
     */
    public double mjd_utc(){
        return this.MJD_UTC;
    }
    
    /** Return the current UTC time in JD format.
     * @return JD of current UTC time.
     */
    public double jd_utc(){
        return this.MJD_UTC+2400000.5;
    }
    
    /** Return the current TT time in JD format.
     * @return JD of current TT time.
     */
    public double jd_tt(){
        return this.MJD_TT+2400000.5;
    }
    
    /** Return the current UT1 time in MJD format.
     * @return MJD of current UT1 time.
     */
    public double mjd_ut1(){
        return this.MJD_UT1;
    }
    
    /** Computes the mean obliquity of the ecliptic. Uses Mjd_TT (Terrestrial Time).
     *@return  Mean obliquity of the ecliptic
     */
    public double MeanObliquity() {
        double T = (this.MJD_TT - MJD_J2000)/36525.0;
        //double T = (this.MJD_TDB - MJD_J2000)/36525.0; //* used prior to IERS1996 convention
//      TODO
        return MathUtils.DEG2RAD *( 23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
        //* Debug - updated the epsilon-0 value according to the Supplemental Astro Almanac (old above)
        //return MathUtils.DEG2RAD *( 23.4392911111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
        //return Constants.arcsec2rad *( 84381.448-(46.8150+(0.00059-0.001813*T)*T)*T );
    }
    
    /** Transformation of equatorial to ecliptical coordinates. Uses Mjd_TT (Terrestrial Time).
     *   @return  Transformation matrix
     */
    public RotationMatrix EclMatrix() {
        RotationMatrix out = new RotationMatrix(1, MeanObliquity());
        return out;
    }
    
    /** Precession transformation of equatorial coordinates. Uses Mjd_TT (Terrestrial Time).
     *  @return  Precession transformation matrix
     */
    public Matrix PrecMatrix() {
        double Mjd_1 = MJD_J2000;  // Epoch given (Modified Julian Date TT)
        double Mjd_2 = this.MJD_TT;   // Epoch to precess to (Modified Julian Date TT)
        //double Mjd_2 = this.MJD_TDB;   //* used prior to IERS1996 convention
        
        // Constants
        final double Arcs = MathUtils.ARCSEC2RAD;
        final double T  = (Mjd_1-MJD_J2000)/36525.0;
        final double dT = (Mjd_2-Mjd_1)/36525.0;
        
        // Variables
        double zeta,z,theta;
        
        // Precession angles
        zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+
        ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT*Arcs;
        z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT*Arcs; 
        theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-
        ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT*Arcs;
        
        //***Debug - Vallado's method
//        double Deg2Rad = 0.01745329251994;
//        double TTDB = dT;
//        double TTDB2 = TTDB*dT;
//        double TTDB3 = TTDB2*dT;
//        zeta = 0.6406161 * TTDB + 0.0000839 * TTDB2 + 0.0000050 * TTDB3;
//        z = 0.6406161 * TTDB + 0.0003041 * TTDB2 + 0.0000051 * TTDB3;
//        theta = 0.5567530 * TTDB - 0.0001185 * TTDB2 - 0.0000116 * TTDB3;
//        zeta *= Deg2Rad;
//        z *= Deg2Rad;
//        theta *= Deg2Rad;


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
    public void NutAngles() {
        double Mjd_TT = this.MJD_TT;
        //double Mjd_TT = this.MJD_TDB; //* used prior to IERS1996 convention
        
        // Constants
        
        final double T  = (Mjd_TT-MJD_J2000)/36525.0;
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

        //*** DEBUG - Vallado
//        double rr = 360.0;   /*deg*/
//        double TTdb = T;
//        double Ttdb2 = T2;
//        double Ttdb3 = T3;
//        double Deg2Rad = 0.01745329251994;
//        l = 134.9629814 + (1325 * rr + 198.8673981) * TTdb + 0.0086972 * Ttdb2 +
//            0.00001778 * Ttdb3;
//        lp = 357.5277233 +
//             (99 * rr + 359.05034) * TTdb - 0.00016028 * Ttdb2 - 0.00000333 * Ttdb3;
//        F = 93.2719103 + (1342 * rr + 82.0175381) * TTdb - 0.0036825 * Ttdb2 +
//            0.00000306 * Ttdb3;
//        D = 297.8503631 + (1236 * rr + 307.111480) * TTdb - 0.00191417 * Ttdb2 +
//            0.00000528 * Ttdb3;
//        Om = 125.0445222 - (5 * rr + 134.1362608) * TTdb + 0.0020708 * Ttdb2 +
//      	  0.00000222 * Ttdb3;
//        l = MathUtils.Modulo(l, 360.0) * Deg2Rad;
//        lp = MathUtils.Modulo(lp, 360.0) * Deg2Rad;
//        F = MathUtils.Modulo(F, 360.0) * Deg2Rad;
//        D = MathUtils.Modulo(D, 360.0) * Deg2Rad;
//        Om = MathUtils.Modulo(Om, 360.0) * Deg2Rad;
//        
//        // Nutation in longitude and obliquity [rad]
//        final double Arcs = MathUtils.ARCSEC2RAD;
//        
//        deps = 0.0;
//        dpsi = 0.0;
//        for (int i=0; i<N_coeff; i++) {
//            arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om );            
//            dpsi += ( C[i][5]+C[i][6]*T ) * Math.sin(arg);
//            deps += ( C[i][7]+C[i][8]*T ) * Math.cos(arg);
//        }
//        
//        // TODO
//        //* Delta Delta corrections
////        double ddpsi = -0.0534;
////        double ddeps = -0.00472;
////        this.dpsi = this.dpsi + ddpsi;
////        this.deps = this.deps + ddeps;
//        
//        this.dpsi = 1.0E-5 * dpsi*Arcs;
//        this.deps = 1.0E-5 * deps*Arcs;
//        this.Om = Om;
        
    }
    
    /** Transformation from mean to true equator and equinox. Uses Mjd_TT (Terrestrial Time).
     *   @return  Nutation matrix
     */
    public Matrix NutMatrix() {
        double Mjd_TT = this.MJD_TT;
        //double Mjd_TT = this.MJD_TDB; //* used prior to IERS1996 convention
        
        // Mean obliquity of the ecliptic
        double eps = MeanObliquity();
        
        // Nutation in longitude and obliquity
        NutAngles();
        
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
//    public Matrix NutMatrixSimple() {
//        double Mjd_TT = this.MJD_TT;
//        
//        // Constants
//        final double Arcs = 1.0/MathUtils.ARCSEC2RAD;
//        final double T  = (Mjd_TT-MJD_J2000)/36525.0;
//        
//        // Variables
//        double  ls, D, F, N;
//        double  eps;
//        
//        // Mean arguments of luni-solar motion
//        ls = pi2*MathUtils.Frac(0.993133+  99.997306*T);   // mean anomaly Sun
//        D  = pi2*MathUtils.Frac(0.827362+1236.853087*T);   // diff. longitude Moon-Sun
//        F  = pi2*MathUtils.Frac(0.259089+1342.227826*T);   // mean argument of latitude
//        N  = pi2*MathUtils.Frac(0.347346-   5.372447*T);   // longit. ascending node
//        
//        // Nutation angles
//        dpsi = ( -17.200*Math.sin(N)   - 1.319*Math.sin(2*(F-D+N)) - 0.227*Math.sin(2*(F+N))
//        + 0.206*Math.sin(2*N) + 0.143*Math.sin(ls) ) / Arcs;
//        deps = ( + 9.203*Math.cos(N)   + 0.574*Math.cos(2*(F-D+N)) + 0.098*Math.cos(2*(F+N))
//        - 0.090*Math.cos(2*N)                 ) / Arcs;
//        
//        // Mean obliquity of the ecliptic
//        eps  = 0.4090928-2.2696E-4*T;
//        
//        RotationMatrix r1 = new RotationMatrix(1, (-eps-deps));
//        RotationMatrix r2 = new RotationMatrix(3, -dpsi);
//        RotationMatrix r3 = new RotationMatrix(1, eps);
//        
//        Matrix out = r1.times(r2.times(r3));
//        return  out;
//    }
    
    /** Computation of the equation of the equinoxes.
     * Uses Mjd_TT (Terrestrial Time).
     * Notes: The equation of the equinoxes dpsi*Math.cos(eps) is the right ascension of the
     * mean equinox referred to the true equator and equinox and is equal to the
     * difference between apparent and mean sidereal time.
     * @return Equation of the equinoxes
     */
    public double EqnEquinox() {
        // Nutation in longitude and obliquity
        NutAngles();
        double arcs = Constants.arcsec2rad;
        // Equation of the equinoxes
//      TODO
//        return  dpsi * Math.cos( MeanObliquity() );
        //* Montenbruck
//        return  dpsi * Math.cos( MeanObliquity() ) 
//        			+ (0.002649*Math.sin(Om)-0.000013*Math.cos(Om))*arcs;
        //* Astro Almenac 96 modified
        return  dpsi * Math.cos( MeanObliquity() )
        			+ (0.002649*Math.sin(Om)+0.000063*Math.sin(2*Om))*arcs;
        //* Astro Almenac 96
//        return  dpsi * Math.cos( MeanObliquity() )
//					+ (0.00264*Math.sin(Om)+0.000063*Math.sin(2*Om))*arcs;
    }
    
    /** Greenwich Mean Sidereal Time. Uses Mjd_UT1.
     *   @return  GMST in [rad]
     */
    public double GMST() {
        double Mjd_UT1 = this.MJD_UT1;
        
        // Constants
        final double Secs = 86400.0;        // Seconds per day
        
        // Variables
        double Mjd_0,UT1,T_0,T,gmst;
        
        // Mean Sidereal Time
        Mjd_0 = Math.floor(Mjd_UT1);
        UT1   = Secs*(Mjd_UT1-Mjd_0);          // [s]
        T_0   = (Mjd_0  -MJD_J2000)/36525.0;
        T     = (Mjd_UT1-MJD_J2000)/36525.0;
        double T2 = T*T;
        
        gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1
        + (0.093104-6.2e-6*T)*T*T; // [s]
        double tmp = pi2*MathUtils.Frac(gmst/Secs);       // [rad], 0..2pi
        if(tmp < 0) tmp = tmp + pi2;
        return tmp;
      //TODO
        //* debug below IERS 1996
//        double r = 1.002737909350795 + (5.9006e-11 - 5.9e-15 * T_0)*T_0;
//        gmst  = 24110.54841 + 8640184.812866*T_0 + (0.093104-6.2e-6*T_0)*T_0*T_0 + r*UT1; 
//        //double tmp =  gmst*Constants.omega_e;
//        double tmp;
//        tmp = pi2*MathUtils.Frac(gmst/Secs);
//        //tmp = MathUtils.Modulo(tmp,pi2);
//        return tmp;
        
    }
    
    /** Greenwich Apparent Sidereal Time.
     *   @return  GAST in [rad]
     */
    public double GAST() {
        return MathUtils.Modulo( GMST() + EqnEquinox(), pi2 );
    }
    
    /** Transformation from true equator and equinox to Earth equator and
     * Greenwich meridian system.
     * @return Greenwich Hour Angle matrix
     */
    public RotationMatrix GHAMatrix() {
        RotationMatrix out = new RotationMatrix(3, GAST());
        return  out;
    }
    
    /** Transformation from pseudo Earth-fixed to Earth-fixed coordinates
     * for a given date. Uses Mjd_UTC if poles are a function of UTC.
     * @return Pole matrix
     */
    public Matrix PoleMatrix() {
        RotationMatrix r1 = new RotationMatrix(2, -this.x_pole);
        RotationMatrix r2 = new RotationMatrix(1, -this.y_pole);
        Matrix out = r1.times(r2);
//        double[][] tmp = {{1,0,this.x_pole},
//                		  {0,1,-this.y_pole},
//                		  {-this.x_pole,this.y_pole,1}};
//        Matrix out = new Matrix(tmp);
        return  out;
    }
    
//    public Matrix PoleMatrix() {
//        Matrix out = new Matrix(3);
//        out.set(0,2,x_pole);
//        out.set(1,2, -y_pole);
//        out.set(2,0, -x_pole);
//        out.set(2,1, y_pole);
//        return  out;
//    }
    
    /** J2000 to TOD Transformation
     * @return J2000 to ECEF transformation matrix
     */
    
    public Matrix trueOfDate(){
        Matrix P = PrecMatrix();
        Matrix N = NutMatrix();
        Matrix out = N.times(P);
        return out;
    }
    
    /** ECI to ECEF Transformation
     * @return ECI to ECEF transformation matrix
     */
    public Matrix eci2ecef(){
        Matrix T = trueOfDate();
        Matrix G = GHAMatrix();
        //Matrix Pole = PoleMatrix();
        Matrix Pole = new Matrix(3);
        Matrix A = Pole.times(G);
        Matrix E = A.times(T);
        return E;
    }
        
    /** Increment time to allow the computation of reference frames at a new time
     * @param sec Increment to move time in seconds. Can be + or -.
     */
    public void incrementTime(double sec){
        double frac = sec/86400.0;
        this.MJD_TT = this.MJD_TT + frac;
        this.MJD_TDB = TimeUtils.TTtoTDB(this.MJD_TT);
        this.MJD_UTC = this.MJD_UTC + frac;
        this.MJD_UT1 = this.MJD_UT1 + frac;
        this.T = trueOfDate();
        this.E = eci2ecef();
        this.compute_JPL_Sun_Vector();
    }
    
    /** Updates the Earth model a number of seconds since the initial epoch.
     * Since the RungeKutta algorithm requires the ability to look at an 
     * arbitrary point in time, it is necessary to allow arbitrary stepping
     * of time in the Earth model.  
     * 
     * @param sec: seconds since MJD_UTC_START
     */
    public void updateTimeSinceStart(double sec){
        this.sim_time = sec;
        this.MJD_UTC = this.MJD_UTC_START+sec/86400;
        this.MJD_TT = CalDate.UTC2TT(this.MJD_UTC);
        this.MJD_TDB = TimeUtils.TTtoTDB(this.MJD_TT);
        this.MJD_UT1 = this.MJD_UTC + this.UT1_UTC/86400.0;
        this.T = trueOfDate();
        this.E = eci2ecef();
        if(this.use_sun)  compute_JPL_Sun_Vector();
        if(this.use_moon) compute_JPL_Moon_Vector();
    }
    
    public double get_omega_e(){
//      Variables
        double Mjd_0;
        
        // Mean Sidereal Time
        Mjd_0 = Math.floor(this.MJD_UT1);
        double Tu   = (Mjd_0  -MJD_J2000)/36525.0;
        return 7292115.8553e-11 + 4.3e-15 * Tu;
    }
    
    /** Computes the Sun's geocentric position using a low precision analytical series.
     * @return Solar position vector [m] with respect to the mean equator and equinox of J2000 (EME2000, ICRF)
     */
    public VectorN sunVector() {
        double Mjd_TT = this.MJD_TT;
        // Constants
        
        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
        double T   = (Mjd_TT-MJD_J2000)/36525.0;  // Julian cent. since J2000
        
        
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
        // We assume (hope) that the passed in time is the current
        // EarthRef time.
    	
    	// TODO: not using the given Time is an interface violation - this is dangerous
    	
        if(this.use_sun)
            return r_sun;
        else
            return new VectorN(3);
    }
    
    /**
     * Get the geocentric position of the moon [m] in the J2000 (DE405) inertial frame
     * @return Vector 3 [m]
     */
    public VectorN get_JPL_Moon_Vector(){
        if(this.use_moon)
            return r_moon;
        else
            return new VectorN(3);
    }
    
    /**
     * Compute the current JPL vector to the Sun [m]
     */
    private void compute_JPL_Sun_Vector(){
//        r_sun = new VectorN(jpl_ephem.get_Geocentric_Sun_pos(this.jd_tt()));
//        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
//        RotationMatrix R = new RotationMatrix(1, -eps);
//        r_sun = R.times(r_sun);
        //double jd_tdb = Time.TTtoTDB(this.mjd_tt())+2400000.5;
        
        //* July 7 2005 - buggy
//        jpl_ephem.planetary_ephemeris(this.jd_utc());
//        VectorN r_suntmp = jpl_ephem.get_pos(DE405.SUN,this.jd_utc());
//        VectorN r_earth = jpl_ephem.get_pos(DE405.EARTH,this.jd_utc());
//        VectorN r_body = r_suntmp.minus(r_earth);
//        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
//        RotationMatrix R = new RotationMatrix(1, -eps);
//        r_sun = R.times(r_body);
//        r_sun = r_sun.times(1000);
        
        //* June 24 2005 - Works
        if(this.use_sun){
            r_sun = new VectorN(jpl_ephem.get_planet_posvel(DE405_Body.GEOCENTRIC_SUN, this.MJD_TT));
        }
    }
    
    /**
     * Compute the current JPL vector to the Moon [m]
     */
    private void compute_JPL_Moon_Vector(){
        if(this.use_moon){
            r_moon = new VectorN(jpl_ephem.get_planet_posvel(DE405_Body.GEOCENTRIC_MOON, this.MJD_TT));
        }
    }
    
    public void set_use_sun(boolean b){
        this.use_sun = b;
        if(b){
            if(jpl_ephem == null){
                String fs = FileUtil.file_separator();
                String dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
                jpl_ephem = new DE405(dir_in);
            }
            this.compute_JPL_Sun_Vector();
        }
    }
    
    public void set_use_moon(boolean b){
        this.use_moon = b;
        if(b){
            if(jpl_ephem == null){
                String fs = FileUtil.file_separator();
                String dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs;
                jpl_ephem = new DE405(dir_in);
            }
            this.compute_JPL_Moon_Vector();
        }
    }
    
    /** Computes the Moon's geocentric position using a low precision analytical series.
     * @return Lunar position vector [m] with respect to the mean equator and equinox of J2000 (EME2000, ICRF).
     */
    public VectorN moonVector() {
        double Mjd_TT = this.MJD_TT;
        
        double eps = Constants.eps*MathUtils.DEG2RAD;             // Obliquity of J2000 ecliptic
        double T   = (Mjd_TT-MJD_J2000)/36525.0;  // Julian cent. since J2000
        
        
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
     * Converts the MJD expressed in UTC to give the days since Jan 01 00:00
     * @return The number of days since the beginning of the current year.
     */
    public int dayOfYear(){
        GPSTimeFormat time = new GPSTimeFormat(mjd_utc());
        CalDate date = time.calDate();
        return date.doy();
    }
    /**
     * Converts the MJD expressed in UTC to give the seconds since 00:00 UTC
     * @return The number of seconds since the beginning of the current day (UTC).
     */
    public double secOfDay(){
        GPSTimeFormat time = new GPSTimeFormat(mjd_utc());
        CalDate date = time.calDate();
        return date.sec_of_day();
    }
    
    public double get_sim_time(){
        return this.sim_time; // [s]
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#get_spin_rate()
     */
    public double get_spin_rate(Time t) {
        // TODO Auto-generated method stub
        return get_omega_e();
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#get_mean_radius()
     */
    public double get_mean_radius() {
        // TODO Auto-generated method stub
        return Constants.R_Earth;
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#get_grav_const()
     */
    public double get_grav_const() {
        // TODO Auto-generated method stub
        return Constants.GM_Earth;
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#inertial_to_body(double, double)
     */
    public Matrix inertial_to_body(Time t) {
        // TODO Auto-generated method stub
        return eci2ecef();
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#body_to_inertial(double, double)
     */
    public Matrix body_to_inertial(Time t) {
        // TODO Auto-generated method stub
        return eci2ecef().transpose();
    }

    /* (non-Javadoc)
     * @see jat.spacetime.BodyRef#trueOfDate(jat.spacetime.Time)
     */
    public Matrix trueOfDate(Time t) {
        // TODO Auto-generated method stub
        return trueOfDate();
    }
   
    /**
     * Returns a translater to translate into other reference frames.
     * @param other another reference frame
     * @param t time at which translation will be done
     * @return translater object or null if does not know how
     * to translate
     */
    public ReferenceFrameTranslater getTranslater(ReferenceFrame other, Time t)
    {
      // Currently does not translate anything (but itself).
      return (other instanceof EarthRef ?
          new ReferenceFrameTranslater() : null);
    }

}
