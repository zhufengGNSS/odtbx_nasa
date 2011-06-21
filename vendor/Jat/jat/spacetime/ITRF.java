package jat.spacetime;
/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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

import jat.math.*;
import jat.matvec.data.*;
/**
 * ITRF provides a representation of the Internation Terrestrial Reference Frame.
 * It is used to compute the transformation between ECI and ECEF frames.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ITRF {
	private double omega = 0.0;
	private double dpsi = 0.0;
	private double deps = 0.0;
	private Matrix nutMatrix;
	private Matrix precMatrix;
	private Matrix poleMatrix;
	private Matrix gstMatrix;
	private double mjd_tt;
	private double mjd_ut1;
	private double _xp;
	private double _yp;
	private static final double JAN11997 = TimeUtils.JDtoMJD(2450449.5);
	private double t;
	private double t2;
	private double t3;
	private double t4;
	private double eps;
	
	/**
	 * Create an instance of ITRF
	 * @param time a Time Object
	 * @param xp double containing the x-pole location in arcsec (from IERS Bulletin B)
	 * @param yp double containing the y-pole location in arcsec (from IERS Bulletin B)
	 * @param ut1_utc double containing the difference between UT1 and UTC in seconds (from IERS Bulletin B)
	 */
	public ITRF(Time time, double xp, double yp, double ut1_utc){
		// convert polar motion data to radians
		_xp = xp * MathUtils.ARCSEC2RAD;
		_yp = yp * MathUtils.ARCSEC2RAD;
//		System.out.println("xp = "+_xp+" yp = "+_yp);
		// get the Terrestrial Time
		mjd_tt = time.mjd_tt();
		time.set_UT1_UTC(ut1_utc);
		t = centuriesSinceJ2000(mjd_tt);
		t2 = t*t;
		t3 = t2*t;
		t4 = t3*t;
		mjd_ut1 = time.mjd_ut1();
		
		// make the computations
    	this.precession();
    	this.nutation();
    	this.sidereal();
    	this.polarMotion();		
	}
	
	/**
	 * Return Julian centuries since J2000 (Terrestrial Time)
	 * @param mjd_tt double containing TT in MJD format
	 * @return Julian centuries since J2000 (Terrestrial Time)
	 */
	private static double centuriesSinceJ2000(double mjd) {
		return (mjd - TimeUtils.MJD_J2000) / 36525.0;
	}
	
	/**
	 * Return Greenwich Mean Sidereal Time in radians
	 * @return GMST in radians
	 */
	private double GMST(){
    	// compute time since J2000 in days
		double T = mjd_ut1 - TimeUtils.MJD_J2000;
		// compute GMST in radians (CSR-02-01, p.21, MSODP Implementation)
		double gmst = 4.894961212823058751375704430 + T*(6.30038809898489355227651372 + T*(5.075209994113591478053805523E-15 - 9.25309756819433560067190688E-24*T));
		// quadrant check (return a value between 0 and 2PI)
		gmst = MathUtils.radiancheck(gmst);
//		System.out.println("gmst = "+gmst);
		return gmst;
	}
	
	/**
	 * Return mean obliquity in radians
	 * @return GMST in radians
	 */
	private void meanObliquity (){
		// compute meanObliquity in arcsec
		eps = 84381.448 - 46.8150*t - 0.00059*t2 + 0.001813*t3;
		// convert to radians
		eps = eps*MathUtils.ARCSEC2RAD;
		eps = MathUtils.radiancheck(eps);
//		System.out.println("eps = "+eps);
	}
	
    /** Computes Nutation in longitude and obliquity using the IAU 1980 nutation theory.
     */
    private void nutation() {
        meanObliquity();
        // Constants
                
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
        double  l, lp, F, D;
        double  arg;
        
        // Mean arguments of luni-solar motion in degrees
        //
        //   l   mean anomaly of the Moon
        //   l'  mean anomaly of the Sun
        //   F   mean argument of latitude
        //   D   mean longitude elongation of the Moon from the Sun
        //   omega  mean longitude of the ascending node
        
        l  = 134.96340251 + (1717915923.2178*t + 31.8792*t2 + 0.051635*t3 - 0.0002447*t4)/3600.0; 
        lp = 357.52910918 + (129596581.0481*t - 0.5532*t2 + 0.000136*t3 - 0.00001149*t4)/3600.0;
        F  = 93.27209062 + (1739527262.8478*t - 12.7512*t2 - 0.001037*t3 + 0.00000417*t4)/3600.0;
        D  = 297.85019547 + (1602961601.2090*t - 6.3706*t2 + 0.006593*t3 - 0.00003169*t4)/3600.0;
        omega = 125.04455501 + (-6962890.2665*t + 7.4722*t2 + 0.007702*t3 - 0.00005939*t4)/3600.0;
        
        // Quadrant checks and convert to radians
        l = MathUtils.degreecheck(l)*MathUtils.DEG2RAD;
        lp = MathUtils.degreecheck(lp)*MathUtils.DEG2RAD;
        F = MathUtils.degreecheck(F)*MathUtils.DEG2RAD;       
        D = MathUtils.degreecheck(D)*MathUtils.DEG2RAD;
        omega = MathUtils.degreecheck(omega)*MathUtils.DEG2RAD;
               
        // Nutation in longitude and obliquity     
        for (int i=0; i<N_coeff; i++) {
            arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*omega );
            dpsi += ( C[i][5]+C[i][6]*t ) * Math.sin(arg);
            deps += ( C[i][7]+C[i][8]*t ) * Math.cos(arg);
        }

        // convert to arcsec (1.0E-05) then to radians        
        dpsi = 1.0E-5 * dpsi * MathUtils.ARCSEC2RAD;
        deps = 1.0E-5 * deps * MathUtils.ARCSEC2RAD;
//        System.out.println("dpsi = "+dpsi+" deps = "+deps);
        
        // compute nutation matrix
        RotationMatrix r1 = new RotationMatrix(1, -eps);
        RotationMatrix r2 = new RotationMatrix(3, dpsi);
        RotationMatrix r3 = new RotationMatrix(1, eps+deps);        
        nutMatrix = r1.times(r2).times(r3);
//        nutMatrix.print("nut");
    }
    /** Computes precession matrix
     *
     * @return precession matrix
     */
    private void precession(){
    	double zeta = MathUtils.ARCSEC2RAD*(2306.2181*t + 0.30188*t2 + 0.017998*t3);
    	double theta = MathUtils.ARCSEC2RAD*(2004.3109*t - 0.42665*t2 - 0.041833*t3);
    	double z = MathUtils.ARCSEC2RAD*(2306.2181*t + 1.09468*t2 + 0.018203*t3);
    	
    	zeta = MathUtils.radiancheck(zeta);
    	theta = MathUtils.radiancheck(theta);
    	z = MathUtils.radiancheck(z);
    	
//    	System.out.println("zeta = "+zeta+" theta = "+theta+" z = "+z);
        RotationMatrix r1 = new RotationMatrix(3, zeta);
        RotationMatrix r2 = new RotationMatrix(2, -theta);
        RotationMatrix r3 = new RotationMatrix(3, z);
        
        precMatrix = r1.times(r2.times(r3));
//        precMatrix.print("prec");
    }
	/**
	 * Computes Greenwich Sidereal Time in radians
	 * Note: nutation must be called before this is called
	 */    
    private double GST(){
    	double gmst = GMST();
    	meanObliquity();
    	double eqnEquinoxes = dpsi * Math.cos(eps);
    	double gst = gmst + eqnEquinoxes;
    	// if date is after Jan 1, 1997 use extra terms
    	if (mjd_tt > JAN11997){
    		double term1 = 0.00264 * Math.sin(omega);
    		double term2 = 0.000063 * Math.sin(2.0 * omega);
    		gst = gst + (term1 + term2)*MathUtils.ARCSEC2RAD;
    	}
    	// quadrant check
    	gst = MathUtils.radiancheck(gst);
//    	System.out.println("gst = "+gst);
    	return gst;
    }
    /**
     * Computes the Polar Motion Matrix
     */
    private void polarMotion(){
    	RotationMatrix r1 = new RotationMatrix(1, _yp);
    	RotationMatrix r2 = new RotationMatrix(2, _xp);
    	poleMatrix = r1.times(r2);
//    	poleMatrix.print("pole");
    }
    /**
     * Computes the Matrix to account for Sidereal Time 
     */
    private void sidereal(){
    	double gst = this.GST();
    	gstMatrix = new RotationMatrix(3, -gst);
//    	gstMatrix.print("gstMatrix");
    }
    /**
     * Computes the ECEF to ECI Transformation Matrix
     * @return Matrix containing the ECEF to ECI transformation
     */
    public Matrix ECEF2ECI(){
    	return precMatrix.times(nutMatrix.times(gstMatrix.times(poleMatrix)));
    }
    /**
     * Computes the ECI to ECEF Transformation Matrix
     * @return Matrix containing the ECI to ECEF transformation
     */   
    public Matrix ECI2ECEF() {
    	return ECEF2ECI().transpose();
    }
	 

	/**
	 * Runs the CSR-02-01 Test Case
	 * @param args
	 */
	public static void main(String[] args) {
		Time t = new Time(2001, 11, 26, 18, 30, 24.50);
		t.set_UT1_UTC(-0.088405058232);
		double xpole = -5.299915237216626E-07 * MathUtils.RAD2DEG * 3600.0;
		double ypole = 8.9015513447899459E-07 * MathUtils.RAD2DEG * 3600.0;
		ITRF x = new ITRF(t, xpole, ypole, -0.088405058232);
		Matrix ecef2eci = x.ECEF2ECI();
		ecef2eci.print("ecef2eci");
		Matrix eci2ecef = x.ECI2ECEF();
		eci2ecef.print("eci2ecef");
		
//		Time t2 = new Time(2004, 6, 1, 12, 0, 0.0);
//		ITRF test = new ITRF(t2, 0.2251760, 0.3927010, -0.0256666);
//		eci2ecef = test.ECI2ECEF();
//		VectorN iss_eci = new VectorN(-4453.783586, -5038.203756, -426.384456);
//		VectorN iss_ecef = eci2ecef.times(iss_eci);
//		iss_ecef.print("iss_ecef");
//		VectorN ss_eci = new VectorN(-2290.301063, -6379.471940, 0.0);
//		VectorN ss_ecef = eci2ecef.times(ss_eci);
//		ss_ecef.print("ss_ecef");
//		VectorN gps_eci = new VectorN(5525.33668, -15871.18494, -20998.992446);
//		VectorN gps_ecef = eci2ecef.times(gps_eci);
//		gps_ecef.print("gps_ecef");
//		VectorN mol_eci = new VectorN(-1529.894287, -2672.877357, -6150.115340);
//		VectorN mol_ecef = eci2ecef.times(mol_eci);
//		mol_ecef.print("mol_ecef");
//		VectorN geo_eci = new VectorN(36607.358256, -20921.723703, 0.0);
//		VectorN geo_ecef = eci2ecef.times(geo_eci);
//		geo_ecef.print("geo_ecef");	
	}

}
