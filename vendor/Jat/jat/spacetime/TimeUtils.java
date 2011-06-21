/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2003 United States Government as represented by the
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
 */

package jat.spacetime;

import jat.math.MathUtils;
import jat.util.FileUtil;

/**
 * A collection of utility functions for converting between various time 
 * scales and basic calendar arithmetic. 
 * Conversions involving UTC depend on a data file of leap seconds.
 */
public class TimeUtils {
	
    /**
     * Modified Julian Date of the J2000 Epoch in Terrestrial Time (TT).
     * (J2000 = Jan 1, 2000, 12:00 TT)
     */
    public static final double MJD_J2000 = 51544.5;
    
    /**
     * Fraction of an SI day per SI second.
     */
    public static final double sec2days = 1.0/86400;
    
    /**
     * Number of SI seconds per SI day.
     */
    public static final double days2sec = 86400;
    
    /**
     * Constant used for conversion to Terrestrial Time (TT) from
     * International Atomic Time (TAI), in SI seconds:
     *  TT = TAI + TT_TAI.
     */
    public static final double TT_TAI = 32.184;  // constant

    /**
     * Constant used for conversion to GPS time from TAI, in
     * SI seconds:
     * TAI = GPS + TAI_GPS.
     */
    public static final int TAI_GPS = 19;  // constant
    
    /**
     * The UTC-TAI leap second table that underlies this class.
     * Note, this is static, which means that there is one set of leap seconds
     * for the TimeUtils class.
     */
    protected static LeapSecondTable leapsectbl;
    
    /**
     * Static initialization block.
     */
    static {
		String directory;
		
    	try {
		    directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("TimeUtils could not load its leap second data file, leapsec.dat.");
			System.err.println("This means all JAT UTC leap second conversions are incorrect.");
			System.err.println("DO NOT use JAT until this error has been corrected!");
		}
    }
    
    /**
     * Determines if the given year is a leap year or a common year.  
     * Applies all the leap year rules for years divisible by "4", 
     * "100", and "400", which are required for proper "century" handling 
     * (e.g. 1700 (commmon), 1800 (commmon), 1900 (commmon), 2000 (leap), 
     * 2100 (commmon), 2200 (commmon), etc.).
     * @param year The year (assumed Gregorian Calendar)
     * @return True if it is a leap year, false if is a common year (not a leap year).
     */
    public static boolean isLeapYear(int year) {
    	/* Leap year rules:
    	 * http://en.wikipedia.org/wiki/Leap_year
    	 * 1) "Most" years evenly divisible by 4,
    	 * 2) Except those years evenly divisible by 100
    	 * 3) However, years evenly divisible by 400 are leap years, "century leap years"
    	 * Doesn't include John Herschel's recommendation for 4000 not being a 
    	 * leap year, which is standard practice.
    	 */
    	if ( ( ((year%4) == 0) && ((year%100) != 0) ) || ((year%400) == 0) ){
    		return true;
    	}
    	else {
    		return false;
    	}
    }

    /** Convert Day of Month to Day of Year.
     * @param year Year.
     * @param month Month.
     * @param mday Day of month.
     * @return Day of year.
     */
    public static int day2doy(int year, int month, int mday){
        int [] regu_month_day = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
        int [] leap_month_day = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
        int yday = 0;

        // check for leap year
        if (isLeapYear(year)) {
            yday = leap_month_day[month - 1] + mday;
        }
        else {
            yday = regu_month_day[month - 1] + mday;
        }
        return yday;
    }

    /** Compute the month from day of year. From GPS Toolkit code.
     * @param year Year.
     * @param doy Day of year.
     * @return Month.
     */
    public static int doy2month(int year, int doy){
        int [] regu_month_day = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
        int [] leap_month_day = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};

        int month = 0;
//        int mday = 0;

        int guess = (int)(doy*0.032);
        int more = 0;

        // check for leap year
        if (isLeapYear(year)) {
            if ((doy - leap_month_day[guess+1]) > 0) more = 1;
            month = guess + more + 1;
//            mday = doy - leap_month_day[guess+more];
        }
        else {
            if ((doy - regu_month_day[guess+1]) > 0) more = 1;
            month = guess + more + 1;
//            mday = doy - regu_month_day[guess+more];
        }
        return month;
    }

    /** Compute the day of month from day of year. From GPS Toolkit code.
     * @param year Year.
     * @param doy Day of year.
     * @return Day of month.
     */
    public static int doy2day(int year, int doy){
        int [] regu_month_day = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
        int [] leap_month_day = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};

//        int month = 0;
        int mday = 0;

        int guess = (int)(doy*0.032);
        int more = 0;

        // check for leap year
        if (isLeapYear(year)) {
            if ((doy - leap_month_day[guess+1]) > 0) more = 1;
//            month = guess + more + 1;
            mday = doy - leap_month_day[guess+more];
        }
        else {
            if ((doy - regu_month_day[guess+1]) > 0) more = 1;
//            month = guess + more + 1;
            mday = doy - regu_month_day[guess+more];
        }
        return mday;
    }
	
	
    /** Return the difference between TAI and UTC (known as accumulated leap 
     * seconds), where:<br>
     * TAI = UTC + accumLeapSecs.<br>
     * Returns zero if before the first leap second. Holds the value of the 
     * last specified leap second for dates after the last leap second.<br>
     * 
     * Note, TAI is a continuous time but UTC is not.  When crossing leap
     * seconds the same UTC time can result from different TAI times, as in:<br>
     * (at 0.25 sec steps, all times MJD)<br>
     * 42048.0001244213 TAI,   42047.999985532406 UTC, dATutc=12, dATtai=12<br>
     * 42048.000127314815 TAI, 42047.99998842592  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013020833 TAI,  42047.99999131944  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013310185 TAI,  42047.99999421296  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013599537 TAI,  42047.999997106475 UTC, dATutc=12, dATtai=12<br>
     * (begin leap second)<br>
     * 42048.00013888889 TAI,  42048.0            UTC, dATutc=13, dATtai=12<br>
     * 42048.00014178241 TAI,  42048.00000289352  UTC, dATutc=13, dATtai=12<br>
     * 42048.00014467593 TAI,  42048.000005787035 UTC, dATutc=13, dATtai=12<br>
     * 42048.000147569444 TAI, 42048.00000868055  UTC, dATutc=13, dATtai=12<br>
     * (end leap second)<br>
     * 42048.00015046296 TAI,  42048.0            UTC, dATutc=13, dATtai=13<br>
     * 42048.00015335648 TAI,  42048.00000289352  UTC, dATutc=13, dATtai=13<br>
     * 42048.00015625 TAI,     42048.000005787035 UTC, dATutc=13, dATtai=13<br>
     * 42048.00015914352 TAI,  42048.00000868056  UTC, dATutc=13, dATtai=13<br>
     * 42048.00016203704 TAI,  42048.00001157408  UTC, dATutc=13, dATtai=13<br>
     * 42048.000164930556 TAI, 42048.000014467594 UTC, dATutc=13, dATtai=13<br>
     * where dATutc is the results of tai_utc(UTC) and <br>
     * dATtai = tai_utc_from_tai(TAI).<br>
     * 
     * Since tai_utc() is driven by UTC the result will always assume the leap 
     * second has already occurred.  This function cannot explicitly show the leap 
     * seconds. 
     * @param mjd Modified Julian Date (UTC)
     * @return number of accumulated leap seconds
     */
    public static int tai_utc(double mjd_utc){
    	return leapsectbl.getAccumLeapSecs(mjd_utc);
    }
    
    /** Return the difference between TAI and UTC (known as accumulated leap 
     * seconds), where:<br>
     * TAI = UTC + accumLeapSecs.<br>
     * Returns zero if before the first leap second. Holds the value of the 
     * last specified leap second for dates after the last leap second.<br>
     * 
     * Note, TAI is a continuous time but UTC is not.  When crossing leap
     * seconds the same UTC time can result from different TAI times, as in:<br>
     * (at 0.25 sec steps, all times MJD)<br>
     * 42048.0001244213 TAI,   42047.999985532406 UTC, dATutc=12, dATtai=12<br>
     * 42048.000127314815 TAI, 42047.99998842592  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013020833 TAI,  42047.99999131944  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013310185 TAI,  42047.99999421296  UTC, dATutc=12, dATtai=12<br>
     * 42048.00013599537 TAI,  42047.999997106475 UTC, dATutc=12, dATtai=12<br>
     * (begin leap second)<br>
     * 42048.00013888889 TAI,  42048.0            UTC, dATutc=13, dATtai=12<br>
     * 42048.00014178241 TAI,  42048.00000289352  UTC, dATutc=13, dATtai=12<br>
     * 42048.00014467593 TAI,  42048.000005787035 UTC, dATutc=13, dATtai=12<br>
     * 42048.000147569444 TAI, 42048.00000868055  UTC, dATutc=13, dATtai=12<br>
     * (end leap second)<br>
     * 42048.00015046296 TAI,  42048.0            UTC, dATutc=13, dATtai=13<br>
     * 42048.00015335648 TAI,  42048.00000289352  UTC, dATutc=13, dATtai=13<br>
     * 42048.00015625 TAI,     42048.000005787035 UTC, dATutc=13, dATtai=13<br>
     * 42048.00015914352 TAI,  42048.00000868056  UTC, dATutc=13, dATtai=13<br>
     * 42048.00016203704 TAI,  42048.00001157408  UTC, dATutc=13, dATtai=13<br>
     * 42048.000164930556 TAI, 42048.000014467594 UTC, dATutc=13, dATtai=13<br>
     * where dATutc is the results of tai_utc(UTC) and <br>
     * dATtai = tai_utc_from_tai(TAI).<br>
     * 
     * This means that this function is different from tai_utc() in two ways:<br>
     * 1) The argument is TAI so this function has to check for leap second
     * crossings, and 2) this function result will accurately show the leap second,
     * unlike tai_utc(), as a duplicated UTC time.  During a leap second the 
     * results from tai_utc() will not match this function.
     * @param mjd Modified Julian Date (TAI)
     * @return number of accumulated leap seconds
     */
    public static int tai_utc_from_tai(double mjd_tai) {
    	/*
    	 * TAI and UTC aren't too different.  Use the accumulated leap seconds
    	 * based on TAI date as a first cut, then check to see if we're 
    	 * crossing a leap second boundary.
    	 */
    	// get the approximate leap second data based on TAI date: 
    	LeapSecondTable.LeapSecData leapsec = leapsectbl.getData(mjd_tai);
    	
    	// UTC, assuming not crossing a leap sec boundary:
    	double mjd_utc = mjd_tai - (leapsec.accumLeapSecs*sec2days);
    	
    	if (mjd_utc < leapsec.mjd_utc_start) {
    		// UTC and TAI straddle a leap second, get another leap second 
    		// based on our UTC 
    		leapsec = leapsectbl.getData(mjd_utc);
    	}
    	return leapsec.accumLeapSecs;
    }
	   
    /** Return the value of UTC time from the given GPS time in 
     * Modified Julian Date, corrected for leap seconds.  This function treats
     * GPS time as a linear timescale with a constant bias from TAI and 
     * assumes the GPS timescale extends linearly backwards before the true GPS 
     * epoch.
     * @param mjd_gps Modified Julian Date in GPS time
     * @return Modified Julian Date in UTC
     */
    public static double gps2utc(double mjd_gps){
    	double mjd_tai = mjd_gps + (TAI_GPS)*sec2days;
    	int leapsecs = tai_utc_from_tai(mjd_tai);
    	return mjd_tai - (leapsecs*sec2days);
    }
	   
	   /**
	    * Return Julian centuries since J2000 (Terrestrial Time)
	    * @param mjd_tt double containing TT in MJD format
	    * @return Julian centuries since J2000 (Terrestrial Time)
		*/
		public static double centuriesSinceJ2000(double mjd_tt){
			return (mjd_tt - TimeUtils.MJD_J2000)/36525.0;
		}
		
		/**
	     * Converts modified julian date to julian date.
	     * @param MJD Modified Julian Date
	     * @return JD Julian Date
	     */
	    public static double MJDtoJD(double MJD){
	        return MJD+2400000.5;
	    }
	    
	    /**
	     * Converts julian date to modified julian date.
	     * @param JD Julian Date
	     * @return MJD Modified Julian Date
	     */
	    public static double JDtoMJD(double JD){
	        return JD-2400000.5;
	    }
	    
	    /** Convert UTC time to TT.  Since this is driven by UTC inputs it 
	     * always assumes the leap second has occurred (e.g. produces a 
	     * discontinuity in the TT results).
	     * @param mjd_utc MJD of Current UTC time
	     * @return MJD of current TT.
	     */
	    public static double UTCtoTT(double mjd_utc){
	    	// TT = UTC + dAT + TT_TAI
	    	return mjd_utc + ((double)tai_utc(mjd_utc) + TT_TAI)*sec2days;
	    }
	    
	    /** Convert TT time to UTC.  This function produces duplicate UTC 
	     * times when crossing a leap second. 
	     * @param mjd_tt MJD in TT
	     * @return mjd_utc MJD in UTC
	     */
	    public static double TTtoUTC(double mjd_tt){	    	
	    	// get the number of leap seconds using TAI:
	    	int dAT = tai_utc_from_tai(mjd_tt - (TT_TAI*sec2days));
	    	// UTC = TT - dAT - TT_TAI
	        return mjd_tt - ((double)(dAT) + TT_TAI)*sec2days;
	    }	    
	    
	    /** Convert TDB time to TT
	     * Converted from NOVAS C code. Expressions used in this version are approximations resulting
	     * in accuracies of about 20 microseconds.
	     * 
	     * @param mjd_tbd
	     * @return mjd_tt
	     */	    
	    public static double TDBtoTT (double mjd_tdb) {
	    	double jd_tdb = MJDtoJD(mjd_tdb);
	    	double ecc = 0.01671022;
	    	double rev = 1296000.0;
	    	double T0 = 2451545.0;
	    	double RAD2SEC = 1.0/MathUtils.ARCSEC2RAD;
	    	double tdays = jd_tdb - T0;
	    	double m = ( 357.51716 + 0.985599987 * tdays ) * 3600.0;
	    	double l = ( 280.46435 + 0.985609100 * tdays ) * 3600.0;
	    	double lj = ( 34.40438 + 0.083086762 * tdays ) * 3600.0;
	    	m = MathUtils.mod(m,rev) / RAD2SEC;
	    	l = MathUtils.mod(l,rev) / RAD2SEC;
	    	lj = MathUtils.mod(lj,rev) / RAD2SEC;
	    	double e = m + ecc * Math.sin (m) + 0.5 * ecc * ecc * Math.sin (2.0 * m);
	    	double secdiff = 1.658e-3 * Math.sin (e) + 20.73e-6 * Math.sin (l - lj);
	    	// System.out.println("secdiff = "+secdiff);
	    	double jd_tt = jd_tdb - secdiff / 86400.0;
	    	double mjd_tt = JDtoMJD(jd_tt);
	    	return mjd_tt;
	    }
	    
	    /**
	     * Convert from the TT to TDB time system. 
	     * Accuracy is claimed to be 10 ns (see comments at end of file).
	     * @param MJD_TT Modified Julian Date of Terrestrial Time
	     * @return time in TDB expressed as a modified julian date
	     * 
	     */
	   public static double TTtoTDB (double MJD_TT){//long MJDint_in, double MJDfr_in) {
	     long MJDint = (long)Math.floor(MJD_TT);
	     double MJDfr = MJD_TT - MJDint;
	     //public TimeSys ts ; //* not used in this function

	     double tdbtdt =0;
	     double tdbtdtdot =0;
	     long oldmjd = 0 ;
	     long l ;

	     while ( MJDfr >= 1.0 ) {
	       MJDint++ ;
	       MJDfr-- ;
	     }
	     while ( MJDfr < 0.0 ) {
	       MJDint-- ;
	       MJDfr++ ;
	     }

	     if ( MJDint != oldmjd ) {
	       oldmjd = MJDint ;
	       l = oldmjd + 2400001 ;

	       tdbtdt = ctatv (l, 0.0) ;
	       tdbtdtdot = ctatv (l, 0.5) - ctatv (l, -0.5) ;
	     }
	     double TDB_minus_TT = ( tdbtdt + (MJDfr - 0.5) * tdbtdtdot );
	     return MJD_TT+TDB_minus_TT/86400.0;
	   }

	   //*** See bottom of this file for original header
	   /**
	    * Computes the cumulative relativistic time correction to
	    * earth-based clocks, TDB-TDT, for a given time. Routine
	    * furnished by the Bureau des Longitudes, modified by
	    * removal of terms much smaller than 0.1 microsecond.
	    * @param jdno. Julian day number of lookup
	    * @param fjdno. Fractional part of Julian day number
	    * @return Time difference TDB-TDT (seconds)
	    */
	   private static double ctatv(long jdno, double fjdno)
	   {

	     double t, tt, t1, t2, t3, t4, t5, t24, t25, t29, t30;
//	     double t31 ;

	         t = ((jdno-2451545) + fjdno)/(365250.0) ;
	         tt = t*t ;

	         t1  =       1656.674564 * Math.sin(  6283.075943033*t + 6.240054195)
	              +        22.417471 * Math.sin(  5753.384970095*t + 4.296977442)
	              +        13.839792 * Math.sin( 12566.151886066*t + 6.196904410)
	              +         4.770086 * Math.sin(   529.690965095*t + 0.444401603)
	              +         4.676740 * Math.sin(  6069.776754553*t + 4.021195093)
	              +         2.256707 * Math.sin(   213.299095438*t + 5.543113262)
	              +         1.694205 * Math.sin(    -3.523118349*t + 5.025132748)
	              +         1.554905 * Math.sin( 77713.772618729*t + 5.198467090)
	              +         1.276839 * Math.sin(  7860.419392439*t + 5.988822341)
	              +         1.193379 * Math.sin(  5223.693919802*t + 3.649823730)
	              +         1.115322 * Math.sin(  3930.209696220*t + 1.422745069)
	              +         0.794185 * Math.sin( 11506.769769794*t + 2.322313077)
	              +         0.600309 * Math.sin(  1577.343542448*t + 2.678271909)
	              +         0.496817 * Math.sin(  6208.294251424*t + 5.696701824)
	              +         0.486306 * Math.sin(  5884.926846583*t + 0.520007179)
	              +         0.468597 * Math.sin(  6244.942814354*t + 5.866398759)
	              +         0.447061 * Math.sin(    26.298319800*t + 3.615796498)
	              +         0.435206 * Math.sin(  -398.149003408*t + 4.349338347)
	              +         0.432392 * Math.sin(    74.781598567*t + 2.435898309)
	              +         0.375510 * Math.sin(  5507.553238667*t + 4.103476804) ;

	         t2  =          0.243085 * Math.sin(  -775.522611324*t + 3.651837925)
	              +         0.230685 * Math.sin(  5856.477659115*t + 4.773852582)
	              +         0.203747 * Math.sin( 12036.460734888*t + 4.333987818)
	              +         0.173435 * Math.sin( 18849.227549974*t + 6.153743485)
	              +         0.159080 * Math.sin( 10977.078804699*t + 1.890075226)
	              +         0.143935 * Math.sin(  -796.298006816*t + 5.957517795)
	              +         0.137927 * Math.sin( 11790.629088659*t + 1.135934669)
	              +         0.119979 * Math.sin(    38.133035638*t + 4.551585768)
	              +         0.118971 * Math.sin(  5486.777843175*t + 1.914547226)
	              +         0.116120 * Math.sin(  1059.381930189*t + 0.873504123)
	              +         0.101868 * Math.sin( -5573.142801634*t + 5.984503847)
	              +         0.098358 * Math.sin(  2544.314419883*t + 0.092793886)
	              +         0.080164 * Math.sin(   206.185548437*t + 2.095377709)
	              +         0.079645 * Math.sin(  4694.002954708*t + 2.949233637)
	              +         0.075019 * Math.sin(  2942.463423292*t + 4.980931759)
	              +         0.064397 * Math.sin(  5746.271337896*t + 1.280308748)
	              +         0.063814 * Math.sin(  5760.498431898*t + 4.167901731)
	              +         0.062617 * Math.sin(    20.775395492*t + 2.654394814)
	   	   +         0.058844 * Math.sin(   426.598190876*t + 4.839650148)
	              +         0.054139 * Math.sin( 17260.154654690*t + 3.411091093) ;

	         t3  =          0.048373 * Math.sin(   155.420399434*t + 2.251573730)
	              +         0.048042 * Math.sin(  2146.165416475*t + 1.495846011)
	              +         0.046551 * Math.sin(    -0.980321068*t + 0.921573539)
	              +         0.042732 * Math.sin(   632.783739313*t + 5.720622217)
	              +         0.042560 * Math.sin(161000.685737473*t + 1.270837679)
	              +         0.042411 * Math.sin(  6275.962302991*t + 2.869567043)
	              +         0.040759 * Math.sin( 12352.852604545*t + 3.981496998)
	              +         0.040480 * Math.sin( 15720.838784878*t + 2.546610123)
	              +         0.040184 * Math.sin(    -7.113547001*t + 3.565975565)
	              +         0.036955 * Math.sin(  3154.687084896*t + 5.071801441)
	              +         0.036564 * Math.sin(  5088.628839767*t + 3.324679049)
	              +         0.036507 * Math.sin(   801.820931124*t + 6.248866009)
	              +         0.034867 * Math.sin(   522.577418094*t + 5.210064075)
	              +         0.033529 * Math.sin(  9437.762934887*t + 2.404714239)
	              +         0.033477 * Math.sin(  6062.663207553*t + 4.144987272)
	              +         0.032438 * Math.sin(  6076.890301554*t + 0.749317412)
	              +         0.032423 * Math.sin(  8827.390269875*t + 5.541473556)
	              +         0.030215 * Math.sin(  7084.896781115*t + 3.389610345)
	              +         0.029862 * Math.sin( 12139.553509107*t + 1.770181024)
	              +         0.029247 * Math.sin(-71430.695617928*t + 4.183178762) ;

	         t4  =          0.028244 * Math.sin( -6286.598968340*t + 5.069663519)
	   	   +         0.027567 * Math.sin(  6279.552731642*t + 5.040846034)
	              +         0.025196 * Math.sin(  1748.016413067*t + 2.901883301)
	              +         0.024816 * Math.sin( -1194.447010225*t + 1.087136918)
	              +         0.022567 * Math.sin(  6133.512652857*t + 3.307984806)
	              +         0.022509 * Math.sin( 10447.387839604*t + 1.460726241)
	              +         0.021691 * Math.sin( 14143.495242431*t + 5.952658009)
	              +         0.020937 * Math.sin(  8429.241266467*t + 0.652303414)
	              +         0.020322 * Math.sin(   419.484643875*t + 3.735430632)
	              +         0.017673 * Math.sin(  6812.766815086*t + 3.186129845)
	              +         0.017806 * Math.sin(    73.297125859*t + 3.475975097)
	              +         0.016155 * Math.sin( 10213.285546211*t + 1.331103168)
	              +         0.015974 * Math.sin( -2352.866153772*t + 6.145309371)
	              +         0.015949 * Math.sin(  -220.412642439*t + 4.005298270)
	              +         0.015078 * Math.sin( 19651.048481098*t + 3.969480770)
	              +         0.014751 * Math.sin(  1349.867409659*t + 4.308933301)
	              +         0.014318 * Math.sin( 16730.463689596*t + 3.016058075)
	              +         0.014223 * Math.sin( 17789.845619785*t + 2.104551349)
	              +         0.013671 * Math.sin(  -536.804512095*t + 5.971672571)
	              +         0.012462 * Math.sin(   103.092774219*t + 1.737438797) ;

	         t5  =          0.012420 * Math.sin(  4690.479836359*t + 4.734090399)
	              +         0.011942 * Math.sin(  8031.092263058*t + 2.053414715)
	              +         0.011847 * Math.sin(  5643.178563677*t + 5.489005403)
	              +         0.011707 * Math.sin( -4705.732307544*t + 2.654125618)
	              +         0.011622 * Math.sin(  5120.601145584*t + 4.863931876)
	              +         0.010962 * Math.sin(     3.590428652*t + 2.196567739)
	              +         0.010825 * Math.sin(   553.569402842*t + 0.842715011)
	              +         0.010396 * Math.sin(   951.718406251*t + 5.717799605)
	              +         0.010453 * Math.sin(  5863.591206116*t + 1.913704550)
	              +         0.010099 * Math.sin(   283.859318865*t + 1.942176992)
	              +         0.009858 * Math.sin(  6309.374169791*t + 1.061816410)
	              +         0.009963 * Math.sin(   149.563197135*t + 4.870690598)
	              +         0.009370 * Math.sin(149854.400135205*t + 0.673880395) ;

	         t24 = t * (  102.156724 * Math.sin(  6283.075849991*t + 4.249032005)
	              +         1.706807 * Math.sin( 12566.151699983*t + 4.205904248)
	              +         0.269668 * Math.sin(   213.299095438*t + 3.400290479)
	              +         0.265919 * Math.sin(   529.690965095*t + 5.836047367)
	              +         0.210568 * Math.sin(    -3.523118349*t + 6.262738348)
	              +         0.077996 * Math.sin(  5223.693919802*t + 4.670344204) ) ;

	         t25 = t * (    0.059146 * Math.sin(    26.298319800*t + 1.083044735)
	              +         0.054764 * Math.sin(  1577.343542448*t + 4.534800170)
	              +         0.034420 * Math.sin(  -398.149003408*t + 5.980077351)
	              +         0.033595 * Math.sin(  5507.553238667*t + 5.980162321)
	              +         0.032088 * Math.sin( 18849.227549974*t + 4.162913471)
	              +         0.029198 * Math.sin(  5856.477659115*t + 0.623811863)
	              +         0.027764 * Math.sin(   155.420399434*t + 3.745318113)
	              +         0.025190 * Math.sin(  5746.271337896*t + 2.980330535)
	              +         0.024976 * Math.sin(  5760.498431898*t + 2.467913690)
	   	   +         0.022997 * Math.sin(  -796.298006816*t + 1.174411803)
	              +         0.021774 * Math.sin(   206.185548437*t + 3.854787540)
	              +         0.017925 * Math.sin(  -775.522611324*t + 1.092065955)
	              +         0.013794 * Math.sin(   426.598190876*t + 2.699831988)
	              +         0.013276 * Math.sin(  6062.663207553*t + 5.845801920)
	              +         0.012869 * Math.sin(  6076.890301554*t + 5.333425680)
	              +         0.012152 * Math.sin(  1059.381930189*t + 6.222874454)
	              +         0.011774 * Math.sin( 12036.460734888*t + 2.292832062)
	              +         0.011081 * Math.sin(    -7.113547001*t + 5.154724984)
	              +         0.010143 * Math.sin(  4694.002954708*t + 4.044013795)
	              +         0.010084 * Math.sin(   522.577418094*t + 0.749320262)
	              +         0.009357 * Math.sin(  5486.777843175*t + 3.416081409) ) ;

	         t29 = tt * (   0.370115 * Math.sin(                     4.712388980)
	              +         4.322990 * Math.sin(  6283.075849991*t + 2.642893748)
	              +         0.122605 * Math.sin( 12566.151699983*t + 2.438140634)
	              +         0.019476 * Math.sin(   213.299095438*t + 1.642186981)
	              +         0.016916 * Math.sin(   529.690965095*t + 4.510959344)
	              +         0.013374 * Math.sin(    -3.523118349*t + 1.502210314) ) ;

	         t30 = t * tt * 0.143388 * Math.sin( 6283.075849991*t + 1.131453581) ;

	         return (t1+t2+t3+t4+t5+t24+t25+t29+t30) * 1.0e-6 ;
	   }

}
//*** Original .h file header
/* RCS: $Id: bary.h,v 3.4 2001/06/28 21:47:02 arots Exp arots $ */
/*-----------------------------------------------------------------------
 *
 *  bary.h
 *
 *  Date: 19 November 1997
 *  Unified FITS version
 *
 *  Arnold Rots, USRA
 *
 *  This is the header file for bary.c and related source files.
 *
 *  Do not ever use DE405 on a FITS file that has:
 *      TIMESYS='TDB'
 *    but not:
 *      RADECSYS='ICRS' or PLEPHEM='JPL-DE405'
 *  TIMESYS='TDB' _with_    RADECSYS='ICRS' and/or PLEPHEM='JPL-DE405'
 *                             should use denum=405.
 *  TIMESYS='TDB' _without_ RADECSYS='ICRS' or     PLEPHEM='JPL-DE405'
 *                             should use denum=200.
 *
 *  Time is kept in three possible ways:
 *  The basic convention is:
 *    MJDTime t ;  The MJDTime struct is defined below;
 *                   a C++ class would be better.
 *  The JPL ephemeris functions use:
 *    double time[2] ;  Where:
 *      time[0]   Integer part of JD
 *      time[1]   Fractional part of JD; -0.5 <= time[1] < 0.5
 *  MET uses a single double:
 *    double time ;  Where:
 *      time      Seconds or days since MJDREF
 *
 * We shall use the following convention:
 *
 *   If the FITS file has   Then we shall adopt
 *   TIMESYS  TIMEUNIT      refTS  fromTS  mjdRef
 *
 *     TT        s           TT     TT     MJDREF
 *     TT        d           TT     TT     MJDREF
 *     UTC       s           UTC    TT     toTT(MJDREF)
 *     UTC       d           UTC    UTC    MJDREF
 *     TDB       s           TDB    TDB    MJDREF
 *     TDB       d           TDB    TDB    MJDREF
 *
 * System files (JPLEPH, JPLEPH.200, JPLEPH.405, psrtime.dat, psrbin.dat,
 * tai-utc.dat, tdc.dat) are first searched for in $TIMING_DIR,
 * then $LHEA_DATA ($ASC_DATA, depending on what LHEADIR is set to).
 *
 *----------------------------------------------------------------------*/

//*** Original .c file header
/*-----------------------------------------------------------------------
*
*  bary.c
*
*  Date: 20 November 1997
*
*  Arnold Rots, USRA
*
*  bary contains a set of functions around barycorr, calculating barycenter
*  corrections for, in principle, any observations.
*  Required:
*    bary.h
*    dpleph.c
*    cfitsio
*
*  Externally referenced:
*    int baryinit (enum Observatory, char *, char *, double, double, char *, int) ;
*    double barycorr (MJDTime *, double *, double *) ;
*    int timeparms (char *, MJDTime *, double, int, PsrTimPar *, PsrBinPar *, int) ;
*    double absphase (MJDTime *,  PsrTimPar *, int, PsrBinPar *,
*                     int, double *, char *, int *) ;
*    double xabsphase (double, PsrTimPar *, int, PsrBinPar *,
*                     int, double *, char *, int *) ;
*    fitsfile *openFFile (const char *name)
*    double *scorbit (char *, MJDTime *, int *oerror)
*    void met2mjd (double, MJDTime *) ;
*    double mjd2met (MJDTime *) ;
*    const char *convertTime (MJDTime *) ;
*    const char *fitsdate ()
*    double ctatv (long, double) ;
*    void c200to405 (double *, int) ;
*
*  MET equivalents:
*    double xbarycorr (double, double *, double *) ;
*    double *xscorbit (char *, double, int *) ;
*
*  Internal:
*    double TTtoTDB (MJDTime *) ;
*    double TTtoTCG (MJDTime *) ;
*    double TDBtoTCB (MJDTime *) ;
*    double UTCtoTAI (MJDTime *) ;
*    double toTT (MJDTime *) ;
*    double toUTC (MJDTime *) ;
*    double binorbit (MJDTime *, PsrBinPar *) ;
*
*  It is assumed that the environment variable $TIMING_DIR is defined
*  and that the ephemeris file $TIMING_DIR/ephem.fits exists.
*  It is also assumed that $TIMING_DIR/psrtime.dat and
*  $TIMING_DIR/psrbin.dat exist.
*  However, all these names are defined as mocros in bary.h.
*
*  Certain parts are adapted from software
*  by Joseph M. Fierro, Stanford University, EGRET project,
*  and from software provided with the JPL ephemeris.
*
*----------------------------------------------------------------------*/

//*** original TTtoTDB function header
/*-----------------------------------------------------------------------
*
*  TTtoTDB calculates TDB-TT at time TT.
*  double TTtoTDB (long jdTTint, double jdTTfr)
*   mjdTT:   MJD (TT) day
*   return:  TDB - TT (s)
*
*  It uses the coefficients from Fairhead & Bretagnon 1990,
*  A&A 229, 240, as provided by ctatv.
*  The accuracy is better than 100 ns.
*
*  The proper way to do all this is to abandon TDB and use TCB.
*
*  The way this is done is as follows: TDB-TT and its derivative are
*  calculated for the integer part of the Julian Day (when needed).
*  The precise value is derived from fractional part and derivative.
*
*----------------------------------------------------------------------*/

//*** original ctatv function header
/*  $Id: ctatv.c,v 3.0 1998/08/25 19:26:49 arots Exp arots $
 *----------------------------------------------------------------------
 *
 *     Routine name: CTATV
 *
 *     Programmer and Completion Date:
 *        Lloyd Rawley - S.T.X. - 07/28/89
 *        (Retyped by Masaki Mori - 04/19/96)
 *        (Converted to C for bary, optimized,
 *         corrected following Faithead & Bretagnon 1990, A&A 229, 240
 *         by Arnold Rots - 1997-11-20)
 *
 *     Function: Computes the cumulative relativistic time correction to
 *               earth-based clocks, TDB-TDT, for a given time. Routine
 *               furnished by the Bureau des Longitudes, modified by
 *               removal of terms much smaller than 0.1 microsecond.
 *
 *     Calling Sequence:  tdbdt = ctacv(jdno, fjdno)
 *        Argument   Type   I/O                Description
 *        --------   ----   ---  ----------------------------------------
 *        jdno      long     I   Julian day number of lookup
 *        fjdno     double   I   Fractional part of Julian day number
 *        tdbtdt    double   O   Time difference TDB-TDT (seconds)
 *
 *     Called by:  TTtoTDB
 *
 *     Calls:  none
 *
 *     COMMON use:  none
 *
 *     Significant Local Variables:
 *        Variable   Type   Ini. Val.        Description
 *        --------   ----   ---------  ----------------------------------
 *          t       double      -      Time since 2000.0 (millennia)
 *          tt      double      -      Square of T
 *
 *     Logical Units Used:  none
 *
 *     Method:
 *        Convert input time to millennia since 2000.0
 *        For each sinusoidal term of sufficient amplitude
 *           Compute value of term at input time
 *        End for
 *        Add together all terms and convert from microseconds to seconds
 *     End CTATV
 *
 *     Note for the retyped version:
 *        Results of this routine has been confirmed up to (1E-10)microsecond
 *        level compared with the calculation by the IBM machine: this seems
 *        to be within the 64-bit precision allowance, but still the original
 *        hardcopy should be kept as the right one. (M.M.)
 *
 *     Note for the C version: the accuracy is guaranteed to 100 ns.
 *
 *---------------------------------------------------------------------------*/

