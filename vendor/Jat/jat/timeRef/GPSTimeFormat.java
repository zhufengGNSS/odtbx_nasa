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
import jat.math.*;
import jat.spacetime.*;

/**
 * <P>
 * The GPSTimeFormat Class provides methods of converting between the GPS Time and MJD and CalDate time formats.
 * Reference: Satellite Orbits by Montenbruck and Gill.
 * Translated code from the NOAA GPS Toolbox (www.ngs.noaa.gov/gps-toolbox/index.html)
 *
 * @deprecated
 * @see jat.spacetime.GPSTimeFormat
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GPSTimeFormat {

    private long GPS_Week;

    private double GPS_SOW;

    /** MJD of GPS epoch (Jan 6, 1980).
     */
    public static final long JAN61980 = 44244;

    /** MJD of Jan 1, 1901
     */
    public static final long JAN11901 = 15385;

    /** Seconds per day.
     */
    public static final double SEC_PER_DAY = 86400.0;

    public static final double TAI_TT = 32.184;  // constant

    public static final int TAI_GPS = 19;  // constant

    /** Create a GPSTimeFormat object using GPS Week and Seconds of the Week.
     * @param gps_week Long GPS Week number.
     * @param gps_sow Double GPS seconds of the week.
     */
    public GPSTimeFormat(long gps_week, double gps_sow){
        this.GPS_Week = gps_week;
        this.GPS_SOW = gps_sow;
    }

    /** Create a clone of a GPSTimeFormat object.
     * @param gps GPSTimeFormat object to be cloned.
     */
    public GPSTimeFormat(GPSTimeFormat gps){
        this.GPS_Week = gps.GPS_Week;
        this.GPS_SOW = gps.GPS_SOW;
    }

    /** Creates a new instance of GPSTimeFormat
     * @param year Integer year.
     * @param doy Integer day of year.
     * @param hour Integer hours.
     * @param minute Integer minutes.
     * @param second Integer seconds.
     */
    public GPSTimeFormat(int year, int doy, int hour, int minute, double second) {
        long mjd;
        double fmjd;

        mjd = ((year - 1901)/4)*1461 + ((year - 1901)%4)*365 + doy - 1 + JAN11901;
        fmjd = ((second/60.0 + minute)/60.0 + hour)/24.0;

        this.GPS_Week = (mjd - JAN61980)/7;
        this.GPS_SOW = ( (mjd - JAN61980) - GPS_Week*7 + fmjd )*SEC_PER_DAY;
    }

    /** Create a GPSTimeFormat object from a CalDate object.
     * @param cal CalDate object.
     */
    public GPSTimeFormat(CalDate cal){
        long mjd;
        double fmjd;
        int year = cal.Year;
        int doy = cal.DOY;
        int hour = cal.Hour;
        int minute = cal.Min;
        double second = cal.Sec;

        mjd = ((year - 1901)/4)*1461 + ((year - 1901)%4)*365 + doy - 1 + JAN11901;
        fmjd = ((second/60.0 + minute)/60.0 + hour)/24.0;

        this.GPS_Week = (mjd - JAN61980)/7;
        this.GPS_SOW = ( (mjd - JAN61980) - GPS_Week*7 + fmjd )*SEC_PER_DAY;
    }

    /** Create a GPSTimeFormat object from an MJD.
     * @param mjd Modified Julian Date.
     */
    public GPSTimeFormat(double mjd){
        long lmjd = (long)mjd;
        double fmjd = MathUtils.Frac(mjd);
        this.GPS_Week = (lmjd - JAN61980)/7;
        this.GPS_SOW = ( (lmjd - JAN61980) - GPS_Week*7 + fmjd )*SEC_PER_DAY;
    }

    /** Returns the modified Julian date.
     * @return modified Julian date.
     */
    public double mjd(){
        double out = (double)(GPS_Week*7 + GPS_SOW/SEC_PER_DAY + JAN61980);
        return out;
    }

    /** Returns the time in CalDate format.
     * @return Time in CalDate format.
     */
    public CalDate calDate(){
        CalDate out = new CalDate(this);
        return out;
    }

    /** Returns the GPS Week number.
     * @return GPS week number.
     */
    public long gps_week(){
        return this.GPS_Week;
    }

    /** Returns the GPS Seconds of the Week.
     * @return GPS seconds of the week.
     */
    public double gps_sow(){
        return this.GPS_SOW;
    }

    /**
     * Print out the GPS time
     * @param title String containing a title
     */
    public void print(String title) {
        System.out.println(title + " GPS Time");
        System.out.println("   GPS Week: "+this.GPS_Week);
        System.out.println("   GPS SOW: "+this.GPS_SOW);
        System.out.println("   MJD: "+this.mjd());
    }

	/**
	 * Print out the GPS time in CalDate format
	 * @param title String containing a title
	 */
    public void printCalDate(String title){
        CalDate temp = this.calDate();
        temp.print(title);
    }

    /** Increments time forward or backwards.
     * @param sec Increment to move time in seconds. Positive for forward, negative for backwards.
     */
    public void increment(double sec){
        double temp = this.GPS_SOW + sec;
        double max_sow = 604800.0;
        if (Math.abs(sec) > max_sow){
            System.out.println("Sorry I can't increment time by >= 1 week");
            return;
        }
        // incrementing backwards into last week
        if (temp < 0.0){
            this.GPS_Week = this.GPS_Week - 1;
            this.GPS_SOW = max_sow + temp;
        }
        // incrementing forward into next week
        if (temp >= max_sow){
            this.GPS_Week = this.GPS_Week + 1;
            this.GPS_SOW = temp - max_sow;
        }
        // incrementing within the week
        if ((temp >= 0.0)&&(temp < max_sow)){
            this.GPS_SOW = temp;
        }
    }

    /** Convert GPS time to UTC time.
     * @return CalDate object with current UTC time.
     */
    public CalDate GPS2UTC(){
        // get the current mjd
        double mjd = this.mjd();

        // convert the current time to CalDate format
        CalDate out = new CalDate(this);

        // compute the difference between GPS and UTC
        int utc_gps = TAI_GPS - TimeUtils.tai_utc(mjd);

        // convert from GPS time to UTC time
        out.increment(utc_gps);
        return out;
    }

    /** Convert GPS time to TT.
     * @return mjd with current TT.
     */
    public double GPS2TT(){
        // get the current mjd
        double mjd = this.mjd();

        // compute the difference between Terrestrial and GPS
        double tt_gps = TAI_GPS - TAI_TT;

        // make a copy
        GPSTimeFormat out = new GPSTimeFormat(this);

        // convert from GPS time to Terrestrial time
        out.increment(tt_gps);
        return out.mjd();
    }

}
