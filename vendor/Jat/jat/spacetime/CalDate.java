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
import jat.math.*;
import java.io.*;

/**
 * <P>
 * The CalDate Class provides methods of converting between the CalDate and MJD and GPSTimeFormat formats.
 * Reference: Satellite Orbits by Montenbruck and Gill.
 *
 * Note: CalDate is just a time format. The actual time depicted could be GPS time, UTC, TAI, etc.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class CalDate implements Serializable {

    /** Year.
     */
    protected int Year;
    /** Month.
     */
    protected int Month;
    /** Day.
     */
    protected int Day;
    /** Hours.
     */
    protected int Hour;
    /** Minutes.
     */
    protected int Min;
    /** Seconds.
     */
    protected double Sec;

    /** Day of year.
     */
    protected int DOY;

    /** Construct a time object using a calendar date.
     * @param Yr Year.
     * @param Mon Month.
     * @param D Day.
     * @param Hr Hours.
     * @param Mn Minutes.
     * @param S Seconds.
     */
    public CalDate( int Yr, int Mon, int D, int Hr, int Mn, double S ) {
        // accept the inputs, probably should add some checking some day
        this.Year = Yr;
        this.Month = Mon;
        this.Day = D;
        this.Hour = Hr;
        this.Min = Mn;
        this.Sec = S;
        this.DOY = TimeUtils.day2doy(Yr, Mon, D);
    }

    /** Construct a clone of a CalDate object.
     * @param cal CalDate object to be cloned.
     */
    public CalDate( CalDate cal) {
        this.Year = cal.Year;
        this.Month = cal.Month;
        this.Day = cal.Day;
        this.Hour = cal.Hour;
        this.Min = cal.Min;
        this.Sec = cal.Sec;
        this.DOY = cal.DOY;
    }

    /** Create a CalDate object using MJD. From Montenbruck C++ code.
     * @param mjd Modified Julian Date.
     */
    public CalDate(double mjd){
        long    a,b,c,d,e,f;
        double  Hours,x;

        // Convert Julian day number to calendar date
        a = (long)(mjd+2400001.0);

        if ( a < 2299161 ) {  // Julian calendar
            b = 0;
            c = a + 1524;
        }
        else {                // Gregorian calendar
            b = (long)((a-1867216.25)/36524.25);
            c = a +  b - (b/4) + 1525;
        }

        d     = (long)( (c-122.1)/365.25 );
        e     = 365*d + d/4;
        f     = (long)( (c-e)/30.6001 );

        long temp = (long)(30.6001*f);
        this.Day   = (int)(c - e - temp);
        temp = (long)(f/14);
        this.Month = (int)(f - 1 - 12*temp);
        temp = (long)((7+Month)/10);
        this.Year  = (int)(d - 4715 - temp);

        Hours = 24.0*(mjd-Math.floor(mjd));

        this.Hour = (int)Hours;
        x = (Hours-Hour)*60.0;
        this.Min = (int) x;
        this.Sec = (x-Min)*60.0;
        this.DOY = TimeUtils.day2doy(this.Year, this.Month, this.Day);
    }

    /** Create a CalDate object using GPS Time. Translated code from the
     * NOAA GPS Toolbox (www.ngs.noaa.gov/gps-toolbox/index.html).
     * @param gps GPSTimeFormat object.
     */
    public CalDate(GPSTimeFormat gps){
        long GPS_Week = gps.gps_week();
        double GPS_SOW = gps.gps_sow();

        long mjd, days_fr_jan1_1901;
        double fmjd;
        long delta_yrs, num_four_yrs, years_so_far, days_left;

        mjd = (long)(GPS_Week*7 + GPS_SOW/GPSTimeFormat.SEC_PER_DAY + GPSTimeFormat.JAN61980);
        fmjd = MathUtils.mod(GPS_SOW, GPSTimeFormat.SEC_PER_DAY)/GPSTimeFormat.SEC_PER_DAY;

        days_fr_jan1_1901 = mjd - GPSTimeFormat.JAN11901;
        num_four_yrs = days_fr_jan1_1901/1461;
        years_so_far = 1901 + 4*num_four_yrs;
        days_left = days_fr_jan1_1901 - 1461*num_four_yrs;
        delta_yrs = days_left/365 - days_left/1460;

        this.Year = (int)(years_so_far + delta_yrs);
        this.DOY = (int)(days_left - 365*delta_yrs + 1);
        this.Month = TimeUtils.doy2month(this.Year, this.DOY);
        this.Day = TimeUtils.doy2day(this.Year, this.DOY);
        this.Hour = (int)(fmjd*24.0);
        this.Min = (int)(fmjd*1440.0 - Hour*60.0);
        this.Sec = fmjd*86400.0 - Hour*3600.0 - Min*60.0;
    }

    /** Returns the modified Julian date. From Montenbruck C++ code.
     * @return modified Julian date.
     */
    public double mjd(){
        long    MjdMidnight;
        double  FracOfDay;
        int     b;
        int year = this.Year;
        int month = this.Month;
        int day = this.Day;
        int hour = this.Hour;
        int min = this.Min;
        double sec = this.Sec;

        // compute the MJD
        if (month<=2) {
            month += 12;
            --year;
        }

        if ( (10000L*year+100L*month+day) <= 15821004L ){
            b = -2 + ((year+4716)/4) - 1179;     // Julian calendar
        }
        else {
            b = (year/400)-(year/100)+(year/4);  // Gregorian calendar
        }

        int temp = (int)(30.6001*(month+1));

        MjdMidnight = 365L*year - 679004L + b + temp + day;
        FracOfDay   = (hour+min/60.0+sec/3600.0) / 24.0;

        double mjd = MjdMidnight + FracOfDay;
        return mjd;
    }

    /** Return the time in GPSTimeFormat format. From GPS Toolkit code.
     * @return Time in GPSTimeFormat format.
     */
    public GPSTimeFormat gpsTime(){
        GPSTimeFormat out = new GPSTimeFormat(this);
        return out;
    }


    /** Increments time forward or backwards.
     * @param sec Increment to move time in seconds. Positive for forward, negative for backwards.
     */
    public void increment(double sec){
        GPSTimeFormat gps = new GPSTimeFormat(this);
        gps.increment(sec);
        CalDate out = new CalDate(gps);
        this.Year = out.Year;
        this.Month = out.Month;
        this.Day = out.Day;
        this.Hour = out.Hour;
        this.Min = out.Min;
        this.Sec = out.Sec;
        this.DOY = out.DOY;
    }

    /** Returns the year.
     * @return Year.
     */
    public int year(){
        return this.Year;
    }

    /** Returns the month.
     * @return Month.
     */
    public int month(){
        return this.Month;
    }

    /** Return the day of the month.
     * @return Day of month.
     */
    public int day(){
        return this.Day;
    }

    /** Return the hour.
     * @return Hour.
     */
    public int hour(){
        return this.Hour;
    }

    /** Return the minutes.
     * @return Minute.
     */
    public int min(){
        return this.Min;
    }

    /** Return the second.
     * @return Second.
     */
    public double sec(){
        return this.Sec;
    }

    /** Return the day of year.
     * @return Day of year.
     */
    public int doy(){
        return this.DOY;
    }
    /** Return the seconds of the day.
     * @return seconds of the day.
     */
    public double sec_of_day(){
        return (this.Hour*3600.0 + this.Min*60.0 + this.Sec);
    }


    /** Convert UTC time to GPS time.
     * @return CalDate object with current UTC time.
     */
    public GPSTimeFormat UTC2GPS(){
        // get the current mjd
        double mjd = this.mjd();

        // compute the difference between GPS and UTC
        int gps_utc = TimeUtils.tai_utc(mjd) - TimeUtils.TAI_GPS;

        // convert the current time to GPSTimeFormat format
        GPSTimeFormat out = new GPSTimeFormat(this);

        // convert from UTC time to GPS time
        out.increment(gps_utc);
        return out;
    }

    /** Convert UTC time to TT.
     * @param mjd_utc MJD of Current UTC time
     * @return MJD of current TT.
     */
    public static double UTC2TT(double mjd_utc){

        // compute the difference between TT and UTC
        double tt_utc = (double)(TimeUtils.tai_utc(mjd_utc) + TimeUtils.TT_TAI);
        double out = mjd_utc + tt_utc/86400.0;
        return out;
    }

    /**
     * Send it to string
     * @return String containing the calendar date
     */
    public String toString(){
        String out = this.Year+","+this.Month+","+this.Day+" "+this.Hour+":"+this.Min+":"+this.Sec;
        return out;
    }

    /**
     * Print out the date
     * @param title String containing a title
     */
    public void print(String title){
        System.out.println(title);
        System.out.println(this.Year+","+this.Month+","+this.Day+" "+this.Hour+":"+this.Min+":"+this.Sec);
    }
}
