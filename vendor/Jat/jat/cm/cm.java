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

package jat.cm;

import jat.math.*;
import java.util.*;
import jat.matvec.data.*;
import java.io.*;

/**
 * Methods for Celestial Mechanics
 * 
 * @author Tobias Berthold
 *  Date        :   3-29-2002
 * @version 1.0
 */
public class cm
{
	Constants c = new Constants();

	// Constants
	double AU = 149597870.;
	double TU = 58.132440906;
	public static double sun_radius = 1391980.;
	public static double mercury_radius = 4880.;
	public static double mercury_a = 57900000.; // mean semimajor axis of mercury orbit around sun
	public static double venus_radius = 6051.8;
	public static double venus_a = 108210000.; // mean semimajor axis of mercury orbit around sun
	public static double earth_radius = 6378.;
	public static double earth_a = 149597892.; // mean semimajor axis of earth orbit around sun
	public static double mars_radius = 3397.;
	public static double mars_a = 227939186.; // mean semimajor axis of Mars orbit around sun
	public static double jupiter_radius = 71492;
	public static double jupiter_a = 778412010;
	public static double moon_radius = 1738.;
	public static double moon_a = 384746.924;
	public static double moon_obl = 23.3439; // obliquity of the ecliptic in degrees
	public static double moon_period = 2360591.47; // in seconds

	// Variables
	public static double mu = 3.986e5;
	//static double mu=1.;

	public cm(double mu)
	{
		cm.mu = mu;
	}

	/**
	 * @author Tobias Berthold
	 * Simple date class. The Java Calendar class is typically used for storing dates, 
	 * but it has a large overhead for some applications.  
	 *
	 */
	public static class Date
	{
		public long year, month, day, hour, minute, second, timezone;
	}

	/** Orbital elements for the planets from "Explanatory supplement to the Astronomical 
	 * Almanac" by Kenneth Seidelmann 
	 */
	public static final KeplerElements mars_elements =
		new KeplerElements(227939186, 0.0934006199474, Rad(1.85061), Rad(49.57854), Rad(336.04084), 355.45332);
	public static final KeplerElements earth_moon_elements =
		new KeplerElements(149598023, 0.01671022, Rad(0.00005), Rad(-11.26064), Rad(102.94719), 100.46435);

	/**
	 * Convert Calendar date format to Julian date
	 * @param c Date in Calendar format as in java.util.Calendar 
	 * @return Julian date
	 */
	public static double juliandate(Calendar c)
	{
		Date d = new Date();
		d.year = c.get(Calendar.YEAR);
		d.month = c.get(Calendar.MONTH) + 1;
		d.day = c.get(Calendar.DAY_OF_MONTH);
		d.hour = c.get(Calendar.HOUR_OF_DAY);
		d.minute = c.get(Calendar.MINUTE);
		d.second = c.get(Calendar.SECOND);
		return juliandate(d);
	}

	/**
	 * Convert the Julian day number to date in internal format
	 * 
	 * @param JD
	 */
	public static Date Date(double JD)
	{
		long I, A, B, C, D, E, G;
		double decimalday, decimalhour, decimalminute;
		Date d = new Date();

		double JD1 = JD + 0.5; // Add 0.5 to JD
		I = (long)Math.floor(JD1); // set I to integer part
		double F = MathUtils.Frac(JD1); // set F to fractional part
		if (I > 2299160)
		{
			A = (long)Math.floor((I - 1867216.25) / 36524.25); // set A to ..
			B = I + 1 + A - (long)Math.floor(A / 4);
		} else
			B = I;
		C = B + 1524;
		D = (long)Math.floor((C - 122.1) / 365.25);
		E = (long)Math.floor(365.25 * D);
		G = (long)Math.floor((C - E) / 30.6001);
		decimalday = (C - E - (long)Math.floor(30.6001 * G)) + F;
		d.day = (long)Math.floor(decimalday);
		if (G < 14)
			d.month = G - 1;
		else
			d.month = G - 13;
		if (d.month < 3)
			d.year = D - 4715;
		else
			d.year = D - 4716;

		decimalhour = 24 * F;
		d.hour = (long)decimalhour;
		decimalminute = 60 * (decimalhour - d.hour);
		d.minute = (long)decimalminute;
		d.second = (long) (60 * (decimalminute - d.minute));

		return d;
	}

	/**
	 * Convert Julian date to Calendar date format
	 * Algortihm from book "Practical Astronomy with your Calculator" by Peter Duffett-Smith 
	 * @param JD Julian date
	 * @return Date in Calendar format as in java.util.Calendar
	 */
	public static Calendar Calendar(double JD)
	{
		long I, A, B, C, D, E, G;
		double decimalday, decimalhour, decimalminute;
		Date d = new Date();

		double JD1 = JD + 0.5; // Add 0.5 to JD
		I = (long)Math.floor(JD1); // set I to integer part
		double F = MathUtils.Frac(JD1); // set F to fractional part
		if (I > 2299160)
		{
			A = (long)Math.floor((I - 1867216.25) / 36524.25); // set A to ..
			B = I + 1 + A - (long)Math.floor(A / 4);
		} else
			B = I;
		C = B + 1524;
		D = (long)Math.floor((C - 122.1) / 365.25);
		E = (long)Math.floor(365.25 * D);
		G = (long)Math.floor((C - E) / 30.6001);
		decimalday = (C - E - (long)Math.floor(30.6001 * G)) + F;
		d.day = (long)Math.floor(decimalday);
		if (G < 14)
			d.month = G - 1;
		else
			d.month = G - 13;
		if (d.month < 3)
			d.year = D - 4715;
		else
			d.year = D - 4716;

		decimalhour = 24 * F;
		d.hour = (long)decimalhour;
		decimalminute = 60 * (decimalhour - d.hour);
		d.minute = (long)decimalminute;
		d.second = (long) (60 * (decimalminute - d.minute));

		Calendar cal =
			new GregorianCalendar((int)d.year, (int)d.month - 1, (int)d.day, (int)d.hour, (int)d.minute, (int)d.second);
		return cal;
	}

	public static double juliandate(Date d)
	{
		return juliandate(d.year, d.month, d.day, d.hour, d.minute, d.second);
	}

	public static double juliandate(long year, long month, long day, long hour, long minute, long second)
	{
		/* after Oct 15th, 1582 */
		long j_year = year;
		long j_month = month;
		long A, B, C, D;

		if (month == 1 || month == 2)
		{
			j_month = month + 12;
			j_year = year - 1;
		}
		A = (long) (j_year / 100);
		B = 2 - A + (long) (A / 4);
		C = (long) (365.25 * j_year);
		D = (long) (30.6001 * (j_month + 1));
		return B + C + D + decimal_day(day, hour, minute, second) + 1720994.5;
	}

	/*
	horizon_coordinates equatorial_to_horizon(double latitude, double declination, double H)
	{
	    double phi,del,Ha;
	    double altemp1,altemp2,aztemp1,aztemp2,aztemp3;
	    phi=deg_to_rad(latitude);
	    del=deg_to_rad(declination);
	    Ha=deg_to_rad(H*15);
	    
	    
	    hc=new horizon_coordinates();
	    //hc.azimuth=java.lang.Math.sin(declination);
	    //hc.azimuth=Math.sin(del)*Math.sin(phi)+Math.cos(del)*Math.cos(phi)*Math.cos(Ha);
	    //hc.altitude=rad_to_deg(Math.asin(hc.azimuth));
	    altemp1=Math.sin(del)*Math.sin(phi)+Math.cos(del)*Math.cos(phi)*Math.cos(Ha);
	    altemp2=Math.asin(altemp1);
	    hc.altitude=rad_to_deg(altemp2);
	    aztemp1=(Math.sin(del)-Math.sin(phi)*altemp1)/(Math.cos(phi)*Math.cos(altemp2));
	    aztemp2=rad_to_deg(Math.acos(aztemp1));
	    if(Math.sin(Ha)<0)
	        aztemp3=aztemp2;
	    else
	        aztemp3=360-aztemp2;
	    hc.azimuth=aztemp3;
	
	    return hc;
	}
	*/
	//horizon_coordinates hc=new horizon_coordinates();

	public static double hour_angle(double LST, double RA)
	{
		return get_range(LST - RA);
	}

	public static double GST_decimal(
		long year,
		long month,
		long day,
		long hour,
		long minute,
		long second,
		long timezone)
	{
		double S, T, T0, UT1, GST;

		S = juliandate(year, month, day, 0, 0, 0) - 2451545.0;
		T = S / 36525;
		T0 = get_range(6.697374558 + 2400.051336 * T + 0.000025862 * T * T);
		UT1 = UT_decimal(hour, minute, second, timezone) * 1.002737909;
		GST = get_range(T0 + UT1);
		return GST;
	}

	public static double LST_decimal(
		long year,
		long month,
		long day,
		long hour,
		long minute,
		long second,
		double latitude,
		double longitude,
		long timezone)
	{

		double LST;

		LST = GST_decimal(year, month, day, hour, minute, second, timezone) + longitude / 15;

		return LST;
	}

	public static double decimal_day(long day, long hour, long minute, long second)
	{
		double temp = day + decimal_hour(hour, minute, second) / 24;
		return temp;
	}

	public static double get_range(double number)
	{
		double factor = java.lang.Math.abs((long)number / 24);
		if (number > 24)
			number = number - 24 * factor;
		if (number < 0)
			number = number + 24 * (factor + 1);
		return number;
	}

	public static double decimal_hour(long hour, long minute, long second)
	{
		double temp = (double)hour + (double)minute / 60 + (double)second / 3600;
		return temp;
	}

	public static double UT_decimal(long hour, long minute, long second, long timezone)
	{
		return get_range(decimal_hour(hour - timezone, minute, second));
	}

	public static Matrix e(double[] x_) // Eccentricity vector
	{
		Matrix e = new Matrix(3, 1);
		double x = x_[0], y = x_[1], z = x_[2], vx = x_[3], vy = x_[4], vz = x_[5];
		double r = Math.sqrt(x * x + y * y + z * z);

		e.set(0, 0, -x / r + (vy * vy * x + vz * vz * x - vx * vy * y - vx * vz * z) / mu);
		e.set(1, 0, -y / r + (-vx * vy * x + vx * vx * y + vz * vz * y - vy * vz * z) / mu);
		e.set(2, 0, -z / r + (-vx * vz * x - vy * vz * y + vx * vx * z + vy * vy * z) / mu);

		return e;
	}

	/**
	 * Eccentricity
	 * 
	 * @param r_ Position vector
	 * @param v_ Velocity vector
	 */
	public static Matrix e(Matrix r_, Matrix v_) // Eccentricity vector
	{
		Matrix e = new Matrix(3, 1);
		double x = r_.get(0, 0), y = r_.get(1, 0), z = r_.get(2, 0);
		double vx = v_.get(0, 0), vy = v_.get(1, 0), vz = v_.get(2, 0);
		double r = Math.sqrt(x * x + y * y + z * z);

		e.set(0, 0, -x / r + (vy * vy * x + vz * vz * x - vx * vy * y - vx * vz * z) / mu);
		e.set(1, 0, -y / r + (-vx * vy * x + vx * vx * y + vz * vz * y - vy * vz * z) / mu);
		e.set(2, 0, -z / r + (-vx * vz * x - vy * vz * y + vx * vx * z + vy * vy * z) / mu);

		return e;
	}

	public static Matrix h(Matrix r_, Matrix v_) // Angular momentum vector
	{
		Matrix h = new Matrix(3, 1);
		double x = r_.get(0, 0), y = r_.get(1, 0), z = r_.get(2, 0);
		double vx = v_.get(0, 0), vy = v_.get(1, 0), vz = v_.get(2, 0);

		h.set(0, 0, vz * y - vy * z);
		h.set(1, 0, -vz * x + vx * z);
		h.set(2, 0, vy * x - vx * y);

		return h;
	}

	public static Matrix h(double[] x_) // Angular momentum vector
	{
		Matrix h = new Matrix(3, 1);
		double x = x_[0], y = x_[1], z = x_[2], vx = x_[3], vy = x_[4], vz = x_[5];

		h.set(0, 0, vz * y - vy * z);
		h.set(1, 0, -vz * x + vx * z);
		h.set(2, 0, vy * x - vx * y);

		return h;
	}

	public static double p(double a, double e) // Semilatus Rectum
	{
		return a * (1 - e * e);
	}

	public static double Degree(double Rad) // radians -> degrees
	{
		return Rad * 180. / Math.PI;
	}

	/**
	 * Convert radians to degrees
	 * 
	 * @param Degree
	 */
	public static double Rad(double Degree) // degrees -> radians
	{
		return Degree * Math.PI / 180.;
	}

	public static Matrix Rot1(double angle) // Rotation matrix about x-axis
	{
		Matrix Rot1 = new Matrix(3, 3);
		Rot1.set(0, 0, 1.);
		Rot1.set(1, 0, 0.);
		Rot1.set(2, 0, 0.);
		Rot1.set(0, 1, 0.);
		Rot1.set(1, 1, Math.cos(angle));
		Rot1.set(2, 1, -Math.sin(angle));
		Rot1.set(0, 2, 0.);
		Rot1.set(1, 2, Math.sin(angle));
		Rot1.set(2, 2, Math.cos(angle));
		return Rot1;
	}

	public static Matrix Rot3(double angle) // Rotation matrix about z-axis
	{
		Matrix Rot3 = new Matrix(3, 3);
		Rot3.set(0, 0, Math.cos(angle));
		Rot3.set(1, 0, -Math.sin(angle));
		Rot3.set(2, 0, 0.);
		Rot3.set(0, 1, Math.sin(angle));
		Rot3.set(1, 1, Math.cos(angle));
		Rot3.set(2, 1, 0.);
		Rot3.set(0, 2, 0.);
		Rot3.set(1, 2, 0.);
		Rot3.set(2, 2, 1.);
		return Rot3;
	}

	public static Matrix r_from_el(double a, double e, double i, double Om, double om, double tau)
	{
		Matrix r_pqw = new Matrix(3, 1);
		double p = p(a, e);
		r_pqw.set(0, 0, p * Math.cos(tau) / (1 + e * Math.cos(tau)));
		r_pqw.set(1, 0, p * Math.sin(tau) / (1 + e * Math.cos(tau)));
		r_pqw.set(2, 0, 0.);

		Matrix r_ijk = Rot3(-Om).times(Rot1(-i).times(Rot3(-om).times(r_pqw)));

		return r_ijk;
	}

	public static Matrix v_from_el(double a, double e, double i, double Om, double om, double tau)
	{
		Matrix v_pqw = new Matrix(3, 1);
		double p = p(a, e);
		v_pqw.set(0, 0, -Math.sqrt(mu / p) * Math.sin(tau));
		v_pqw.set(1, 0, Math.sqrt(mu / p) * (e + Math.cos(tau)));
		v_pqw.set(2, 0, 0.);
		//myutil.print_Matrix(v_pqw,"v_pqw");

		Matrix v_ijk = Rot3(-Om).times(Rot1(-i).times(Rot3(-om).times(v_pqw)));

		return v_ijk;
	}

	public static Matrix r_from_el(Matrix orbit_el)
	{
		double a = orbit_el.get(0, 0), e = orbit_el.get(1, 0), i = orbit_el.get(2, 0);
		double Om = orbit_el.get(3, 0), om = orbit_el.get(4, 0), tau = orbit_el.get(5, 0);
		Matrix r_pqw = new Matrix(3, 1);
		double p = p(a, e);
		r_pqw.set(0, 0, p * Math.cos(tau) / (1 + e * Math.cos(tau)));
		r_pqw.set(1, 0, p * Math.sin(tau) / (1 + e * Math.cos(tau)));
		r_pqw.set(2, 0, 0.);

		Matrix r_ijk = Rot3(-Om).times(Rot1(-i).times(Rot3(-om).times(r_pqw)));

		return r_ijk;
	}

	public static Matrix v_from_el(Matrix orbit_el)
	{
		double a = orbit_el.get(0, 0), e = orbit_el.get(1, 0), i = orbit_el.get(2, 0);
		double Om = orbit_el.get(3, 0), om = orbit_el.get(4, 0), tau = orbit_el.get(5, 0);
		Matrix v_pqw = new Matrix(3, 1);
		double p = p(a, e);
		v_pqw.set(0, 0, -Math.sqrt(mu / p) * Math.sin(tau));
		v_pqw.set(1, 0, Math.sqrt(mu / p) * (e + Math.cos(tau)));
		v_pqw.set(2, 0, 0.);
		//myutil.print_Matrix(v_pqw,"v_pqw");

		Matrix v_ijk = Rot3(-Om).times(Rot1(-i).times(Rot3(-om).times(v_pqw)));

		return v_ijk;
	}

	public static void print_orbit(Matrix orbit_el, PrintWriter p_out, int steps)
	{
		double a = orbit_el.get(0, 0), e = orbit_el.get(1, 0), i = orbit_el.get(2, 0);
		double Om = orbit_el.get(3, 0), om = orbit_el.get(4, 0);
		double tau, delta_tau = 360. / steps;
		Matrix r;

		//System.out.println("x       y        z      ");
		for (tau = 0.; tau < 361; tau += delta_tau)
		{
			r = r_from_el(a, e, i, Om, om, Rad(tau));
			p_out.println("" + r.get(0, 0) + " " + r.get(1, 0) + " " + r.get(2, 0));
		}
	}

	public static void print_orbit(double a, double e, double i, double Om, double om, PrintWriter p_out, int steps)
	{
		double tau, delta_tau = 360. / steps;
		Matrix r;

		//System.out.println("x       y        z      ");
		for (tau = 0.; tau < 361; tau += delta_tau)
		{
			r = r_from_el(a, e, i, Om, om, Rad(tau));
			p_out.println("" + r.get(0, 0) + " " + r.get(1, 0) + " " + r.get(2, 0));
		}
	}

}

// Example from Vallado
/*
r_ijk=cm.r_from_el(5.664141370600592,0.83285,cm.Rad(87.87),cm.Rad(227.89),cm.Rad(53.38),cm.Rad(92.335));
v_ijk=cm.v_from_el(5.664141370600592,0.83285,cm.Rad(87.87),cm.Rad(227.89),cm.Rad(53.38),cm.Rad(92.335));
myutil.print_Matrix(r_ijk,"r_ijk");
myutil.print_Matrix(v_ijk,"v_ijk");
*/
/*
long RA_second;
long RA_minute;
long RA_hour;
double declination;
double RA;
double H;
double longitude;
double latitude;
double LST;
double GST;
double JD;
long timezone;
double UT;
*/
//long local_day,local_month,local_year,local_hour,local_minute,local_second;
/*
    // Mars Polar Lander
	jultime_depart = cm.juliandate(1999,1,3,12,0,0);
	jultime_arrive = cm.juliandate(1999,12,3,12,0,0);
    // Mars Climate Orbiter
	jultime_depart = cm.juliandate(1998,12,11,12,0,0);
	jultime_arrive = cm.juliandate(1999,9,23,12,0,0);
    // Mars Observer
	jultime_depart = cm.juliandate(1992,9,16,12,0,0);
	jultime_arrive = cm.juliandate(1993,8,4,12,0,0);
    // Mars Global Surveyor
	jultime_depart = cm.juliandate(1996,11,6,12,0,0);
	jultime_arrive = cm.juliandate(1997,9,11,12,0,0);
*/
/*
        x0[7]=1.7;
        x0[8]=-489.;
        x0[9]=0.015;
        x0[10]=14.17;
        x0[11]=-56.8;
        x0[12]=81.07;
        x0[13]=1.;      // lm
*/
/*
        x0[7]=-2;	// lx
        x0[8]=2;	// ly
        x0[9]=1;	// lz
        x0[10]=4;	// lvx
        x0[11]=5;	// lvy
        x0[12]=6;	// lvz
        x0[13]=-19.890212800451618;	// lm
*/