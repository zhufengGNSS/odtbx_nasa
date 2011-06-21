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

package jat.math;

import java.io.*;

/**
 * <P>
 * The MathUtils Class provides the mathematical constants and functions that
 * really should be in java.Math.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class MathUtils implements Serializable {
	/** PI.
	 */
	public final static double PI = Math.acos(-1.0);
	/** 2*PI.
	 */
	public final static double TWOPI = 2.0*Math.acos(-1.0);

	/** Conversion factor for degrees to radians.
	 */
	public final static double DEG2RAD = PI / 180.0;

	/** Conversion factor for arcsec to radians.
	 */
	public final static double ARCSEC2RAD = PI / 648000.0;

	/** Conversion factor for radians to degrees.
	 */
	public final static double RAD2DEG = 180.0 / PI;

	/** Machine precision for doubles.
	 */
	public final static double MACHEPS = 1.1102230246251565E-16;

	/** Creates a new instance of math */
	public MathUtils() {
	}

	/** Fractional part of a number (y=x-[x])
	 * @return fractional part of a number
	 */
	public static double Frac(double x) {
		return x - Math.floor(x);
	}

	/** Modulo function
	 * @return the Modulo
	 */
	public static double Modulo(double x, double y) {
		return y * Frac(x / y);
	}

	/** mod function
	 * @param x Numerator.
	 * @param y Denominator.
	 * @return floating point remainder of (x/y)
	 */
	public static double mod(double x, double y) {
		int temp = (int) (x / y);
		double out = x - temp * y;
		return out;
	}

	/** Convert Degrees:Minutes:Seconds to Degrees.
	 * @param degrees Integer degrees.
	 * @param minutes Integer minutes.
	 * @param seconds seconds.
	 * @return real number of degrees.
	 */
	public static double dms2deg(int degrees, int minutes, double seconds) {
		double out = degrees + (minutes / 60.0) + (seconds / 3600.0);
		return out;
	}

	/** Convert Degrees:Minutes:Seconds to Radians.
	 * @param degrees Integer degrees.
	 * @param minutes Integer minutes.
	 * @param seconds seconds.
	 * @return real number of radians.
	 */
	public static double dms2rad(int degrees, int minutes, double seconds) {
		double out = dms2deg(degrees, minutes, seconds) * DEG2RAD;
		return out;
	}

	/** Returns sign of double number.
	 * @param z double number.
	 * @return int sign.
	 */
	public static int signum(double z) {
		return (z < 0) ? -1 : 1;
	}

	/** Returns absolute value of double number.
	 * @param z double number.
	 * @return int sign.
	 */
	public static double abs(double z) {
		return (z < 0) ? -z : z;
	}

	/** Returns absolute value of a times the sign of b.
	 * @param a double
	 * @param b double number.
	 * @return absolute value of a times the sign of b.
	 */

	public static double sign(double a, double b) {
		double amag = Math.abs(a);
		double out = amag;
		if (b < 0.0)
			out = -1.0 * amag;
		return out;
	}

	/** Returns absolute value of a times the sign of b.
	 * @param a integer
	 * @param b integer.
	 * @return absolute value of a times the sign of b.
	 */
	public static int sign(int a, int b) {
		int amag = Math.abs(a);
		int out = amag;
		if (b < 0.0)
			out = -1 * amag;
		return out;
	}

	/** Returns a double rounded to the nearest integer.
	 * @param a double to be rounded
	 * @return double containing the rounded number.
	 */
	public static double round(double a) {
		return Math.floor(a + 0.5);
	}

	/** Returns a double rounded to the nearest integer.
	 * @param a double to be rounded
	 * @return int containing the rounded number.
	 */
	public static int roundInt(double a) {
		return (int) Math.floor(a + 0.5);
	}

	/** Returns a double rounded to the nearest integer.
	 * @param a double to be rounded
	 * @return long containing the rounded number.
	 */
	public static long roundLong(double a) {
		return (long) Math.floor(a + 0.5);
	}

	/** Returns n!
	 * @param n an int
	 * @return n!
	 */
	public static int factorial(int n) {
		int returnValue;
		if (n < 0) {
			throw new IllegalArgumentException("Illegal value: " + n);
		} else if ((n == 0)||(n == 1)) {
			returnValue = 1;
		} else	{
			returnValue = n * factorial(n - 1);
		}
		return returnValue;
	}
	/** Returns n!
	 * @param n a long
	 * @return n!
	 */
	public static double factorial(double n) {
		double returnValue;
		if (n < 0.0) {
			throw new IllegalArgumentException("Illegal value: " + n);
		} else if ((n == 0)||(n == 1)) {
			returnValue = 1.0;
		} else	{
			returnValue = n * factorial(n - 1.0);
		}
		return returnValue;
	}
	
	/** Returns x+y assuming that both x and y have length: x.length
	 * @param x double array
	 * @param y double array
	 * @return x+y double array
	 */
	public static double[] plus(double[] x, double[] y){
	    double[] out = new double[x.length];
	    for(int i=0; i<x.length; i++){
	        out[i] = x[i]+y[i];
	    }
	    return out;
	}
	/** Returns an angle between 0 and 2 pi
	 * @param x double containing angle in radians
	 * @return an angle between 0 and 2 pi
	 */
	
	public static double radiancheck(double x){
		x = Modulo(x, TWOPI);
		if (x < 0.0) x = x + TWOPI;
		return x;
	}
	/** Returns an angle between 0 and 360
	 * @param x double containing angle in degrees
	 * @return an angle between 0 and 360
	 */
	
	public static double degreecheck(double x){
		x = Modulo(x, 360.0);
		if (x < 0.0) x = x + 360.0;
		return x;
	}

	/** Returns a quadrant-checked angle between 0 and 2PI
	 * @param x double containing length of opposite side
	 * @param r double containing length of hypotenuse
	 * @return an angle between 0 and 2 pi in radians
	 */	
	public static double arcsin(double x, double r){
		double phi = Math.asin(x/r);
		if (x < 0.0) phi = MathUtils.TWOPI - phi;
		return phi;		
	}
	
	public static void main (String[] args) {
		double rem = MathUtils.mod(18.5, 4.2);
		System.out.println(rem);
		
	}
}