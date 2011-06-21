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

import jat.math.MathUtils;

/**
 * @author Tobias Berthold
 * internally uses radians
 */
public class Angle
{
	public final static int RADIANS = 1, DEGREES = 2, DECIMALHOURS = 3;
	public final static int ARCDEGREES = 4, HOURANGLE = 5, SHA=6;
	public double radians;
	double degrees;
	double hours;
	HourAngle h = new HourAngle();
	ArcDegrees a = new ArcDegrees();
	ArcDegrees sha = new ArcDegrees();

	public Angle(double angle, int mode)
	{
		if (mode == RADIANS)
			radians = angle;
		if (mode == DEGREES)
			radians = Math.toRadians(angle);
		if (mode == DECIMALHOURS)
			radians = Math.toRadians(angle * 15.);
		convert();
	}

	public Angle(boolean positive, int hours_arcdegrees, int minutes, double seconds, int mode)
	{
		h.positive = positive;
		a.positive = positive;
		if (mode == ARCDEGREES)
			radians = Math.toRadians(hours_arcdegrees + minutes / 60. + seconds / 3600.);
		if (mode == HOURANGLE)
			radians = Math.toRadians(15. * (hours_arcdegrees + minutes / 60. + seconds / 3600.));
		if (!positive)
			radians = -radians;
		convert();
	}

	/**
	 * assumes that angle is store in radians
	 */
	void convert()
	{
		double tmpdegrees;
		double shadegrees;
		
		degrees = Math.toDegrees(radians);
		hours = degrees / 15.;

		if (radians > 0.)
		{
			a.positive=true;
			sha.positive=true;
			h.positive=true;
			tmpdegrees = degrees;
			shadegrees = 360.-degrees;
		}
		else
		{
			a.positive=false;
			sha.positive=false;
			h.positive=false;
			tmpdegrees = -degrees;
			shadegrees = -degrees;
		}
		// Angle expressed in degrees, minutes, seconds
		a.degrees = (int) tmpdegrees;
		a.minutes = (int) (MathUtils.Frac(tmpdegrees) * 60.0);
		a.seconds = (tmpdegrees - a.degrees - a.minutes / 60.) * 3600.;

		// Angle expressed as SHA in degrees, minutes, seconds
		//double shadegrees=360.-degrees;
		sha.degrees = (int) shadegrees;
		sha.minutes = (int) (MathUtils.Frac(shadegrees) * 60.0);
		sha.seconds = (shadegrees - sha.degrees - sha.minutes / 60.) * 3600.;
		
		// Angle expressed as hour angle in hours, minutes, seconds
		double tmphours=tmpdegrees/15.;
		h.hours = (int) tmphours;
		//System.out.println(decimalhours+ " decimal hours");
		double tmpminutes = 60. * MathUtils.Frac(tmphours);
		//System.out.println(decimalminutes+ " decimal minutes");
		h.minutes = (int) (60. * MathUtils.Frac(tmphours));
		h.seconds = 60. * MathUtils.Frac(tmpminutes);
	}

	public Angle add(Angle b)
	{
		Angle result = new Angle(radians + b.radians, RADIANS);
		result.convert();
		return result;
	}

	public double sin(Angle a)
	{
		return Math.sin(a.radians);
	}

	public double tan(Angle a)
	{
		return Math.tan(a.radians);
	}

	public double cos(Angle a)
	{
		return Math.cos(a.radians);
	}

	public void print()
	{
		System.out.println(radians + " Radians");
		System.out.println(degrees + " Decimal degrees");
		if (a.positive)
			System.out.print("+");
		else
			System.out.print("-");
		System.out.println(a.degrees + " degrees " + a.minutes + " minutes " + a.seconds + " seconds");
		if (h.positive)
			System.out.print("+");
		else
			System.out.print("-");
		System.out.println(h.hours + " hours " + h.minutes + " minutes " + h.seconds + " seconds");
		//System.out.println(decimalhours + " Decimal hours");
	}

	public void println()
	{
		print();
		System.out.println();
	}

	public void println(String title)
	{
		System.out.println(title);
		println();
	}

	public void println(String title, int mode)
	{
		System.out.print(title+" = ");
		switch (mode)
		{
			case RADIANS :
				System.out.println(radians + " Radians");
				break;
			case DEGREES :
				System.out.println(degrees + " Decimal degrees");
				break;
			case DECIMALHOURS :
				System.out.println(hours + " Decimal hours");
				break;
			case ARCDEGREES :
				if (a.positive)
					System.out.print("+");
				else
					System.out.print("-");
				System.out.println(a.degrees + " Degrees " + a.minutes + " Minutes " + a.seconds + " Seconds");
				break;
			case SHA :
				if (sha.positive)
					System.out.print("+");
				else
					System.out.print("-");
				System.out.println(sha.degrees + " Degrees " + sha.minutes + " Minutes " + sha.seconds + " Seconds");
				break;
			case HOURANGLE :
				if (h.positive)
					System.out.print("+");
				else
					System.out.print("-");
				System.out.println(h.hours + " Hours " + h.minutes + " Minutes " + h.seconds + " Seconds");
				break;
		}
	}

	public static void main(String[] args)
	{
		System.out.println("Angle class example");
		//Angle alpha = new Angle(1.5, Angle.RADIANS);
		Angle alpha = new Angle(59.51, Angle.DEGREES);
		alpha.println();

		System.out.println(1002.23 + " arcseconds ");
		Angle beta = new Angle(true, 0, 0, 1002.23, Angle.ARCDEGREES);
		//		Angle beta = new Angle(0,16,42.23, Angle.ARCDEGREES);
		//		Angle beta = new Angle(0,0,20.04, Angle.ARCDEGREES);
		beta.println();

		// Example from Duffett-Smith
		System.out.println(18.52417 + " decimal hours");
		Angle gamma = new Angle(18.52417, Angle.DECIMALHOURS);
		gamma.println();

		System.out.println("18 h 31 m 27.012 s");
		Angle delta = new Angle(true, 18, 31, 27.012, Angle.HOURANGLE);
		delta.println();

		System.out.println("");
		Angle n_1960 = new Angle(true, 0, 0, 1.33612, Angle.HOURANGLE);
		n_1960.println();

		Angle a1 = new Angle(40., Angle.DEGREES);
		Angle a2 = new Angle(45., Angle.DEGREES);
		a1.add(a2).println();
	}
}
