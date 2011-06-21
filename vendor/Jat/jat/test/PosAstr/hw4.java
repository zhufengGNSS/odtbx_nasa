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


package jat.test.PosAstr;

import jat.cm.*;

public class hw4
{
	public static void main(String[] args)
	{
		Angle m1960, n1960, m1980, n1980, alpha0, delta0;
		double[] tau_interval = new double[2]; // tau on bracketting days
		Angle[] mu = new Angle[2]; // mu alpha, mu delta
		Angle[] ABCDE_interval = new Angle[10]; // A, B, C, D, E on bracketting days
		Angle[] fgGhHi_interval=new Angle[12];

		// m and n for 1960 and 1980
		m1960 = new Angle(true, 0, 0, 3.07346, Angle.HOURANGLE);
		n1960 = new Angle(true, 0, 0, 1.33612, Angle.HOURANGLE);

		// Alpha Coronae Borealis July 1, 1960
		System.out.println("Alpha Coronae Borealis July 1, 1960, Besselian Day Numbers");
		System.out.println("----------------------------------------------------------");
		alpha0 = new Angle(true, 15, 32, 59.556, Angle.HOURANGLE);
		delta0 = new Angle(true, 26, 50, 53.77, Angle.ARCDEGREES);
		tau_interval[0] = 0.4974;
		tau_interval[1] = 0.5001;
		mu[0] = new Angle(true, 0, 0, 0.009, Angle.HOURANGLE);
		mu[1] = new Angle(false, 0, 0, 0.09, Angle.ARCDEGREES);
		ABCDE_interval[0] = new Angle(true, 0, 0, 8.786, Angle.ARCDEGREES); // A on July 1
		ABCDE_interval[1] = new Angle(true, 0, 0, 8.818, Angle.ARCDEGREES); // A on July 2
		ABCDE_interval[2] = new Angle(true, 0, 0, 9.580, Angle.ARCDEGREES); // B on July 1
		ABCDE_interval[3] = new Angle(true, 0, 0, 9.572, Angle.ARCDEGREES); // B on July 2
		ABCDE_interval[4] = new Angle(true, 0, 0, 2.994, Angle.ARCDEGREES); // C on July 1
		ABCDE_interval[5] = new Angle(true, 0, 0, 3.302, Angle.ARCDEGREES); // C on July 2
		ABCDE_interval[6] = new Angle(false, 0, 0, 20.199, Angle.ARCDEGREES); // D on July 1
		ABCDE_interval[7] = new Angle(false, 0, 0, 20.140, Angle.ARCDEGREES); // D on July 2
		ABCDE_interval[8] = new Angle(false, 0, 0, 0.0004, Angle.HOURANGLE); // E on July 1
		ABCDE_interval[9] = new Angle(false, 0, 0, 0.0005, Angle.HOURANGLE); // E on July 2
		mean_to_apparent_interval_Besselian(m1960, n1960, alpha0, delta0, 0.871, tau_interval, mu, ABCDE_interval);

	}

	private static void mean_to_apparent_interval_Besselian(
		Angle m,
		Angle n,
		Angle alpha0,
		Angle delta0,
		double factor,
		double[] tau_interval,
		Angle[] mu,
		Angle[] ABCDE_interval)
	{
		Angle[] ABCDE = new Angle[5];
		double tau;

		// Interpolate table values
		ABCDE[0] = new Angle(ABCDE_interval[0].radians + factor * (ABCDE_interval[1].radians - ABCDE_interval[0].radians), Angle.RADIANS);
		ABCDE[1] = new Angle(ABCDE_interval[2].radians + factor * (ABCDE_interval[3].radians - ABCDE_interval[2].radians), Angle.RADIANS);
		ABCDE[2] = new Angle(ABCDE_interval[4].radians + factor * (ABCDE_interval[5].radians - ABCDE_interval[4].radians), Angle.RADIANS);
		ABCDE[3] = new Angle(ABCDE_interval[6].radians + factor * (ABCDE_interval[7].radians - ABCDE_interval[6].radians), Angle.RADIANS);
		ABCDE[4] = new Angle(ABCDE_interval[8].radians + factor * (ABCDE_interval[9].radians - ABCDE_interval[8].radians), Angle.RADIANS);
		tau = tau_interval[0] + factor * (tau_interval[1] - tau_interval[0]);
		mean_to_apparent_Besselian(m, n, alpha0, delta0, tau, mu, ABCDE);
	}

	static void mean_to_apparent_Besselian(Angle m, Angle n, Angle alpha0, Angle delta0, double tau, Angle[] mu, Angle[] ABCDE)
	{
		Angle alpha, delta, A, B, C, D, E;
		Angle Aa, Bb, Cc, Dd, Aaprime, Bbprime, Ccprime, Ddprime;
		double a, b, c, d;
		double aprime, bprime, cprime, dprime;
		double taumualpha, taumudelta;
		double sina0 = Math.sin(alpha0.radians);
		double cosa0 = Math.cos(alpha0.radians);
		double sind0 = Math.sin(delta0.radians);
		double cosd0 = Math.cos(delta0.radians);
		double tand0 = Math.tan(delta0.radians);
		Angle epsilon = new Angle(true, 23, 26, 45., Angle.ARCDEGREES);

		taumualpha = tau * mu[0].radians;
		taumudelta = tau * mu[1].radians;
		A = ABCDE[0];
		B = ABCDE[1];
		C = ABCDE[2];
		D = ABCDE[3];
		E = ABCDE[4];
		a = m.radians / n.radians + sina0 * Math.tan(delta0.radians);
		b = cosa0 * Math.tan(delta0.radians);
		c = cosa0 * 1. / cosd0;
		d = sina0 * 1. / cosd0;
		aprime = cosa0;
		bprime = -sina0;
		cprime = Math.tan(epsilon.radians) * cosd0 - sina0 * sind0;
		dprime = cosa0 * sind0;
		Aa = new Angle(A.radians * a, Angle.RADIANS);
		Bb = new Angle(B.radians * b, Angle.RADIANS);
		Cc = new Angle(C.radians * c, Angle.RADIANS);
		Dd = new Angle(D.radians * d, Angle.RADIANS);
		Aaprime = new Angle(A.radians * aprime, Angle.RADIANS);
		Bbprime = new Angle(B.radians * bprime, Angle.RADIANS);
		Ccprime = new Angle(C.radians * cprime, Angle.RADIANS);
		Ddprime = new Angle(D.radians * dprime, Angle.RADIANS);
		alpha = new Angle(alpha0.radians + taumualpha + Aa.radians + Bb.radians + Cc.radians + Dd.radians + E.radians, Angle.RADIANS);
		delta = new Angle(delta0.radians + taumudelta + Aaprime.radians + Bbprime.radians + Ccprime.radians + Ddprime.radians, Angle.RADIANS);

		System.out.println("tau = " + tau);
		A.println("A", Angle.ARCDEGREES);
		B.println("B", Angle.ARCDEGREES);
		C.println("C", Angle.ARCDEGREES);
		D.println("D", Angle.ARCDEGREES);
		E.println("E", Angle.HOURANGLE);
		System.out.println("a = " + a);
		System.out.println("b = " + b);
		System.out.println("c = " + c);
		System.out.println("d = " + d);
		System.out.println("a' = " + aprime);
		System.out.println("b' = " + bprime);
		System.out.println("c' = " + cprime);
		System.out.println("d' = " + dprime);
		Aa.println("Aa", Angle.HOURANGLE);
		Bb.println("Bb", Angle.HOURANGLE);
		Cc.println("Cc", Angle.HOURANGLE);
		Dd.println("Dd", Angle.HOURANGLE);
		Aaprime.println("Aa'", Angle.ARCDEGREES);
		Bbprime.println("Bb'", Angle.ARCDEGREES);
		Ccprime.println("Cc'", Angle.ARCDEGREES);
		Ddprime.println("Dd'", Angle.ARCDEGREES);
		alpha0.println("alpha0 = ", Angle.HOURANGLE);
		alpha.println("alpha  = ", Angle.HOURANGLE);
		delta0.println("delta0 = ", Angle.ARCDEGREES);
		delta.println("delta  = ", Angle.ARCDEGREES);
	}
}
