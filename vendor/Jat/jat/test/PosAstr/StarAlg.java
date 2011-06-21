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

/**
 * Star Algorithms
 * @author Tobias Berthold
 *
 */
public class StarAlg
{
	/** Reduction from mean to apparent place by Besselian day numbers
	 * @param m
	 * @param n
	 * @param alpha0
	 * @param delta0
	 * @param tau
	 * @param mua
	 * @param mud
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @param E
	 * @param print
	 * @return New position
	 */
	public static Position mean_to_apparent_Besselian(
		Angle m,
		Angle n,
		Angle alpha0,
		Angle delta0,
		double tau,
		Angle mua,
		Angle mud,
		Angle A,
		Angle B,
		Angle C,
		Angle D,
		Angle E,
		boolean print)
	{
		Angle alpha, delta;
		Angle Aa, Bb, Cc, Dd, Aaprime, Bbprime, Ccprime, Ddprime;
		double a, b, c, d;
		double aprime, bprime, cprime, dprime;
		Angle taumualpha, taumudelta;
		double sina0 = Math.sin(alpha0.radians);
		double cosa0 = Math.cos(alpha0.radians);
		double sind0 = Math.sin(delta0.radians);
		double cosd0 = Math.cos(delta0.radians);
		double tand0 = Math.tan(delta0.radians);
		Angle epsilon = new Angle(true, 23, 26, 34.06, Angle.ARCDEGREES);

		// Computations
		taumualpha = new Angle(tau * mua.radians, Angle.RADIANS);
		taumudelta = new Angle(tau * mud.radians, Angle.RADIANS);
		a = m.radians / n.radians + sina0 * tand0;
		b = cosa0 * tand0;
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
		alpha = new Angle(alpha0.radians + taumualpha.radians + Aa.radians + Bb.radians + Cc.radians + Dd.radians + E.radians, Angle.RADIANS);
		delta = new Angle(delta0.radians + taumudelta.radians + Aaprime.radians + Bbprime.radians + Ccprime.radians + Ddprime.radians, Angle.RADIANS);

		//Printing
		if (print)
		{
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
			taumualpha.println("tau*mualpha = ", Angle.HOURANGLE);
			taumudelta.println("tau*mudelta = ", Angle.ARCDEGREES);
			alpha0.println("alpha0 = ", Angle.HOURANGLE);
			alpha.println("alpha  = ", Angle.HOURANGLE);
			delta0.println("delta0 = ", Angle.ARCDEGREES);
			delta.println("delta  = ", Angle.ARCDEGREES);
		}

		Position pos = new Position(alpha, delta);
		return pos;
	}

	/** Reduction from mean to apparent place by independent day numbers
	 * @param alpha0
	 * @param delta0
	 * @param tau
	 * @param mua
	 * @param mud
	 * @param f
	 * @param g
	 * @param G
	 * @param h
	 * @param H
	 * @param i
	 * @param print
	 * @return new position
	 */
	public static Position mean_to_apparent_Independent(
		Angle alpha0,
		Angle delta0,
		double tau,
		Angle mua,
		Angle mud,
		Angle f,
		Angle g,
		Angle G,
		Angle h,
		Angle H,
		Angle i,
		boolean print)
	{
		Angle alpha, delta;
		Angle gaterm, haterm, gdterm, hdterm, idterm;
		Angle taumualpha, taumudelta;
		double sina0 = Math.sin(alpha0.radians);
		double cosa0 = Math.cos(alpha0.radians);
		double sind0 = Math.sin(delta0.radians);
		double cosd0 = Math.cos(delta0.radians);
		double tand0 = Math.tan(delta0.radians);
		//Angle epsilon = new Angle(true, 23, 26, 34.06, Angle.ARCDEGREES);

		// Computations
		taumualpha = new Angle(tau * mua.radians, Angle.RADIANS);
		taumudelta = new Angle(tau * mud.radians, Angle.RADIANS);
		gaterm = new Angle(g.radians * Math.sin(G.radians + alpha0.radians) * tand0, Angle.RADIANS);
		haterm = new Angle(h.radians * Math.sin(H.radians + alpha0.radians) / cosd0, Angle.RADIANS);
		gdterm = new Angle(g.radians * Math.cos(G.radians + alpha0.radians), Angle.RADIANS);
		hdterm = new Angle(h.radians * Math.cos(H.radians + alpha0.radians) * sind0, Angle.RADIANS);
		idterm = new Angle(i.radians * cosd0, Angle.RADIANS);
		alpha = new Angle(alpha0.radians + taumualpha.radians + f.radians + gaterm.radians + haterm.radians, Angle.RADIANS);
		delta = new Angle(delta0.radians + taumudelta.radians + gdterm.radians + hdterm.radians + idterm.radians, Angle.RADIANS);

		//Printing
		if (print)
		{
			System.out.println("tau = " + tau);
			f.println("f", Angle.HOURANGLE);
			g.println("g", Angle.ARCDEGREES);
			G.println("G", Angle.HOURANGLE);
			h.println("h", Angle.ARCDEGREES);
			H.println("H", Angle.HOURANGLE);
			i.println("i", Angle.ARCDEGREES);
			gaterm.println("g*sin(G+alpha0)*tan(delta0)", Angle.HOURANGLE);
			haterm.println("h*sin(H+alpha0)*sec(delta0)", Angle.HOURANGLE);
			gdterm.println("g*cos(G+alpha0)            ", Angle.ARCDEGREES);
			hdterm.println("h*cos(H+alpha0)*sin(delta0)", Angle.ARCDEGREES);
			idterm.println("i*cos(delta0)              ", Angle.ARCDEGREES);
			taumualpha.println("tau*mualpha = ", Angle.HOURANGLE);
			taumudelta.println("tau*mudelta = ", Angle.ARCDEGREES);
			alpha0.println("alpha0 = ", Angle.HOURANGLE);
			alpha.println("alpha  = ", Angle.HOURANGLE);
			delta0.println("delta0 = ", Angle.ARCDEGREES);
			delta.println("delta  = ", Angle.ARCDEGREES);
		}
		Position pos = new Position(alpha, delta);
		return pos;
	}

	public static Position annual_precession(Angle m, Angle n, Angle alpha0, Angle delta0, boolean print)
	{
		double delta_alpha, delta_delta;
		Angle alpha, delta;

		delta_alpha = m.radians + n.radians * Math.sin(alpha0.radians) * Math.tan(delta0.radians);
		delta_delta = n.radians * Math.cos(alpha0.radians);
		//new Angle(delta_alpha, Angle.RADIANS).println("Delta alpha");
		//new Angle(delta_delta, Angle.RADIANS).println("Delta delta");
		alpha = new Angle(alpha0.radians + delta_alpha, Angle.RADIANS);
		delta = new Angle(delta0.radians + delta_delta, Angle.RADIANS);

		//Printing
		if (print)
		{
			alpha0.println("alpha0", Angle.HOURANGLE);
			alpha.println("alpha ", Angle.HOURANGLE);
			delta0.println("delta0", Angle.ARCDEGREES);
			delta.println("delta ", Angle.ARCDEGREES);
		}
		Position pos = new Position(alpha, delta);
		return pos;

	}
}