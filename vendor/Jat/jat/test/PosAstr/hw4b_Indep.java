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

public class hw4b_Indep
{
	public static void main(String[] args)
	{
		boolean plus=true,minus=false;
		Angle alpha0, delta0;
		double tau1, tau2, tau;
		double factor;
		Angle mua, mud; // mu alpha, mu delta
		Angle f1, g1, G1, h1, H1, i1; //f, etc on July 1
		Angle f2, g2, G2, h2, H2, i2; //f, etc on July 2
		Angle f, g, G, h, H, i; // interpolated value

		// Alpha Tauri
		System.out.println("Aldebaran June 3, 1980 Independent Day Numbers");
		System.out.println("----------------------------------------------");
		alpha0 = new Angle(true, 4, 34, 46., Angle.HOURANGLE);
		delta0 = new Angle(true, 16, 28, 12., Angle.ARCDEGREES);
		mua = new Angle(true, 0, 0, 0.439, Angle.HOURANGLE);
		mud = new Angle(false, 0, 0, 18.97, Angle.ARCDEGREES);

		// Independent day number table values
		tau1 = 0.4211;
		tau2 = 0.4239;
		//factor = 0.871;
		factor = 0.;
		f1 = new Angle(plus, 0, 0, 0.6251, Angle.HOURANGLE);
		f2 = new Angle(plus, 0, 0, 0.6368, Angle.HOURANGLE);
		g1 = new Angle(plus, 0, 0, 8.934, Angle.ARCDEGREES);
		g2 = new Angle(plus, 0, 0, 8.943, Angle.ARCDEGREES);
		G1 = new Angle(plus, 4,11,8., Angle.HOURANGLE);
		G2 = new Angle(plus, 4,9,3, Angle.HOURANGLE);
		h1 = new Angle(plus, 0,0,20.345, Angle.ARCDEGREES);
		h2 = new Angle(plus, 0, 0, 20.362, Angle.ARCDEGREES);
		H1 = new Angle(plus, 13,4,15, Angle.HOURANGLE);
		H2 = new Angle(plus, 13,0,42., Angle.HOURANGLE);
		i1 = new Angle(minus, 0, 0, 2.441, Angle.ARCDEGREES);
		i2 = new Angle(minus, 0, 0, 2.311, Angle.ARCDEGREES);

		// Interpolate table values
		f = new Angle(f1.radians + factor * (f2.radians - f1.radians), Angle.RADIANS);
		g = new Angle(g1.radians + factor * (g2.radians - g1.radians), Angle.RADIANS);
		G = new Angle(G1.radians + factor * (G2.radians - G1.radians), Angle.RADIANS);
		h = new Angle(h1.radians + factor * (h2.radians - h1.radians), Angle.RADIANS);
		H = new Angle(H1.radians + factor * (H2.radians - H1.radians), Angle.RADIANS);
		i = new Angle(i1.radians + factor * (i2.radians - i1.radians), Angle.RADIANS);
		tau = tau1 + factor * (tau2 - tau1);

		Position pos = StarAlg.mean_to_apparent_Independent(alpha0, delta0, tau, mua, mud, f, g, G, h, H, i, true);
		alpha0.println("alpha0  = ", Angle.SHA);
		pos.RA.println("alpha  = ", Angle.SHA);
		//pos.RA.println("alpha  = ", Angle.HOURANGLE);
		//pos.dec.println("delta  = ", Angle.ARCDEGREES);
	}
}
