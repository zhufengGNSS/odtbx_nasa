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

public class hw4b_Bessel
{
	public static void main(String[] args)
	{
		boolean plus=true,minus=false;
		Angle m, n, alpha0, delta0;
		Angle A1, B1, C1, D1, E1; //A, etc on June 3
		Angle A2, B2, C2, D2, E2; //A, etc on June 3
		Angle A, B, C, D, E; // interpolated value
		double tau1, tau2, tau;
		//double factor = 0.871;
		double factor = 0.;
		Angle mua, mud; // mu alpha, mu delta

		// m and n
		m = new Angle(true, 0, 0, 3.07383, Angle.HOURANGLE);
		n = new Angle(true, 0, 0, 1.33600, Angle.HOURANGLE);

		System.out.println("Aldebaran June 3, 1980 Besselian Day Numbers");
		System.out.println("---------------------------------------------");
		alpha0 = new Angle(true, 4, 34, 46., Angle.HOURANGLE);
		delta0 = new Angle(true, 16, 28, 12., Angle.ARCDEGREES);
		mua = new Angle(true, 0, 0, 0.439, Angle.HOURANGLE);
		mud = new Angle(false, 0, 0, 18.97, Angle.ARCDEGREES);

		// Besselian day number table values
		tau1 = 0.4211;
		tau2 = 0.4239;
		A1 = new Angle(plus, 0, 0, 4.086, Angle.ARCDEGREES);
		A2 = new Angle(true, 0, 0, 4.162, Angle.ARCDEGREES);
		B1 = new Angle(true, 0, 0, 7.945, Angle.ARCDEGREES);
		B2 = new Angle(true, 0, 0, 7.915, Angle.ARCDEGREES);
		C1 = new Angle(false, 0, 0, 5.629, Angle.ARCDEGREES);
		C2 = new Angle(false, 0, 0, 5.330, Angle.ARCDEGREES);
		D1 = new Angle(false, 0, 0, 19.551, Angle.ARCDEGREES);
		D2 = new Angle(false, 0, 0, 19.652, Angle.ARCDEGREES);
		E1 = new Angle(false, 0, 0, 0.0016, Angle.HOURANGLE);
		E2 = new Angle(false, 0, 0, 0.0016, Angle.HOURANGLE);

		// Interpolate table values
		A = new Angle(A1.radians + factor * (A2.radians - A1.radians), Angle.RADIANS);
		B = new Angle(B1.radians + factor * (B2.radians - B1.radians), Angle.RADIANS);
		C = new Angle(C1.radians + factor * (C2.radians - C1.radians), Angle.RADIANS);
		D = new Angle(D1.radians + factor * (D2.radians - D1.radians), Angle.RADIANS);
		E = new Angle(E1.radians + factor * (E2.radians - E1.radians), Angle.RADIANS);
		tau = tau1 + factor * (tau2 - tau1);

		Position pos = StarAlg.mean_to_apparent_Besselian(m, n, alpha0, delta0, tau, mua, mud, A, B, C, D, E, true);
		alpha0.println("alpha0  = ", Angle.SHA);
		pos.RA.println("alpha  = ", Angle.SHA);
		//pos.dec.println("delta  = ", Angle.ARCDEGREES);
	}
}
