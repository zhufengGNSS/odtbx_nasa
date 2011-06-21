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

public class hw4a_Bessel
{
	public static void main(String[] args)
	{
		Angle m1960, n1960, alpha0, delta0;
		Angle A1, B1, C1, D1, E1; //A, etc on July 1
		Angle A2, B2, C2, D2, E2; //A, etc on July 2
		Angle A, B, C, D, E; // interpolated value
		double tau1, tau2, tau;
		double factor = 0.871;
		Angle mua, mud; // mu alpha, mu delta

		// m and n for 1960 and 1980
		m1960 = new Angle(true, 0, 0, 3.07346, Angle.HOURANGLE);
		n1960 = new Angle(true, 0, 0, 1.33612, Angle.HOURANGLE);

		// Alpha Coronae Borealis July 1, 1960
		System.out.println("Alpha Coronae Borealis July 1, 1960, Besselian Day Numbers");
		System.out.println("----------------------------------------------------------");
		alpha0 = new Angle(true, 15, 32, 59.556, Angle.HOURANGLE);
		delta0 = new Angle(true, 26, 50, 53.77, Angle.ARCDEGREES);
		mua = new Angle(true, 0, 0, 0.009, Angle.HOURANGLE);
		mud = new Angle(false, 0, 0, 0.09, Angle.ARCDEGREES);

		// Besselian day number table values
		tau1 = 0.4974;
		tau2 = 0.5001;
		A1 = new Angle(true, 0, 0, 8.786, Angle.ARCDEGREES);
		A2 = new Angle(true, 0, 0, 8.818, Angle.ARCDEGREES);
		B1 = new Angle(true, 0, 0, 9.580, Angle.ARCDEGREES);
		B2 = new Angle(true, 0, 0, 9.572, Angle.ARCDEGREES);
		C1 = new Angle(true, 0, 0, 2.994, Angle.ARCDEGREES);
		C2 = new Angle(true, 0, 0, 3.302, Angle.ARCDEGREES);
		D1 = new Angle(false, 0, 0, 20.199, Angle.ARCDEGREES);
		D2 = new Angle(false, 0, 0, 20.140, Angle.ARCDEGREES);
		E1 = new Angle(false, 0, 0, 0.0004, Angle.HOURANGLE);
		E2 = new Angle(false, 0, 0, 0.0005, Angle.HOURANGLE);

		// Interpolate table values
		A = new Angle(A1.radians + factor * (A2.radians - A1.radians), Angle.RADIANS);
		B = new Angle(B1.radians + factor * (B2.radians - B1.radians), Angle.RADIANS);
		C = new Angle(C1.radians + factor * (C2.radians - C1.radians), Angle.RADIANS);
		D = new Angle(D1.radians + factor * (D2.radians - D1.radians), Angle.RADIANS);
		E = new Angle(E1.radians + factor * (E2.radians - E1.radians), Angle.RADIANS);
		tau = tau1 + factor * (tau2 - tau1);

		Position pos=StarAlg.mean_to_apparent_Besselian(m1960, n1960, alpha0, delta0, tau, mua, mud, A, B, C, D, E,true);
		//pos.RA.println("alpha  = ", Angle.HOURANGLE);
		//pos.dec.println("delta  = ", Angle.ARCDEGREES);
	}
}
