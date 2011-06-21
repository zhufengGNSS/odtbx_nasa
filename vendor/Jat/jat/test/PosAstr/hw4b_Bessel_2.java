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

public class hw4b_Bessel_2
{
	public static void main(String[] args)
	{
		boolean plus=true,minus=false;
		Angle m, n, alpha0, delta0;
		Angle A, B, C, D, E; // tabulated day numbers
		double tau;
		Angle mua, mud; // mu alpha, mu delta

		// m and n
		m = new Angle(true, 0, 0, 3.07383, Angle.HOURANGLE);
		n = new Angle(true, 0, 0, 1.33600, Angle.HOURANGLE);

		System.out.println("Aldebaran June 3, 1980 Besselian Day Numbers");
		System.out.println("---------------------------------------------");
		alpha0 = new Angle(plus, 4, 34, 46., Angle.HOURANGLE);
		delta0 = new Angle(plus, 16, 28, 12., Angle.ARCDEGREES);
		mua = new Angle(plus, 0, 0, 0.00439, Angle.HOURANGLE);
		mud = new Angle(minus, 0, 0, .1897, Angle.ARCDEGREES);
		//mud = new Angle(minus, 0, 0, 18.97, Angle.ARCDEGREES);

		// Besselian day number table values
		tau = 0.4211;
		A = new Angle(plus, 0, 0, 4.086, Angle.ARCDEGREES);
		B = new Angle(plus, 0, 0, 7.945, Angle.ARCDEGREES);
		C = new Angle(minus, 0, 0, 5.629, Angle.ARCDEGREES);
		D = new Angle(minus, 0, 0, 19.551, Angle.ARCDEGREES);
		E = new Angle(minus, 0, 0, 0.0016, Angle.HOURANGLE);


		Position pos = StarAlg.mean_to_apparent_Besselian(m, n, alpha0, delta0, tau, mua, mud, A, B, C, D, E, true);
		System.out.println("---------------------------------------------");
		alpha0.println("alpha0  = ", Angle.SHA);
		pos.RA.println("alpha  = ", Angle.SHA);
		delta0.println("delta0  = ", Angle.ARCDEGREES);
		pos.dec.println("delta  = ", Angle.ARCDEGREES);
	}
}
