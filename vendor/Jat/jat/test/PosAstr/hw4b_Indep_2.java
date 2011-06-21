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

public class hw4b_Indep_2
{
	public static void main(String[] args)
	{
		boolean plus=true,minus=false;
		Angle alpha0, delta0;
		double tau;
		Angle mua, mud; // mu alpha, mu delta
		Angle f, g, G, h, H, i; // tabulated day numbers

		// Alpha Tauri
		System.out.println("Aldebaran June 3, 1980 Independent Day Numbers");
		System.out.println("----------------------------------------------");
		alpha0 = new Angle(true, 4, 34, 46., Angle.HOURANGLE);
		delta0 = new Angle(true, 16, 28, 12., Angle.ARCDEGREES);
		mua = new Angle(true, 0, 0, 0.00439, Angle.HOURANGLE);
		mud = new Angle(false, 0, 0, 0.1897, Angle.ARCDEGREES);

		// Independent day number table values
		tau = 0.4211;
		f = new Angle(plus, 0, 0, 0.6251, Angle.HOURANGLE);
		g = new Angle(plus, 0, 0, 8.934, Angle.ARCDEGREES);
		G = new Angle(plus, 4,11,8., Angle.HOURANGLE);
		h = new Angle(plus, 0,0,20.345, Angle.ARCDEGREES);
		H = new Angle(plus, 13,4,15, Angle.HOURANGLE);
		i = new Angle(minus, 0, 0, 2.441, Angle.ARCDEGREES);

		Position pos = StarAlg.mean_to_apparent_Independent(alpha0, delta0, tau, mua, mud, f, g, G, h, H, i, true);
		System.out.println("---------------------------------------------");
		alpha0.println("alpha0  = ", Angle.SHA);
		pos.RA.println("alpha  = ", Angle.SHA);
		delta0.println("delta0  = ", Angle.ARCDEGREES);
		pos.dec.println("delta  = ", Angle.ARCDEGREES);
	}
}
