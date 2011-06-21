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

public class AnnualPrecession
{
	public static void main(String[] args)
	{
		boolean plus = true, minus = false;
		Angle m1971, n1971, m1974, n1974, m1973, n1973, m1980, n1980;
		Angle alpha0, delta0;

		// m and n
		m1971 = new Angle(true, 0, 0, 3.07366, Angle.HOURANGLE);
		m1973 = new Angle(true, 0, 0, 3.07370, Angle.HOURANGLE);
		m1974 = new Angle(true, 0, 0, 3.07372, Angle.HOURANGLE);
		m1980 = new Angle(true, 0, 0, 3.07383, Angle.HOURANGLE);
		n1971 = new Angle(true, 0, 0, 1.33605, Angle.HOURANGLE);
		n1973 = new Angle(true, 0, 0, 1.33604, Angle.HOURANGLE);
		n1974 = new Angle(true, 0, 0, 1.33603, Angle.HOURANGLE);
		n1980 = new Angle(true, 0, 0, 1.33600, Angle.HOURANGLE);
		//m.println("Annual precession in R.A.");
		//n.println("Annual precession in Dec.");

		// Star Theta Oct 1971
		System.out.println("Star Theta Oct 1971 -> 1972");
		alpha0 = new Angle(plus, 0, 0, 7.7, Angle.HOURANGLE);
		delta0 = new Angle(minus, 77, 13, 32., Angle.ARCDEGREES);
		StarAlg.annual_precession(m1973, n1973, alpha0, delta0, true);
		System.out.println();

		// Star Theta Oct 1973
		System.out.println("Star Theta Oct 1973 -> 1974");
		alpha0 = new Angle(plus, 0, 0, 13.8, Angle.HOURANGLE);
		delta0 = new Angle(minus, 77, 12, 53., Angle.ARCDEGREES);
		StarAlg.annual_precession(m1973, n1973, alpha0, delta0, true);
		System.out.println();

		// Star 33 PsC 1973
		System.out.println("Star 33 PsC 1973 -> 1974");
		alpha0 = new Angle(plus, 0, 3, 57.1, Angle.HOURANGLE);
		delta0 = new Angle(minus, 5, 51, 31., Angle.ARCDEGREES);
		StarAlg.annual_precession(m1973, n1973, alpha0, delta0, true);
		System.out.println();

		// Star 33 PsC 1974
		System.out.println("Star 33 PsC 1974 -> 1975");
		alpha0 = new Angle(plus, 0, 4, 0.2, Angle.HOURANGLE);
		delta0 = new Angle(minus, 5, 51, 11., Angle.ARCDEGREES);
		StarAlg.annual_precession(m1974, n1974, alpha0, delta0, true);
		System.out.println();

		// Star 33 PsC 1980
		System.out.println("Star 33 PsC 1980 -> 1981");
		alpha0 = new Angle(plus, 0, 4, 18.6, Angle.HOURANGLE);
		delta0 = new Angle(minus, 5, 49, 10., Angle.ARCDEGREES);
		StarAlg.annual_precession(m1980, n1980, alpha0, delta0, true);
	}

}
