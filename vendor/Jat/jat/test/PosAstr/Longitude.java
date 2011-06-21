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

import jat.matvec.data.*;
import jat.cm.*;

public class Longitude
{
	// Austin 
	// Latitude   30 o 16' ''
	// Longitude  97 o 45' ''
	// Longitude h m s

	public static void main(String[] args)
	{
		boolean plus = true, minus = false;
		int i;
		//double[] t_array = {1,2,3,4};
		//double[] y_array = {1,4,9,16 };

		System.out.println("Latitude and Longitude measurement");
		System.out.println("March 7, 2004 in Austin, Texas");
		System.out.println("Table values:");
		System.out.println("Latitude   30 o 16'");
		System.out.println("Longitude  97 o 45'");
		System.out.println();

		// Find maximum of observation data
		VectorN x=parabola_fit(ObsData.t_C,ObsData.y_C);
		//VectorN x=parabola_fit(ObsData.t_C,ObsData.y_C);
		x.print("Parabola fit");			
		double max_t = -x.get(1) / (2. * x.get(2));
		double max_angle = x.get(0) + x.get(1) * max_t + x.get(2) * max_t * max_t;
		//System.out.println("max t     = " + max_t);
		//System.out.println("max angle = " + max_angle);

		// Corrections for errors
		double clockerr = -28. / 3600;
		double eq_of_time = -11. / 60. + 17. / 3600.;
		double t_corr = max_t + clockerr + eq_of_time;
		double angleerr = -8.8 / 60;
		double angle_corr = max_angle + angleerr;
		//System.out.println("corr t    = " + t_corr);
		//System.out.println("corr angle = " + angle_corr);
		System.out.println();

		// Latitude
		double decl = - (5. + 30. / 60.);
		double lat = 90. - (angle_corr / 2. - decl);
		new Angle(lat,Angle.DEGREES).println("Latitude",Angle.ARCDEGREES);
		//System.out.println("Latitude = " + lat);
		System.out.println();

		// Longitude
		Angle LAT = new Angle(max_t, Angle.DECIMALHOURS);
		Angle LMT = new Angle(t_corr, Angle.DECIMALHOURS);
		Angle GMTcorr = new Angle(plus, 6, 0, 0., Angle.HOURANGLE);
		Angle GMT = new Angle(LMT.radians + GMTcorr.radians, Angle.RADIANS);
		Angle noon = new Angle(plus, 12, 00, 0., Angle.HOURANGLE);
		Angle longitude = new Angle(GMT.radians - noon.radians, Angle.RADIANS);
		LAT.println("LAT of transit", Angle.HOURANGLE);
		LMT.println("LMT of transit", Angle.HOURANGLE);
		GMT.println("GMT of transit", Angle.HOURANGLE);
		System.out.println();
		longitude.println("Longitude", Angle.HOURANGLE);
		longitude.println("Longitude", Angle.ARCDEGREES);
	}
	

	static VectorN parabola_fit(double[] t_array, double[] y_array)
	{
		int i;
		VectorN y = new VectorN(y_array);
		VectorN t = new VectorN(t_array);
		int n = t_array.length;
		
		Matrix H = new Matrix(n, 3);
		for (i = 0; i < n; i++)
		{
			H.set(i, 0, 1.);
			H.set(i, 1, t.get(i));
			H.set(i, 2, t.get(i) * t.get(i));
		}

		Matrix HT = H.transpose();
		Matrix HTH = HT.times(H);
		Matrix HTHinv = HTH.inverse();
		VectorN Hy = HT.times(y);
		//		H.print("H");
		//		HT.print("HT");
		//		HTH.print("HTH");
		//		HTHinv.print("HTHinv");
		VectorN x = HTHinv.times(Hy);
		
		return x;		
	}
	
}
