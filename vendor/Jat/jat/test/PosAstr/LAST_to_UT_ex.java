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

public class LAST_to_UT_ex
{
	public static void main(String[] args)
	{
		boolean plus = true, minus = false;
		Angle LAST;			// observed local apparent sidereal time		
		Angle EoE, Long;
		Angle GMST_0hUT; 	// Greenwich mean sidereal time at 0h UT on given date 
		Angle TableVIII; 	// Mean sidereal to mean solar time on given date

		// Book example page 85: 1960 March 7 
		System.out.println("Book example page 85: 1960 March 7"); 
		LAST = new Angle(plus, 13, 5, 37.249, Angle.HOURANGLE);
		EoE = new Angle(plus, 0, 0, 0.046, Angle.HOURANGLE);
		Long = new Angle(plus, 5, 8, 15.75, Angle.HOURANGLE);
		GMST_0hUT = new Angle(plus, 10, 58, 50.971, Angle.HOURANGLE);
		TableVIII = new Angle(plus, 0, 1, 11.269, Angle.HOURANGLE);
		LAST_to_UT(LAST, EoE, Long,GMST_0hUT,TableVIII);
		System.out.println(""); 

		// Mid Term: 1980 July 8 
		System.out.println("Mid Term: 1980 July 8"); 
		LAST = new Angle(plus, 23, 15, 57.881, Angle.HOURANGLE);
		EoE = new Angle(plus, 0, 0, 0.631, Angle.HOURANGLE);
		Long = new Angle(plus, 72, 39, 33.0, Angle.ARCDEGREES);
		GMST_0hUT = new Angle(plus, 19, 04, 24.215, Angle.HOURANGLE);
		TableVIII = new Angle(plus, 0, 1, 28.828, Angle.HOURANGLE);
		LAST_to_UT(LAST, EoE, Long,GMST_0hUT,TableVIII);
	}

	static void LAST_to_UT(Angle LAST, Angle EoE, Angle Long, Angle GMST_0hUT, Angle TableVIII)
	{
		Angle LMST, GMST,MST_int,UT;

		LMST = new Angle(LAST.radians + EoE.radians, Angle.RADIANS);
		GMST = new Angle(LMST.radians + Long.radians, Angle.RADIANS);
		MST_int= new Angle(GMST.radians-GMST_0hUT.radians , Angle.RADIANS);
		UT= new Angle(MST_int.radians-TableVIII.radians , Angle.RADIANS);
		
		LAST.println("Local apparent sidereal time", Angle.HOURANGLE);
		EoE.println("Equation of the Equinoxes", Angle.HOURANGLE);
		LMST.println("Local mean sidereal time", Angle.HOURANGLE);
		Long.println("Longitude", Angle.HOURANGLE);
		GMST.println("Greenwich mean sidereal time", Angle.HOURANGLE);
		GMST_0hUT.println("Greenwich mean sidereal time at 0h UT", Angle.HOURANGLE);
		MST_int.println("Greenwich mean sidereal time interval", Angle.HOURANGLE);
		TableVIII.println("Mean sidereal time to mean solar time", Angle.HOURANGLE);
		UT.println("U.T. of observation", Angle.HOURANGLE);

	}
}
