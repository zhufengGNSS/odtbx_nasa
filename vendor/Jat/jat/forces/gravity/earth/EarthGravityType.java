/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2003 United States Government as represented by the
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
 */

package jat.forces.gravity.earth;

import jat.forces.gravity.GravityModelType;

public class EarthGravityType extends GravityModelType {

	private EarthGravityType(String nm) {
		super(nm);
	}
	  public final static EarthGravityType
	    EGM96 = new EarthGravityType("EGM96"),
	    GEMT1 = new EarthGravityType("GEMT1"),
	    JGM2 = new EarthGravityType("JGM2"),
	    JGM3 = new EarthGravityType("JGM3"),
	    WGS84 = new EarthGravityType("WGS84"),
	    WGS84_EGM96 = new EarthGravityType("WGS84_EGM96");
	  public final static EarthGravityType[] index =  {
	  		EGM96, GEMT1, JGM2, JGM3, WGS84, WGS84_EGM96
	  };
	  
	  public final static int num_models = 8;	

}
