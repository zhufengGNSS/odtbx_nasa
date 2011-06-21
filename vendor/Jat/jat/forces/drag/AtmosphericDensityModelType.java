/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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

package jat.forces.drag;

/**
 * @author David Gaylor
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtmosphericDensityModelType {
	
	private String name;
	  private AtmosphericDensityModelType (String nm) { name = nm; }

	  public String toString() { return name; }


	  public final static AtmosphericDensityModelType
	    CIRA_Exponential = new AtmosphericDensityModelType("CIRA_Exponential"),
	    Exponential = new AtmosphericDensityModelType("Exponential"),
	    HarrisPriester = new AtmosphericDensityModelType("HarrisPriester"),
	    NRLMSISE = new AtmosphericDensityModelType("NRLMSISE");

	  public final static AtmosphericDensityModelType[] index =  {
		  CIRA_Exponential, Exponential, HarrisPriester, NRLMSISE
	  };
	  
	  public final static int num_models = 4;
}
