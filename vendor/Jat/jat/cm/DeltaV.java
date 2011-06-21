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
 *
 * 
 * File Created on May 19, 2003
 */
 
package jat.cm;
import jat.matvec.data.*;
import java.io.*;


/**
* The DeltaV.java Class represents a single Delta-V.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class DeltaV implements Serializable {
	
	/** Epoch time of delta-v in sim time (seconds) */
	public double t;
		
	/** Delta-v vector */
	public VectorN dv;
	
	
	/** Constructor
	 * @param tsim Time of the measurement in sim time (sec).
	 * @param deltav Delta-V vector (m/s)
	 */
	public DeltaV(double tsim, VectorN deltav) {
		this.t = tsim;
		this.dv = deltav;
	}
	

}
