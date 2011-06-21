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
 
package jat.ins;
import jat.matvec.data.*;
import java.io.*;


/**
* The INS_Measurement.java Class represents a single INS measurement.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class INS_Measurement implements Serializable {
	
	/** Epoch time of measurement in sim time (seconds) */
	public double t;
		
	/** Accelerometer measurement */
	public VectorN f;
	
	/** Gyro measurement */
	public VectorN omega;
	
	/** Constructor
	 * @param tsim Time of the measurement in sim time (sec).
	 * @param sf Accelerometer measurement (m/s)
	 * @param w Gyro measurement (rad/s)
	 */
	public INS_Measurement(double tsim, VectorN sf, VectorN w) {
		this.t = tsim;
		this.f = sf;
		this.omega = w;
	}
	

}
