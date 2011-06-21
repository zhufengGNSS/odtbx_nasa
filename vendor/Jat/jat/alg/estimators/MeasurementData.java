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
 
 package jat.alg.estimators;


/**
* The MeasurementData.java Class provides an interface to get
* the measurement data to the extended Kalman filter
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public interface MeasurementData {
	/** Returns the measurement
	 * @param index int containing the measurement index
	 * @return double containing the measurement corresponding to the index.
	 */
	public double z (int index);

	/**
	 * Returns the time of the measurement
	 * @param index int containing the measurement index
	 * @return double containing time of the measurement corresponding to the index.
	 */
	public double time (int index);

	/**
	 * Returns the predicted measurement based on the current state
	 * @param index measurement index
	 * @param t time of the measurement
	 * @param xref VectorN with the current state at the measurement time
	 */

	public boolean hasNext(int index);
	

}
