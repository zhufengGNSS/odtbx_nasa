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

package jat.forces.gravity;
import jat.eph.*;
import jat.matvec.data.*;
import jat.spacetime.*;

/**
 * @author dgaylor
 *
 */

public class ThirdBodyGravity {
	
	private double GM;
	private DE405 de405 = null;
	private DE405_Body centralBody = null;
	private DE405_Body thirdBody = null;

	/**
	 * Create the third body point mass gravity model
	 * @param de DE405 ephemeris object
	 * @param central DE405_Body of the Central Body
	 * @param third DE405_Body of the Third Body
	 * @param mu double containing the value of GM (determines the units)
	 */	
	
	public ThirdBodyGravity(DE405 de, DE405_Body central, DE405_Body third, double mu){
		de405 = de;
		centralBody = central;
		thirdBody = third;
		GM = mu;
	}
	/**
	 * Compute the perturbational acceleration due to a point mass
	 * @param t Time
	 * @param r satellite position vector wrt Central Body
	 * @return VectorN containing the acceleration
	 */	
	public VectorN acceleration(Time t, VectorN r){
		// get central body position vector in meters
		VectorN r_cb = de405.get_planet_pos(centralBody, t);
		// get third body position vector in meters
		VectorN r_tb = de405.get_planet_pos(thirdBody, t);
		// get third body position wrt central body
		VectorN s = r_tb.minus(r_cb);
		// get position of satellite wrt third body
		VectorN d = s.minus(r);
		// compute acceleration
		double dmag = d.mag();
		double denom1 = 1.0/(dmag*dmag*dmag);
		VectorN term1 = d.times(denom1);
		double smag = s.mag();
		double denom2 = 1.0/(smag*smag*smag);
		VectorN term2 = s.times(denom2);
		VectorN accel = (term1.minus(term2)).times(GM);
		return accel;
	}

}
