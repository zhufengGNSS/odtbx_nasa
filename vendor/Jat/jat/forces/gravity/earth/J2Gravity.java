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
 * File Created on Aug 30, 2003
 */
 
package jat.forces.gravity.earth;
import jat.cm.Constants;
import jat.forces.ForceModel;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.Time;

/**
* The J2Gravity.java Class computes the gravitational acceleration and
* gravity gradient matrix for a J2 gravity field. Values are from JGM-3.
* 
* Reference: Montenbruck
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/ 
public class J2Gravity implements ForceModel {

	private static final long serialVersionUID = 1L;

	private final static double j2 = 0.00108263;
    /** GM in m^3/s^2.
     */
    private static final double mu = Constants.GM_Earth;  // [m^3/s^2]; JGM3
    /** Radius of planet in m.
     */
    private static final double re = Constants.rEarth_STKJGM3;  // Radius Earth [m]; JGM3;
    
    private VectorN grav;
    
    private Matrix gradient;
    
    public J2Gravity(){}
    
    /**
     * Contructor
     * @param ECI position vector
     */
    public J2Gravity(VectorN r){
    	this.compute(r);
    }
    
    /**
     * Returns the local gravity vector
     * @return VectorN containing the local gravity vector
     */
    public VectorN local_gravity() {
    	return this.grav;
    }
    
    /**
     * Returns the gravity gradient matrix
     * @return Matrix containing the gravity gradient
     */
    public Matrix gravityGradient() {
    	return this.gradient;
    }
	
	private void compute (VectorN r){
		// strip out components of r
		double xx = r.x[0];
		double yy = r.x[1];
		double zz = r.x[2];
		
		// compute some intermediate values
		double rmag = r.mag();
		double rcubed = rmag * rmag * rmag;
		double rsq = rmag * rmag;
		double re_r = re / rmag;
		re_r = re_r * re_r;
		double zsq_rsq = (5.0 * zz * zz / rsq) - 1.0;		
		
		// compute accelerations due to gravity
		double ll = -1.0 * (mu * xx / rcubed) * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
		double mm = -1.0 * (mu * yy / rcubed) * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
		double nn = -1.0 * (mu * zz / rcubed) * (1.0 - 1.5 * re_r * j2 * (zsq_rsq - 2.0));
		this.grav = new VectorN(ll, mm, nn);
		
		// compute some intermediate values
		double r5 = rsq * rcubed;
		double mur5 = mu / r5;
		double mur3 = mu / rcubed;
		double sz2r2 = 7.0 * zz * zz / rsq;
		double muxyr5 = mu * xx * yy / r5;
		double muxzr5 = mu * xx * zz / r5;
		double muyzr5 = mu * yy * zz / r5;
		double bracket1 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 1.0);
		double bracket3 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 3.0);
		double bracket5 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 5.0);

		// partials of ll wrt r
		double dldx = ll / xx + mur5 * xx * xx * bracket1;
		double dldy = muxyr5 * bracket1;
		double dldz = muxzr5 * bracket3;
		
		// partials of mm wrt r
		double dmdx = dldy;
		double dmdy = mm / yy + mur5 * yy * yy * bracket1;
		double dmdz = muyzr5 * bracket3;

		// partials of nn wrt r
		double dndx = muxzr5 * bracket3;
		double dndy = dmdz;
		double dndz = nn / zz + mur5 * zz * zz * bracket5;
		
		// set the gravity gradient matrix
		this.gradient = new Matrix(3, 3);
		this.gradient.set(0, 0, dldx);
		this.gradient.set(0, 1, dldy);
		this.gradient.set(0, 2, dldz);
		this.gradient.set(1, 0, dmdx);
		this.gradient.set(1, 1, dmdy);
		this.gradient.set(1, 2, dmdz);
		this.gradient.set(2, 0, dndx);
		this.gradient.set(2, 1, dndy);
		this.gradient.set(2, 2, dndz);		
	}
	
    /** Acceleration, from the ForceModel interface.  This calculation
     * does not depend on the Time or BodyRef argument, only the Spacecraft
     * state.
	 * @param t Time reference class
     * @param bRef Body reference class
     * @param sc Spacecraft parameters and state
     * @return the acceleration [m/s^s]
     */
    public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc) {
    	compute(sc.r());
    	return local_gravity();
    }

}
