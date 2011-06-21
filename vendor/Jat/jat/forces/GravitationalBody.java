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
package jat.forces;

import jat.cm.Constants;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.Time;

/**
 * The GravitationalBody class serves to model a point mass gravitational
 * force from a generic third body.
 * 
 * @author Richard C. Page III
 *
 */
public class GravitationalBody implements ForceModel {

	private static final long serialVersionUID = 1L;
	
	/** Gravitational constant for the body. [m^3/s^2]
	 */
	protected double mu;
	
	/** ECI position of the body. [m]
	 */
	protected VectorN r_body;
	
	/** ECI velocity of the body. [m/s]
	 */
	protected VectorN v_body;
	
	/** Default Constructor
	 */
	public GravitationalBody(){
	    mu = Constants.GM_Earth;
	}
	/** Constructor creating a gravitational body.
	 * @param grav_mass The Gravitational constant of the body. [m^3/s^2]
	 * @param r The ECI position vector of the body. [m]
	 * @param v The ECI velocity vector of the body. [m/s]
	 */
	public GravitationalBody(double grav_mass, VectorN r, VectorN v){
		mu = grav_mass;
		r_body = r;
		v_body = v;
	}
	/**
	 * Constructor: assumes the gravitational body is stationary at the origin.
	 * @param grav_mass The Gravitational constant of the body. [m^3/s^2]
	 */
	public GravitationalBody(double grav_mass){
	    mu = grav_mass;
	    r_body = new VectorN(0,0,0);
	    v_body = new VectorN(0,0,0);
	}
	

	/**
	 * Updates the position of the gravitational body with respect 
	 * to the inertial frame. [m]
	 * @param r
	 */
	public void updatePosition(VectorN r){
	    r_body = r;
	}
		
	/**
	 * Updates the velocity of the gravitational body with respect
	 * to the inertial frame. [m/s]
	 * @param v
	 */
	public void updateVelocity(VectorN v){
	    v_body = v;
	}
	
	/**
	 * Returns the inertial position of the gravitational body.
	 * @return position vector [m]
	 */
	public VectorN getPosition(){
		try{
			return r_body;
		}catch(NullPointerException ne){
			System.err.println("NullPointerException - call to getPosition when not stored");
			return null;
		}
	}
	/**
	 * Returns the inertial position of the gravitational body.
	 * @return velocity vector [m/s]
	 */
	public VectorN getVelocity(){ 
		try{
			return v_body;
		}catch(NullPointerException ne){
			System.err.println("NullPointerException - call to getVelocity when not stored");
			return null;
		}
	}
	
    /** Compute the acceleration due to a point mass gravity field (usually due to a third body).
     * @param r ECI position vector about the body [m].
     * @return Acceleration due to a point mass in m/s^2.
     */
    public VectorN accelPointMass(VectorN r) 
    {
        //  Relative position vector of satellite w.r.t. point mass
    	
        double rmag = r.mag();
        double rcubed = rmag * rmag *rmag;

        VectorN temp1 = r.divide(rcubed);

        VectorN out = temp1.times(-mu);

        // Acceleration

        return  out;

    }
    /** Call the relevant methods to compute the acceleration.
	 * @param t Time reference class
     * @param bRef Body reference class
     * @param sc Spacecraft parameters and state
     * @return the acceleration [m/s^s]
     */
    public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc){
        return accelPointMass(sc.r());
    }
	
}
