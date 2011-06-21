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

import jat.matvec.data.VectorN;
import jat.spacetime.*;

import java.io.Serializable;
import java.util.*;
import jat.spacecraft.Spacecraft;

/**
 * ForceModelList is a countainer for force models.  It allows the simulation to encapsulate
 * a number of different force models and obtain their combined effect on a spacecraft. It 
 * also allows some flags to tell the other models whether or not to use the JPL and IERS
 * algorithms.
 * 
 * @author Richard C. Page III
 *
 */
public class ForceModelList implements Serializable {

	private static final long serialVersionUID = 1L;
	
	/** Sum of all accelerations [m/s^2] 
	 */
	protected VectorN acceleration;
	/** List of forces
	 */
	protected ArrayList forces;
	protected boolean use_iers = true;
	protected boolean use_sun = true;
	protected boolean use_moon = true;
	/**
	 * Tells the integrator whether to update within each timestep.
	 */
	protected boolean use_update = true;
	
	
	/** Default Constructor
	 * 
	 */
	public ForceModelList(){
		acceleration = new VectorN(0,0,0);
		forces = new ArrayList();
	}
	
        
	/**
	 * Adds a generic force to the list
	 * @param f Object which implements the ForceModel interface
	 */
	public void addForce(ForceModel f){
		forces.add(f);
	}
	

	/** Compute the acceleration
	 * @return acceleration [m/s^2]
     */
	public VectorN computeAccel(Time t, BodyRef ref, Spacecraft sc){
		//retrieve/update acceleration
		acceleration = new VectorN(3);
		ForceModel temp;
		VectorN accel;
		for(int i=0; i<forces.size(); i++){
			temp = (ForceModel)forces.get(i);
			accel = temp.acceleration(t,ref,sc);
			acceleration = acceleration.plus(accel);
		}
		return acceleration;
	}

	/**
	 * Return the given force (in order as they were added)
	 * @param i integer index
	 * @return ForceModel (to be cast to the relevant force object)
	 */
	public ForceModel get_Force(int i){
	    return (ForceModel)forces.get(i);
	}
	/**
	 * Return the number of forces
	 * @return Number of forces
	 */
	public int getForceSize(){
	    return forces.size();
	}
	/**
	 * Set whether to update within each timestep;
	 * @param b boolean flag
	 */
	public void set_use_update(boolean b){
	    this.use_update = b;
	}
	/**
	 * Set whether to compute the JPL Sun position
	 * @param b boolean flag
	 */
	public void set_use_sun(boolean b){
	    this.use_sun = b;
	}
	/**
	 * Set whether to compute the JPL Moon position
	 * @param b boolean flag
	 */
	public void set_use_moon(boolean b){
	    this.use_moon = b;
	}
	
	/** Compute Acceleration [m/s^2]
	 * @param t Time
	 * @param bRef Reference Frame
	 * @param sc Spacecraft parameters and state
	 */
	public double[] acceleration(Time t, BodyRef bRef, Spacecraft sc){
		computeAccel(t,bRef,sc);
		return acceleration.x;
	}
}
