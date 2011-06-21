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
package jat.spacetime;

//import jat.spacetime.*;
import java.io.Serializable;
import jat.forces.*;
import jat.spacecraft.Spacecraft;
//import jat.math.*;

/**
 * This class encapsulates the universe model for a simulation.  This includes the
 * following elements: Time, Reference Frame, IERS data, and Force Models.
 * 
 * @author Richard C. Page III
 *
 */
public class UniverseModel implements Serializable {

	private static final long serialVersionUID = 1L;

    /**
     * Reference Frame information, including conversions and relevant parameters 
     * for the body at the center of the frame.
     */
    //** Note: This data member could be implemented better.  In the future, it should
    //** be generalized to represent any central body.  To do this, it should be changed
    //** either to be of type 'jat.spacetime.BodyRef' or include several different
    //** reference frames with a class variable indicating which is the central body. 
    public EarthRef earthRef;
    /**
     * Time information, including conversions between time systems.
     */
    public Time time;
    /**
     * Flag indicating whether or not to use the IERS Bulletin A correction data.
     */
    protected boolean use_iers = false;
    /**
     * International Earth Rotation Service correction data.
     */
    protected FitIERS iers;
    /**
     * List of the force models present in the simulation universe.
     */
    protected ForceModelList forces;
    
    /**
     * Default Constructor.  Initializes to epoch of J2000.
     *
     */
    public UniverseModel(){
        time = new Time(TimeUtils.MJD_J2000);
        earthRef = new EarthRef(time.mjd_ut1(),time.mjd_tt());
        forces = new ForceModelList();
    }
    
    /**
     * Constructor.  Initializes to the given UTC epoch.
     * @param MJD_UTC Modified Julian Date in Universal Coordinated Time
     */
    public UniverseModel(double MJD_UTC){
        time = new Time(MJD_UTC);
        earthRef = new EarthRef(time.mjd_ut1(),time.mjd_tt());
        forces = new ForceModelList();
    }
    
    /**
     * Set the current time to the given value.  Updates time and bodyRef
     * @param mjd_utc Modified Julian Date in Universal Coordinated Time
     */
    public void set_time(double mjd_utc){
        time = new Time(mjd_utc);
        boolean temp1 = earthRef.use_moon;
        boolean temp2 = earthRef.use_sun;
        earthRef = new EarthRef(time.mjd_ut1(),time.mjd_tt());
        earthRef.set_use_moon(temp1);
        earthRef.set_use_sun(temp2);
    }
    
    /**
     * Sets whether the simulation applies the IERS pole wander and UT1-UTC corrections.
     * It then updates the simulation models accordingly.
     * 
     * @param b true to use IERS false to ignore it.
     */
    public void set_use_iers(boolean b){
        use_iers = b;
        if(iers == null) iers = new FitIERS();
        update(time.get_sim_time());
    }
    
    /**
     * Sets whether the simulation needs to calculate the sun's position.  Updates
     * the relevant objects accordingly.
     * @param b Boolean flag.
     */
    public void set_compute_sun(boolean b){
        earthRef.set_use_sun(b);
        update(time.get_sim_time());
    }
    
    /**
     * Sets whether the simulation needs to calculate the moon's position.  Updates
     * the relevant objects accordingly.
     * @param b Boolean flag.
     */
    public void set_compute_moon(boolean b){
        earthRef.set_use_moon(b);
        update(time.get_sim_time());
    }
    
    /**
     * Return the current time in Modified Julian Date of Universal Coordinated Time
     * @return The current time.
     */
    public double get_mjd_utc(){
        return time.mjd_utc();
    }
    
    /**
     * Add a force model to the universe.
     * @see jat.forces.ForceModel
     * @param f The force model.
     */
    public void addForce(ForceModel f){
        forces.addForce(f);
    }
    
    /**
     * Initializes the JPL DE405 ephemerides for use.
     *
     */
    public void initializeMoonEphem(){
    	earthRef.initializeMoonEphem(time.mjd_tt());
    }
    
    /**
     * Initializes the JPL DE405 ephemerides for use.
     *
     */
    public void initializeSunEphem(){
    	earthRef.initializeSunEphem(time.mjd_tt());
    }

    /**
     * Update the models to 't' seconds since the start epoch.
     * @param t [sec]
     */
    public void update(double t){
        if(use_iers){// && MathUtils.mod(t,60)==0){           
        	try{
		    double[] param = iers.search(time.mjd_tt());
		    earthRef.setIERS(param[0],param[1]);
		    time.set_UT1_UTC(param[2]);
        	}catch(Exception e){
        		//e.printStackTrace();
        		System.err.println("Error: Unable to interpret IERS Polar Motion data.");
        		System.err.println("       Polar motion data will be set to zero.");
        		use_iers = false;
        	}
		}else{
			int donothing = 0;
		}
        time.update(t);
        //if(MathUtils.mod(t,60)==0)
            earthRef.update(time.mjd_ut1(),time.mjd_tt());
    }
    
    /**
     * Though UniverseModel does not implement the Derivatives interface, this
     * method has a similar purpose in providing the derivatives of a spacecraft
     * state.  Given the state and relevant properties for a spacecraft, derivs
     * calculates the acceleration due to the current force list and returns
     * the derivaties.
     * 
     * This method is called by SimModel to provide the derivatives to the 
     * RungeKutta integrator.
     * 
     * @see jat.sim.SimModel
     * 
     * @param t Current simulation time in seconds since the initial epoch.
     * @param sc Spacecraft to opperate upon.
     * @return The derivatives of the spacecraft state.
     */
    public double[] derivs(double t, Spacecraft sc) {
        double[] X = sc.toStateVector();
        double[] out = new double[X.length];
        time.update(t); //DMS this may not be necessary since done earlier in update; just doing this here does not update earthRef
        double[] accel = forces.acceleration(time,earthRef,sc);  // DMS this hardcodes earth-based computation
        out[0] = X[3];
        out[1] = X[4];
        out[2] = X[5];
        out[3] = accel[0];
        out[4] = accel[1];
        out[5] = accel[2];
        for(int i=6; i<X.length; i++) out[i] = 0;
        return out;
    }

	/**
	 * Return the given force (in order as they were added)
	 * @param i integer index
	 * @return ForceModel (to be cast to the relevant force object)
	 */
	public ForceModel getForce(int i){
	    return forces.get_Force(i);
	}
	
	/**
	 * Return the number of forces
	 * @return int
	 */
	public int getForceSize(){
	    return forces.getForceSize();
	}

}
