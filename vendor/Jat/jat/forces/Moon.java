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

//import jat.cm.Constants;
import jat.cm.Constants;
import jat.eph.DE405;
import jat.eph.DE405_Body;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;

/**
 * Simple class to model the gravitational attraction of the moon.
 * 
 * @author Richard C. Page III
 *
 */
public class Moon extends GravitationalBody {

	private static final long serialVersionUID = 1L;

    protected DE405 jpl_ephemeris;
    
    /**
     * Default constructor
     */
    public Moon() {
        super();
        this.mu = Constants.GM_Moon; 
        jpl_ephemeris = new DE405();
    }
    
    /**
     * Allows the initialization of the Moon vector using DE405 ephemerides
     * @param mjd_utc Modified Julian Date in Universal Coordinated Time
     */
    public Moon(double mjd_utc){
    	super();
        this.mu = Constants.GM_Moon; 
        jpl_ephemeris = new DE405();
        double mjd_tt = TimeUtils.UTCtoTT(mjd_utc);
        VectorN temp = jpl_ephemeris.get_planet_posvel(DE405_Body.GEOCENTRIC_MOON, mjd_tt);
        this.r_body = temp.get(0, 3);
        this.v_body = temp.get(3, 3);
    }
    
    /**
     * For use with matlab
     *
     */
    public Moon(boolean usematlab) {
        super();
        this.mu = Constants.GM_Moon;      
        jpl_ephemeris = new DE405();
    }

    /**
     * Constructor
     * @param grav_mass [m^3/s^2]
     * @param r position [m]
     * @param v velocity [m/s]
     */
    public Moon(double grav_mass, VectorN r, VectorN v) {
        super(grav_mass, r, v);        
        jpl_ephemeris = new DE405();
    }
    /**
     * Compute the acceleration.
     * @param eRef Earth Reference
     * @param sc Spacecraft parameters and state
     * @return Vector 3 of acceleration [m/s^2]
     */
    public VectorN acceleration(jat.spacetime.EarthRef eRef, Spacecraft sc) {

      //ReferenceFrameTranslater xlater = 
        //new ReferenceFrameTranslater(eRef, moonRef, new Time(eRef.mjd_utc()));

      // Get a vector (in the passed in reference frame) to the
      // spacecraft and to the moon.
      VectorN r_moon = eRef.get_JPL_Moon_Vector(); // m
      VectorN r = sc.r();
      VectorN d = r.minus(r_moon);
      
      
      // Divide the moon vector and the moon-to-spacecraft vector
      // by their magnitude cubed, add them, and scale by mu.
      double dmag = d.mag();
      double dcubed = dmag * dmag *dmag;
      VectorN temp1 = d.divide(dcubed);
      double smag = r_moon.mag();
      double scubed = smag * smag * smag;
      VectorN temp2 = r_moon.divide(scubed);
      VectorN sum = temp1.plus(temp2);
      VectorN accel = sum.times(-mu);

      // We keep track of the position of the moon at the last
      // known time (in the last known reference frame)
      this.r_body = r_moon; // m
        
      return  accel;
    }
    /** Call the relevant methods to compute the acceleration.
	 * @param t Time reference class
     * @param bRef Body reference class
     * @param sc Spacecraft parameters and state
     * @return the acceleration [m/s^s]
     */
    public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc){
        //ReferenceFrameTranslater xlater = 
          //new ReferenceFrameTranslater(bRef, moonRef, t);
        
        // Get a vector (in the passed in reference frame) to the
        // spacecraft and to the moon.        
        VectorN r_moon = new VectorN(jpl_ephemeris.get_planet_pos(DE405_Body.GEOCENTRIC_MOON, t.mjd_tt()));
        //VectorN r_moon = xlater.translatePointBack(new VectorN(3));
        VectorN r = sc.r();
        VectorN d = r.minus(r_moon);
        
        // Divide the moon vector and the moon-to-spacecraft vector
        // by their magnitude cubed, add them, and scale by mu.
        double dmag = d.mag();
        double dcubed = dmag * dmag *dmag;
        VectorN temp1 = d.divide(dcubed);
        double smag = r_moon.mag();
        double scubed = smag * smag * smag;
        VectorN temp2 = r_moon.divide(scubed);
        VectorN sum = temp1.plus(temp2);
        VectorN accel = sum.times(-mu);

        // We keep track of the position of the moon at the last
        // known time (in the last known reference frame) - but we
        // keep it in km.
        this.r_body = r_moon.divide(1000);
        
        return  accel;
    }

      
   
    
}
