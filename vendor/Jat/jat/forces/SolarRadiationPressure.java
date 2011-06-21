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
import jat.eph.DE405_Body;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyCenteredInertialRef;
import jat.spacetime.BodyRef;
import jat.spacetime.ReferenceFrame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;

/**
 * <P>
 * The SolarRadiationPressure class provides a generic model for the influence
 * of solar radiation pressure on a satellite.
 * 
 * @author Richard C. Page III
 *
 */
public class SolarRadiationPressure implements ForceModel {

	private static final long serialVersionUID = 1L;

    //protected DE405 jpl_ephemeris;
    /**
     * Cross sectional (reflective) area [m^2]
     */
	private double area;
	/** Mass [kg] of the satellite
     */
	private double mass;
	/** Satellite coefficient of reflectivity 
     */
	private double CR;
    
    /** All solar pressure radiation forces are determined based on
     * position relative to the sun. */
    private ReferenceFrame sunRef = new BodyCenteredInertialRef(DE405_Body.SUN);
    
    /** A reference frame relative to a planet casting a shadow.  Will be used
     * to determine the spacecraft's position relative to the planet to determine
     * the shadow. */
    private final ReferenceFrame shadowRef;
    
    /** The radius of planet casting shadow.  Or ignored if shadowRef is null. */
    private final double shadowR;
	
	/** Constructor.  Assumes an Earth shadow.
     * @param m New satellite mass value [kg].
     * @param A New satellite area value [m^2].
     * @param coeff Coefficient of reflectivity
	 */
	public SolarRadiationPressure(double m, double A, double coeff){
      this(m, A, coeff, new BodyCenteredInertialRef(DE405_Body.EARTH), Constants.R_Earth);
	}
	
    /** Constructor
     * @param m New satellite mass value [kg].
     * @param A New satellite area value [m^2].
     * @param coeff Coefficient of reflectivity
     * @param bodyRef reference frame that can tell the spacecraft's position
     * relative to a planet that may cast a shadow on the spacecraft.
     * @param bodyRadius radius of the planet in meters
     */
    public SolarRadiationPressure(double m, double A, double coeff, 
          ReferenceFrame bodyRef, double bodyRadius){
        mass = m;
        area = A;
        CR = coeff;
        shadowRef = bodyRef;
        shadowR = bodyRadius;        
    }
    
	/**
	 * Constructor
	 * @param sc Spacecraft Parameters
	 */
	public SolarRadiationPressure(Spacecraft sc){
	    this(sc.mass(), sc.srpArea(), sc.cr());
	}
	
    /** Compute the acceleration due to a solar radiation pressure.
     * @param d Relative position vector of spacecraft w.r.t. Sun 
     * (from the sun to s/c) (meters)
     * @return Acceleration due to solar radiation pressure in m/s^2.
     */
    public VectorN accelSRP(VectorN d) 
    {
        double dmag = d.mag();
        double dcubed = dmag * dmag * dmag;
        //double au2 = Constants.AU * Constants.AU;
        //double r_sun_mag = r_Sun.mag();
        //double au2 = r_sun_mag*r_sun_mag;
        //double P_STK = 4.5344321837439e-06;
        //4.560E-6
        //THE REAL WAY
        double Ls = 3.823e26; //* STK [W]
        double factor = CR * (area/mass) * Ls/ (4*Constants.pi*Constants.c*dcubed);   //* STK HPOP method
        //ALTERNATE
        //double factor = CR * (area/mass) * P_STK * au2 / dcubed;
        //double 
        
        //THE REAL WAY
        VectorN out = d.times(factor);
        //ALTERNATE
        //VectorN out = r_Sun.times(factor);
        
        //VectorN out = d.unitVector().times(P_STK*CR*(area/mass));
//        out.print("accelSRP");

        return  out;
    }
    
    /** Compute the acceleration due to a solar radiation pressure.
     * @param r ECI position vector (meters).
     * @param r_Sun ECI position vector of Sun in (meters).
     * @return Acceleration due to solar radiation pressure in m/s^2.
     */
    public VectorN accelSRP(VectorN r, VectorN r_Sun) 
    {
      // Relative position vector of spacecraft w.r.t. Sun (from the sun to s/c)
      VectorN d = r.minus(r_Sun);
      return accelSRP(d);
    }

    /** Implemented from the ForceModel interface
     * @param t Time reference object
     * @param bRef Earth reference object
     * @param sc Spacecraft parameters
     */
    public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc) {        
        area = sc.srpArea();
    	mass = sc.mass();
    	CR = sc.cr();
        
        // We translate to a sun-centered reference frame to
        // determine the distance from the sun
        ReferenceFrameTranslater xlater = 
          new ReferenceFrameTranslater(bRef, sunRef, t);
        VectorN d = xlater.translatePoint(sc.r());
        
    	return (accelSRP(d).times(partial_illumination_rel(sc.r(), d)));
    }
    
    /** Update satellite mass.
     * @param m New satellite mass value [kg].
     */
    public void updateMass(double m){
    	mass = m;
    }

    /** Update satellite area.
     * @param A New satellite area value [m^2].
     */
    public void updateArea(double A){
    	area = A;
    }
    
    public void updateCR(double Cr){
    	CR = Cr;
    }
    
    /** Determines if the satellite is in sunlight or shadow based on simple cylindrical shadow model.
     * @param r ECI position vector of spacecraft [m].
     * @param r_Sun Sun position vector (geocentric) [m].
     * @return 0.0 if in shadow, 1.0 if in sunlight.
     */
    public double illumination( VectorN r, VectorN r_Sun ) {

        VectorN e_Sun = r_Sun.unitVector();   // Sun direction unit vector
        double s     = r.dotProduct(e_Sun);      // Projection of s/c position
        VectorN temp1 = e_Sun.times(s);
        VectorN temp2 = r.minus(temp1);
        double temp = temp2.mag();

        return ( ( s>0 || temp> Constants.R_Earth ) ?  1.0 : 0.0 );
    }

    /** Determines if the satellite is in sunlight or shadow based on simple cylindrical shadow model.
     * Taken from Montenbruck and Gill p. 80-83
     * @param r ECI position vector of spacecraft [m].
     * @param r_Sun Sun position vector (geocentric) [m].
     * @return 0.0 if in shadow, 1.0 if in sunlight, 0 to 1.0 if in partial shadow
     */
    public double partial_illumination(VectorN r, VectorN r_Sun ){
      VectorN d = r.minus(r_Sun);
      return partial_illumination_rel(r, d);
    }
    
    /** Determines if the satellite is in sunlight or shadow based on simple cylindrical shadow model.
     * Taken from Montenbruck and Gill p. 80-83
     * @param r_earth position vector of spacecraft with respect to Earth (ECI, meters).
     * @param r_sun position vector of spacecraft with respect to the sun (meters).
     * @return 0.0 if in shadow, 1.0 if in sunlight, 0 to 1.0 if in partial shadow
     */
    public double partial_illumination_rel(VectorN r_earth, VectorN r_sun){
        double r_mag = r_earth.mag();
        double R_sun = Constants.R_Sun;
        // TODO: This is hard coded to Earth.  A planet reference
        // frame and radius should be passed in to the SolarRadiationPressure
        // at construction.
        double R_earth = Constants.R_Earth;
        double dmag = -1 * r_sun.mag();
        double sd = -1.0 * r_earth.dotProduct(r_sun);
        double a = Math.asin(R_sun/dmag);
        if(R_earth>r_mag){
        	System.err.println("Error! Collision detected with Earth.");
        	System.err.println("Error! See jat.forces.SolarRadiationPressure");
        	return 0.0;
        	//System.exit(0);
        }
        double b = Math.asin(R_earth/r_mag);
        double c = Math.acos(sd/(r_mag*dmag));
        if((a+b)<=c) 
            return 1.0;
        else if(c < (b-a)) return 0.0;
        else {
            double x = (c*c+a*a-b*b)/(2*c);
            double y = Math.sqrt(a*a-x*x);
            double A = a*a*Math.acos(x/a)+b*b*Math.acos((c-x)/b)-c*y;
            double nu = 1 - A/(Constants.pi*a*a);
            return nu;
        }
        
    }

}
