/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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
* File Created on September 12, 2007
*/
package jat.ground_tracking;
import jat.eph.DE405;
import jat.eph.DE405_Body;
import jat.matvec.data.VectorN;
import jat.spacetime.Time;

/**
 * This class calculates the light time delay due to charged particles.
 * 
 * Source:
 * <b><i>
 * Burkhart, P. D., "Adaptive Orbit Determination For Interplanetary Spacecraft",
 * PhD Dissertation, The University Of Texas At Austin, May 1995, section 2.2.3
 * "Ionosphere Model", pp. 22-27.
 * </i></b>
 */
public class ChargedParticleModel {

	
	final double C0 = 1.5e17;
	
	/** Local electron density at unit heliocentric distance in electrons/cm^2
	 */
	double N1 = 4.3; 
	
	
	/** Provide access to the Sun's Position
	 */
	private DE405 jpl_ephem;
	
	/** Frequency the ground station is using
	 */
	private double freq;
	
	
	public ChargedParticleModel(double frequency){
		
		jpl_ephem = new DE405();
		this.freq = frequency;
		
		
	}
	
	/**
	 * Calculates the delay distance.
	 * @param satPos ECI satellite coordinates (m)
	 * @param t Epoch
	 * @return Signal delay error due to charged particles 
	 * expressed as distance (m)
	 */
	/*Note.  Wanted this to be called directly to avoid calling the 
	 * constructor more than once.  (I thought the DE405 construction
	 * took too long)
	 */
	public double returnDelay(VectorN satPos, Time t)
	{
		double TEC = computeTEC(t,satPos);
		double delay = computeDelay(TEC);

		return delay;
	}
	
	/**
	 * 
	 * @param t Epoch
	 * @param satEarthVec ECI satellite position (m)
	 * @return
	 */
	private double computeTEC(Time t, VectorN satEarthVec)
	{		
		
		//Determine the position of the Sun (m)
		VectorN r_sun = new VectorN(jpl_ephem.get_planet_pos(DE405_Body.GEOCENTRIC_SUN, t.mjd_tt()));
		
		//The problem formulation defines the sun-to-Earth Vector, so flip this around
		r_sun = r_sun.times(-1.0);
		
		//Get the Sun - to - satellite vector using vector addition
		VectorN satSunVec = new VectorN(satEarthVec.plus(r_sun));
		
		//All distances need to be in AU
		satEarthVec = satEarthVec.divide(jat.cm.Constants.AU);
		r_sun		= r_sun.divide(jat.cm.Constants.AU);
		satSunVec   = satSunVec.divide(jat.cm.Constants.AU);
		
		//Compute the Ranges
		double sunNorm = r_sun.mag();
		double satSunVecNorm = satSunVec.mag();
		double satEarthVecNorm = satEarthVec.mag();
		
		//Determine the sun-earth-satellite angle
		double num = satEarthVec.times(-1.0).dotProduct(r_sun);
		double denom = satEarthVecNorm*sunNorm;
		double alpha = Math.acos(num/denom);
		
		//Determine the relative longitude of the spacecraft and Earth
		num = satSunVec.dotProduct(r_sun);
		denom = satSunVecNorm*sunNorm;
		double lambda = Math.acos(num/denom);	
		
//		Compute the TEC for the nominal and boundary conditions
		double TEC = 0;
		if(alpha > 0 && alpha < Math.PI)
		{
			TEC = (C0*N1/sunNorm)*lambda/Math.sin(alpha);		
		}
		else if (alpha == 0)
		{
			TEC = (C0*N1/sunNorm)*(satEarthVecNorm)/(sunNorm - satEarthVecNorm);
		}
		else if(alpha == Math.PI/180)
		{
			TEC = (C0*N1/sunNorm)*(satEarthVecNorm)/(sunNorm + satEarthVecNorm);
		}
		
		return TEC ;
	}
	
	private double computeDelay(double TEC)
	{
		
		double delay = 40.3*TEC/(this.freq*this.freq);
		
		return delay;
	}
	
	public void setN1(double _N1)
	{
		this.N1 = _N1;
	}
	
	public double getN1()
	{
		return this.N1;
	}

	public static void main(String[] args) {
		

		double frequency = 1.6e9;
		Time T = new Time(54626);
		VectorN satPos = new VectorN(2.056020937918350e7,1.411064433075980e7,0.160945891394200e7);

		ChargedParticleModel cpm = new ChargedParticleModel(frequency);
		double delay = cpm.returnDelay(satPos,T);

		System.out.println(delay);

	}


}