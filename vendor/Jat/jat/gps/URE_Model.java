package jat.gps;

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
 * File Created on Jun 20, 2003
 */
 
//import jat.gps.*;
import jat.cm.*;
import jat.matvec.data.*;

/**
 * <P>
 * The URE_Model Class provides a model of the errors in the GPS system
 * due to GPS SV ephemeris and clock errors. For this simulation, they
 * are modeled as random constants with sigmas from:
 * 
 * Reference: J. F. Zumberge and W. I. Bertiger, "Ephemeris and Clock
 * Navigation Message Accuracy", Global Positioning System: Theory and
 * Applications, Volume 1, edited by Parkinson and Spilker. 
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class URE_Model {
	
	
	/** Radial ephemeris error sigma in meters */
	private static final double sigma_r = 1.2;
	
	/** Cross-track ephemeris error sigma in meters */
	private static final double sigma_c = 3.2;
	
	/** Along-track ephemeris error sigma in meters */
	private static final double sigma_a = 4.5;
	
	/** SV Clock error sigma in meters */
	private static final double sigma_t = 1.12E-08 * Constants.c;
	
	private static final double sigma_ure 
		= Math.sqrt(sigma_r*sigma_r + sigma_c*sigma_c + sigma_a*sigma_a + sigma_t*sigma_t);
		
	public static final double correlationTime = 7200.0;
	
	private double qbias;
	
	private double q;
	
	/** Radial ephemeris error vector, one entry per GPS SV */
	private VectorN dr;

	/** Crosstrack ephemeris error vector, one entry per GPS SV */
	private VectorN dc;

	/** Alongtrack ephemeris error vector, one entry per GPS SV */
	private VectorN da;

	/** SV Clock error vector, one entry per GPS SV */
	private VectorN dtc;
	
	/** Size of the GPS Constellation */
	private int size;
	
	/** Constructor
	 * @param n size of the GPS Constellation
	 */
	public URE_Model(int n) {
		this.size = n;
		dr = new VectorN(n);
		dc = new VectorN(n);
		da = new VectorN(n);
		dtc = new VectorN(n);
		RandomNumber rn = new RandomNumber();
		for (int i = 0; i < n; i++) {
			double radial = rn.normal(0.0, sigma_r);
			dr.set(i, radial);
			double crosstrack = rn.normal(0.0, sigma_c);
			dc.set(i, crosstrack);
			double alongtrack = rn.normal(0.0, sigma_a);
			da.set(i, alongtrack);
			double clock = rn.normal(0.0, sigma_t);
			dtc.set(i, clock);
		}
		
		double dt = 1.0;
		double exponent = -2.0*dt/correlationTime;
       	this.qbias = sigma_ure*sigma_ure*(1.0 - Math.exp(exponent));  // in (rad/sec)^2/Hz
		this.q = 2.0 * sigma_ure*sigma_ure / correlationTime;

	}

	/** Constructor
	 * @param n size of the GPS Constellation
	 * @param seed long containing random number seed to be used
	 */
	public URE_Model(int n, long seed) {
		this.size = n;
		dr = new VectorN(n);
		dc = new VectorN(n);
		da = new VectorN(n);
		dtc = new VectorN(n);
		RandomNumber rn = new RandomNumber(seed);
		for (int i = 0; i < n; i++) {
			double radial = rn.normal(0.0, sigma_r);
			dr.set(i, radial);
			double crosstrack = rn.normal(0.0, sigma_c);
			dc.set(i, crosstrack);
			double alongtrack = rn.normal(0.0, sigma_a);
			da.set(i, alongtrack);
			double clock = rn.normal(0.0, sigma_t);
			dtc.set(i, clock);
		}
		
		double dt = 1.0;
		double exponent = -2.0*dt/correlationTime;
       	this.qbias = sigma_ure*sigma_ure*(1.0 - Math.exp(exponent));  // in (rad/sec)^2/Hz
		this.q = 2.0 * sigma_ure*sigma_ure / correlationTime;

	}
	
	/** Compute the User range error due to SV clock and ephemeris errors.
	 * @param i GPS SV index
	 * @param los GPS line of sight vector
	 * @param rGPS GPS SV position vector
	 * @param vGPS GPS SV velocity vector
	 * @return the user range error in meters
	 */
	public double ure (int i, VectorN los, VectorN rGPS, VectorN vGPS) {
		
		// get the transformation from RIC to ECI
		TwoBody orbit = new TwoBody(Constants.GM_Earth, rGPS, vGPS);
		Matrix rot = orbit.RSW2ECI();
		
		// form the ephemeris error vector for the ith GPS SV
		VectorN error = new VectorN(this.dr.x[i], this.da.x[i], this.dc.x[i]);

		// rotate the ephemeris error to the ECI frame
		VectorN errECI = rot.times(error);
		
		// find the magnitude of the projection of the error vector onto the LOS vector
		double drho = errECI.projectionMag(los);
		
		// add the SV clock error
		double out = drho + this.dtc.x[i];
		return out;
	}
	
    /** 
     * Compute the derivatives for the URE state.
     * The URE is modeled as a first order Gauss-Markov process.
     * Used by GPS_INS Process Model.
     * @param ure URE state vector
     * @return the time derivative of the URE
     */
    public VectorN ureProcess(VectorN ure) {
    	double coef = -1.0/correlationTime;
    	VectorN out = ure.times(coef);
    	
    	return out;
    }
	
    /**
     * Return the URE noise strength to be used in
     * the process noise matrix Q.
     * @return URE noise strength
     */
    public double biasQ() {
    	return this.qbias;
    }

    /**
     * Return the URE noise strength to be used in
     * the process noise matrix Q.
     * @return URE noise strength
     */
    public double Q() {
    	return this.q;
    }

    /**
     * Return the URE sigma
     * @return URE sigma
     */
    public double sigma() {
    	return sigma_ure;
    }
		
}
