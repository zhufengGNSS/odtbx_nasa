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
 */

package jat.ins;
import jat.matvec.data.*;
import jat.math.*;

/**
 * <P>
 * The SIGIGyro Class models the SIGI gyro triad.
 * Only gyro bias and measurement noise are included.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class SIGIGyroFilterModel {
        
    /** accelerometer bias correlation time in seconds */
    public static final double correlationTime = 3600.0;
            
    /** gyro bias noise strength */
    private double qbias;
    
    private double q;

    /** Default constructor. Hardcoded with SIGI gyro numbers. */
    public SIGIGyroFilterModel() {
    	
    	VectorN zeroMean = new VectorN(3);
    	
    	// gyro bias parameters for SIGI
    	double biasSigma = 0.0035 * MathUtils.DEG2RAD / 3600.0;
    	double dt = 1.0; // time step
    	double exponent = -2.0*dt/correlationTime;
    	this.qbias = biasSigma*biasSigma*(1.0 - Math.exp(exponent));  // in (rad/sec)^2/Hz
		this.q = 2.0 * biasSigma * biasSigma / correlationTime;
    }

    
    /** 
     * Compute the derivatives for the gyro bias state.
     * The gyro bias is modeled as a first order Gauss-Markov process.
     * Used by GPS_INS Process Model.
     * @param bg gyro bias vector
     * @return the time derivative of the gyro bias
     */
    public VectorN biasProcess(VectorN bg) {
    	double coef = -1.0/correlationTime;
    	VectorN out = bg.times(coef);
    	
    	return out;
    }
    
   
    /**
     * Return the gyro bias noise strength to be used in
     * the process noise matrix Q.
     * @return gyro bias noise strength
     */
    public double Q() {
    	return this.q;
    }    
        	
}
