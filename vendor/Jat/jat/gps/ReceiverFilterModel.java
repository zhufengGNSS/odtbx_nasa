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

package jat.gps;
import jat.matvec.data.*;
import jat.math.*;
//import jat.cm.*;
//import jat.matvec.data.matrixDecompositions.*;
//import jat.gps.*;
/**
 * <P>
 * The ReceiverModel Class provides a GPS receiver model for the Kalman Filter.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ReceiverFilterModel {
        
    private static final double sf = 1.0E-19 * GPS_Utils.c * GPS_Utils.c;
    
    private static final double sg = 4.0E-20 * GPS_Utils.c * GPS_Utils.c * MathUtils.PI * MathUtils.PI;
    	
	private static final double hc = -2.87956633585E-10 * GPS_Utils.c;
	    
    /** Default constructor. */
    public ReceiverFilterModel() {
    }
    
    /** 
     * Compute the derivatives for the clock bias  and drift states.
     * Used by GPS_INS Process Model
     * @param bc clock bias vector
     * @return the time derivative of the gyro bias
     */
    public VectorN biasProcess(VectorN xin) {
    	
		double b = xin.x[0];
		double f = xin.x[1];
    	VectorN out = new VectorN(2);
    	out.set(0, (f + hc));
    	
    	return out;
    }

	
//	public double[] propClock(double dt, double[] in){
//		double bk = in[0];
//		double fk = in[1];
//		double g = 2.0E-11 * k.c;
//		double bkp1 = bk + (fk + this.hc)*dt + 0.5*g*dt*dt;
//		double fkp1 = fk + g*dt;
//		VectorN out = new VectorN(2);
//		out.set(0, bkp1);
//		out.set(1, fkp1);
//    	return out.x;
//	}

    
    /**
     * Return the gyro bias noise strength to be used in
     * the process noise matrix Q. Assumes 1 second time step.
     * Used by GPS_INS Process Model
     * @return gyro bias noise strength
     */
    public Matrix biasQ(double dt) {
    	Matrix out = new Matrix(2,2);
    	out.set(0,0, sf*dt + sg*dt*dt*dt/3.0);
    	out.set(0,1, sg*dt*dt/2.0);
    	out.set(1,0, sg*dt*dt/2.0);
    	out.set(1,1, sg*dt);
    	return out;
    }

    /**
     * Return the noise strength to be used in
     * the process noise matrix Q. Assumes 1 second time step.
     * @return noise strength
     */
    public Matrix Q() {
    	Matrix out = new Matrix(2,2);
    	out.set(0,0, sf);
    	out.set(1,1, sg);
    	return out;
    }
    		   	
}
