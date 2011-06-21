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
import jat.matvec.data.matrixDecompositions.*;
//import jat.gps.*;
/**
 * <P>
 * The ReceiverModel Class provides a GPS receiver model for GPS measurement generation.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ReceiverModel {
        
    //private static final double sf = 1.0E-19 * GPS_Utils.c * GPS_Utils.c;
    private static final double sf = 0.0036 ;
    //private static final double sg = 4.0E-20 * GPS_Utils.c * GPS_Utils.c * MathUtils.PI * MathUtils.PI;
    private static final double sg = 7.106E-05;
    
	private RandomNumber codeMeasNoise;
	private static final double codeMeasSigma = 6.0;
	private RandomNumber cpMeasNoise;
	private static final double cpMeasSigma = 0.02;

	private GaussianVector wgn;
	
	
	private static final double hc = -2.87956633585E-10 * GPS_Utils.c;
	    

    /** Default constructor. */
    public ReceiverModel() {
		VectorN zeroMean = new VectorN(2);
		// receiver clock noise
		VectorN one = new VectorN(2);
		one.set(0, 1.0);
		one.set(1, 1.0);
		this.wgn = new GaussianVector(zeroMean, one);
		
		this.codeMeasNoise = new RandomNumber();
		this.cpMeasNoise = new RandomNumber();		 	
    }
    /**
     * Constructor
     * @param seed long containing seed to be used
     */
    public ReceiverModel(long seed) {
		VectorN zeroMean = new VectorN(2);
		// receiver clock noise
		VectorN one = new VectorN(2);
		one.set(0, 1.0);
		one.set(1, 1.0);
		this.wgn = new GaussianVector(zeroMean, one, seed);
		this.codeMeasNoise = new RandomNumber(seed-1 );
		this.cpMeasNoise = new RandomNumber(seed-2);		 	
		  	
    }    

    /** Outputs the gyro's contribution to the INS error quaternion.
     * @param qref Quaternion containing the current reference quaternion.
     * @param wref Vector3 containing the reference angular velocity vector.
     * @param bg Vector3 containing the current gyro bias state.
     * @return Error quaternion due to gyro errors.
     */
    public double clockError(VectorN los, VectorN v, VectorN vGPS, double clockBias){
		double rr = GPS_Utils.rangeRate(los, v, vGPS);     	
		double error =  (1.0 - (rr/GPS_Utils.c)) * clockBias; 
//		double error = noise; 
//		System.out.println("rr = "+rr+" noise = "+noise); 	
        return error;
    }
    
    /**
     * Return the measurement noise in the pseudorange measurement
     * @return the measurement noise in the pseudorange measurement
     */
    public double codeNoise() {
		double noise = codeMeasNoise.normal(0.0, codeMeasSigma);
    	return noise;
    }
    
    /**
     * Return the measurement noise in the carrier phase range measurement
     * @return the measurement noise in the carrier phase range measurement
     */
   public double cpNoise(){
    	double noise = cpMeasNoise.normal(0.0, cpMeasSigma);
    	return noise;

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
	 * Propagate the clock states
	 * @param dt time step
	 * @param in double[] containing the previous clock states
	 * @return double[] containing the next clock states
	 */
	public double[] propClock(double dt, double[] in){
		VectorN clk = new VectorN(in);
		clk.checkVectorDimensions(2);
		Matrix phi = new Matrix(2);
		//if (dt != 1.0) System.out.println("ReceiverModel.propagate: dt = "+dt);
		phi.set(0, 1, dt);    
    	Matrix q = this.biasQ(dt);
    	CholeskyDecomposition chol = new CholeskyDecomposition(q);
    	Matrix sqrQ = chol.getL();
    	this.wgn.nextSet();
    	VectorN w = sqrQ.times(this.wgn);
    	VectorN out = phi.times(clk);
    	out = out.plus(w);
    	VectorN hvec = new VectorN(2);
    	hvec.set(0, hc*dt);
    	out = out.plus(hvec);
    	return out.x;
	}
    
    /**
     * Return the noise strength to be used in
     * the process noise matrix Q. 
     * Used by GPS_INS Process Model
     * @param dt time step size
     * @return noise strength
     */
    public Matrix biasQ(double dt) {
    	Matrix out = new Matrix(2,2);
    	out.set(0,0, sf*dt + sg*dt*dt*dt/3.0);
    	out.set(0,1, sg*dt*dt/2.0);
    	out.set(1,0, sg*dt*dt/2.0);
    	out.set(1,1, sg*dt);
    	return out;
    }
    
		   	
}
