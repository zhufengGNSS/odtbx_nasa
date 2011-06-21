package jat.alg.estimators;

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
 * File Created on May 7, 2003
 */
import jat.matvec.data.*;

/**
* The ProcessModel.java Class provides an interface to the Extended Kalman Filter
* to allow the user to specify a process model.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public interface ProcessModel {
	
	
	
	/** 
	 * Returns the initial covariance matrix.
	 * @return the initial covariance matrix.
	 */
	public Matrix P0 ();

	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @param x current state vector
	 * @return the process noise matrix.
	 */
	public Matrix Q (double t, double dt, EstSTM x);


	/** 
	 * Returns the initial reference state.
	 * @return the initial reference state.
	 */
	public VectorN xref0 ();
	
	/**
	 * Returns the number of states.
	 * @return the number of states.
	 */	
	public int numberOfStates();
	
	/** 
	 * Propagate the state and state transition matrix to the next measurement time.
	 * @param t0 previous time
	 * @param xin array containing state and state transition matrix at previous time.
	 * @param tf next time
	 */	
	public double[] propagate( double t0, double[] xin, double tf);
	
	/**
	 * Print out the state and covariance data
	 * @param t time
	 * @param state state vector
	 * @param covariance covariance matrix
	 */
	public void print(double t, VectorN state, Matrix covariance);
	
	/**
	 * Print out the residuals
	 * @param t time
	 * @param resid1 residual before the measurement update
	 * @param resid2 residual after the measurement update
	 */
	public void printResiduals(double t, double resid1, double resid2);
	
	/**
	 * Close all open LinePrinters
	 */
	public void closeLinePrinter();

}
