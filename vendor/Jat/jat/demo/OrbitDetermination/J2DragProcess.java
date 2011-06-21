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
 * File Created on May 8, 2003
 */

package jat.demo.OrbitDetermination;

import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

/**
* The J2DragProcess.java Class provides the dynamics model for an orbit 
* including J2 and effects of exponential drag.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class J2DragProcess implements ProcessModel {

	private int n = 12;
	private VectorN xref = new VectorN(n);
	private Matrix phi = new Matrix(n);
	private RungeKutta8 rk8 = new RungeKutta8(1.0);
	private J2DragEOM eom = new J2DragEOM();
	private LinePrinter lp1;
	private LinePrinter lp2;
	
	public J2DragProcess(LinePrinter lp, LinePrinter lp_2){
		this.lp1 = lp;
		this.lp2 = lp_2;
	}

	/** 
	 * Returns the initial covariance matrix.
	 * @return the initial covariance matrix.
	 */
	public Matrix P0() {

		// initialize sigmas
		double[] sigmas = new double[12];
		sigmas[0] = 100.0;
		sigmas[1] = 100.0;
		sigmas[2] = 100.0;
		sigmas[3] = 10.0;
		sigmas[4] = 10.0;
		sigmas[5] = 10.0;
		sigmas[6] = 1.0E+06;
		sigmas[7] = 1.0E-02;
		sigmas[8] = 0.01;
		sigmas[9] = 100.0;
		sigmas[10] = 100.0;
		sigmas[11] = 100.0;

		// square the sigmas
		for (int j = 0; j < 12; j++) {
			sigmas[j] = sigmas[j] * sigmas[j];
		}

		// initialize covariance
		VectorN sig = new VectorN(sigmas);
		Matrix p = new Matrix(sig);
		return p;
	}

	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	public Matrix Q(double t, double dt, EstSTM x) {
		return new Matrix(12, 12);
	}

	/** 
	 * Returns the initial reference state.
	 * @return the initial reference state.
	 */
	public VectorN xref0() {
		VectorN out = new VectorN(12);
		out.x[0] = 4973900.0;
		out.x[1] = -4300600.0;
		out.x[2] = 3486200.0;
		out.x[3] = 2850.0;
		out.x[4] = 5820.0;
		out.x[5] = 3470.0;
		out.x[6] = 3.986004415E+14;
		out.x[7] = 1.1926268E-03;
		out.x[8] = 0.0475;
		out.x[9] = -2517400.0;
		out.x[10] = -4198500.0;
		out.x[11] = 4076500.0;
		return out;
	}

	/**
	 * Returns the number of states.
	 * @return the number of states.
	 */
	public int numberOfStates() {
		return this.n;
	}

	/** 
	 * Propagate the state and state transition matrix to the next measurement time.
	 * @param t0 previous time
	 * @param xin array containing state and state transition matrix at previous time.
	 * @param tf next time
	 */
	public double[] propagate(double t0, double[] x, double tf) {
		double[] out = rk8.step(t0, x, eom);
		return out;
	}

	/**
	 * Returns the state transition matrix.
	 * @return the state transition matrix (after propagation).
	 */
	public Matrix phi() {
		return this.phi;
	}

	/**
	 * Returns the current reference state (after propagation).
	 * @return the current reference state.
	 */
	public VectorN xref() {
		return this.xref;
	}

	public void print(double t, VectorN state, Matrix cov) {
		VectorN sigmas = cov.diagonal();
		sigmas = sigmas.ebeSqrt();
		VectorN printvector = new VectorN(state, sigmas);
		lp1.print(t, printvector.x);

	}
	
	public void printResiduals(double t, double r1, double r2) {
		double[] y = new double[3];
		y[0] = t;
		y[1] = r1;
		y[2] = r2;
		lp2.print(y);
	} 
	
	
	public void closeLinePrinter(){
		lp1.close();
	}

}
