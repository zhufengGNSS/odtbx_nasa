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
 * File Created on May 9, 2003
 */

package jat.alg.estimators;
import jat.matvec.data.*;

/**
* The EstSTM Class contains the state and state transition matrix (STM) for
* estimation problems. It also has some methods to help convert between
* the big array used for propagation and the state vector and STM.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class EstSTM {

	/** state vector */
	private VectorN state;

	/** state transition matrix */
	private Matrix phi;

	/** number of states */
	private int numberOfStates;

	/** total number of variables = n + n^2 */
	private int total;

	/**
	 * Constructor
	 * @param in integer containing number of states
	 */
	public EstSTM(int in) {
		this.numberOfStates = in;
		int n = this.numberOfStates;
		this.total = n + (n * n);
	}

	/**
	 * Constructor
	 * @param in VectorN containing only the state
	 */
	public EstSTM(VectorN in) {
		this.state = in.copy();
		this.phi = new Matrix(in.length);
		this.numberOfStates = in.length;
		int n = this.numberOfStates;
		this.total = n + (n * n);
	}

	/**
	 * Constructor
	 * @param in double[] containing the state and state transition matrix
	 * @param n number of states
	 */
	public EstSTM(double in[], int n) {
		this.numberOfStates = n;
		this.total = n + (n * n);

		if (in.length < this.total) {
			System.out.println("EstSTM: trying to convert a too small vector");
			System.exit(1);
		}
		double[] x = new double[n];
		this.phi = new Matrix(n, n);

		for (int i = 0; i < n; i++) {
			x[i] = in[i];
		}
		this.state = new VectorN(x);

		int k = n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				this.phi.A[i][j] = in[k];
				k = k + 1;
			}
		}
	}

	/**
	 * Stuff the current state and phi into a big array for the integrator
	 * @return double[] containing the current state and phi
	 */
	public double[] longarray() {
		int n = this.numberOfStates;
		double[] out = new double[this.total];
		for (int i = 0; i < n; i++) {
			out[i] = this.state.x[i];
		}

		int k = n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				out[k] = this.phi.A[i][j];
				k = k + 1;
			}
		}
		return out;
	}

	/**
	 * Adds the state correction to the current state vector
	 * @param xhat VectorN containing the state corrections.
	 */
	public void update(VectorN xhat) {
		xhat.checkVectorDimensions(this.numberOfStates);
		this.state = this.state.plus(xhat);
		
		// fix the quaternion, this is a kluge!!!
//		VectorN q = this.get(6, 4);
//		q.unitize();
//		this.setState(6, q);
		
	}
	
	/**
	 * Set the state to a particular value
	 * @param index int containing the index
	 * @param values VectorN containing the values
	 */	
	public void setState(int index, VectorN values) {
		this.state.set(index, values);
	}
	
	/**
	 * Get the state values
	 * @param index int containing the index
	 * @param count int containing the number of values to get
	 */	
	public VectorN get(int index, int count) {
		VectorN out = this.state.get(index, count);
		return out;
	}

	/**
	 * Prints the current state and phi
	 * @param title String containing a title
	 */
	public void print(String title) {
		System.out.println(title);
		this.state.print("state");
		this.phi.print("phi");
	}

	/**
	 * Get the entire state vector
	 * @return VectorN containing the state
	 */
	public VectorN state() {
		return this.state;
	}

	/**
	 * Get the state transition matrix
	 * @return Matrix containing the state transition matrix
	 */
	public Matrix phi() {
		return this.phi;
	}
	
	/**
	 * Get the number of states 
	 * @return int containing the number of states
	 */
	public int numberOfStates(){
		return this.numberOfStates;
	}
	
	/**
	 * Get the number of variables (states + state transition matrix elements)
	 * @return int containing the number of variables
	 */
	public int total(){
		return this.total;
	}
	
	/**
	 * Set the state transition matrix to the Identity matrix
	 */
	public void resetPhi(){
		int n = this.numberOfStates;
		this.phi = new Matrix(n);
	}

}
