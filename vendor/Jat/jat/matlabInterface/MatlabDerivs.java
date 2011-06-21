/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
package jat.matlabInterface;

import jat.alg.integrators.*;
import jat.matvec.data.VectorN;

import com.mathworks.jmi.Matlab;

/**
 * MatlabDerivs allows the user to pass a derivs function implemented in Matlab
 * to JAT numerical integrators. 
 * 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * 
 */
public class MatlabDerivs implements Derivatives {

	private String cmd = null;

	private Object[] inputArgs = null;

	private boolean verbose = false;

	/**
	 * @param command String containing Matlab command for computing the derivatives
	 */
	public MatlabDerivs(String command) {
		cmd = command;
	}

	/**
	 * Set to verbose output to screen 
	 *
	 */
	public void verbose() {
		verbose = true;
	}

	/**
	 * Provide derivative equations to the numerical integrator
	 * @param t independent variable
	 * @param x dependent variables
	 * 
	 * @see jat.alg.integrators.Derivatives#derivs(double, double[])
	 */
	public double[] derivs(double t, double[] x) {
		// Set up the inputs
		int nargs = 2;
		inputArgs = new Object[nargs];
		inputArgs[0] = new Double(t);
		Double[] temp = new Double[x.length];
		for (int i = 0; i < x.length; i++) {
			temp[i] = new Double(x[i]);
		}
		inputArgs[1] = temp;

		// Get the return values from matlab
		double[] returnVals = null;
		if (cmd.equals("JatUniverseJGM2")) {
			/* In here I want to access the EOMs for the JGM2 All Case 
			 * and I want to assign the return from the EOMs to be the 
			 * returnVals. */
		} else if (cmd.equals("JatUniverseJGM3")) {
			/* In here I want to access the EOMs for the JGM3 All Case
			 * and I want to assign the return from the EOMs to be the 
			 * returnVals. */
		} else {
			try {
				if (verbose) {
					returnVals = (double[]) Matlab.mtFevalConsoleOutput(cmd,
							inputArgs, 0);
				} else {
					returnVals = (double[]) Matlab.mtFeval(cmd, inputArgs, 0);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return returnVals;

	}

	public void test(String command) {
		MatlabDerivs d = new MatlabDerivs(command);
		System.out.println("entering test");
		double t = 60;
		double[] x = new double[2];
		x[0] = 1.0;
		x[1] = 2.0;
		System.out.println("calling derivs");
		d.derivs(t, x);
		double[] out = d.derivs(t, x);
		VectorN output = new VectorN(out);
		output.print("final output vector");
	}

}
