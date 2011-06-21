/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.demo.vr.threebody;

import jat.cm.ThreeBody;
//import jat.plot.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;

/**
 * <P>
 * Threebody dynamics simulation
 * @author Tobias Berthold
 * @version 1.0
 */

public class ThreeBody_noplot implements Printable
{
	static double[][] store_states;
	int counter = 0;
	//static VectorN[] V;

	/** Creates a new instance */
	public ThreeBody_noplot()
	{
		store_states = new double[100 + 2][18];
	}

	/** Implements the Printable interface to get the data out of the propagator
	 *  This method is executed by the propagator at each integration step.
	 * @param t Time.
	 * @param y Data array.
	 */
	public void print(double t, double[] y)
	{

		// print to the screen for warm fuzzy feeling
		System.out.println(t + " " + y[0] + " " + y[1] + " " + y[2] + " " + counter);
		for (int i = 0; i < 18; i++)
			store_states[counter][i] = y[i];
		//store_states[1][counter]=y[1];

		counter++;

		//V[0].set(y[0]);

	}

	/** Runs the example.
	 * @param args Arguments.
	 */
	public static void main(String[] args)
	{
		//V=new VectorN[10];

		// create an RungeKutta8 integrator with step-size of 0.1
		RungeKutta8 rk8 = new RungeKutta8(0.1);

		// create an instance
		ThreeBody_noplot tbe = new ThreeBody_noplot();

		// Threebody dynamics
		ThreeBody tb = new ThreeBody(1., 30., 30., 30.);

		// initialize the variables
		double[] x0 = new double[18];
		// Mass 1
		x0[0] = 10.0;
		x0[1] = 0.0;
		x0[2] = 0.0;
		x0[3] = 0.0;
		x0[4] = 1.0;
		x0[5] = 0.0;
		// Mass 2
		x0[6] = -5.0;
		x0[7] = 8.66025;
		x0[8] = 0.0;
		x0[9] = -0.866025;
		x0[10] = -0.5;
		x0[11] = 0.0;
		// Mass 3
		x0[12] = -5.0;
		x0[13] = -8.66025;
		x0[14] = 0.0;
		x0[15] = 0.866025;
		x0[16] = -0.5;
		x0[17] = 0.0;

		// set the final time
		double tf = 10.0;

		// set the initial time to zero
		double t0 = 0.0;

		// integrate the equations
		rk8.integrate(t0, x0, tf, tb, tbe, true);

		//V[0].print();
		VectorN V = new VectorN(store_states[1]);
		V.print();

	}
}
