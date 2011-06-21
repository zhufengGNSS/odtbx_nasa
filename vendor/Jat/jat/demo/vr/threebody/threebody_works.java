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
import jat.plot.*;
import jat.alg.integrators.*;

/**
 * <P>
 * Threebody dynamics simulation
 * @author Tobias Berthold
 * @version 1.0
 */

public class threebody_works implements Derivatives, Printable
{

	TwoPlots traj_plot = new TwoPlots();

	/** Creates a new instance of threebody */
	public threebody_works()
	{
		// set up the trajectory plot
		traj_plot.setTitle("Three body Simulation");
		traj_plot.topPlot.setXLabel("x");
		traj_plot.topPlot.setYLabel("y");
		traj_plot.bottomPlot.setXLabel("x");
		traj_plot.bottomPlot.setYLabel("y");

	}
	/** Compute the derivatives.
	 * @param t    double containing time or the independent variable.
	 * @param x    VectorN containing the required data.
	 * @return      double [] containing the derivatives.
	 */

	/*
	public double[] derivs(double t, double[] x)
	{
	   double [] out = new double[2];
	   out[0] = x[1];
	   out[1] = -x[0];
	   return out;
	}
	*/

	/** Compute the derivatives.
	 * @param t    double containing time or the independent variable.
	 * @param x    VectorN containing the required data.
	 * @return      double [] containing the derivatives.
	 */
	// vector r1 is x[0], x[1], x[2]
	// vector v1 is x[3], x[4], x[5]
	// vector r2 is x[6], x[7], x[8]
	// vector v2 is x[9], x[10], x[11]
	// vector r3 is x[12], x[13], x[14]
	// vector v3 is x[15], x[16], x[17]

	public double[] derivs(double t, double[] x)
	{
		double m1, m2, m3;
		double r1, r2, r3;
		double r12, r23, r13;
		double r12cubed, r23cubed, r13cubed;

		double dxdt[] = new double[18];

		m1 = 30.;
		m2 = 30.;
		m3 = 30.;

		r12 = Math.sqrt((x[6] - x[0]) * (x[6] - x[0]) + (x[7] - x[1]) * (x[7] - x[1]) + (x[8] - x[2]) * (x[8] - x[2]));
		r13 = Math.sqrt((x[12] - x[0]) * (x[12] - x[0]) + (x[13] - x[1]) * (x[13] - x[1]) + (x[14] - x[2]) * (x[14] - x[2]));
		r23 = Math.sqrt((x[12] - x[6]) * (x[12] - x[6]) + (x[13] - x[7]) * (x[13] - x[7]) + (x[14] - x[8]) * (x[14] - x[8]));
		r12cubed = r12 * r12 * r12;
		r13cubed = r13 * r13 * r13;
		r23cubed = r23 * r23 * r23;

		// Derivatives
		dxdt[0] = x[3];
		dxdt[1] = x[4];
		dxdt[2] = x[5];
		dxdt[3] = m2 / r12cubed * (x[6] - x[0]) + m3 / r13cubed * (x[12] - x[0]);
		dxdt[4] = m2 / r12cubed * (x[7] - x[1]) + m3 / r13cubed * (x[13] - x[1]);
		dxdt[5] = m2 / r12cubed * (x[8] - x[2]) + m3 / r13cubed * (x[14] - x[2]);
		dxdt[6] = x[9];
		dxdt[7] = x[10];
		dxdt[8] = x[11];
		dxdt[9] = -m1 / r12cubed * (x[6] - x[0]) + m3 / r23cubed * (x[12] - x[6]);
		dxdt[10] = -m1 / r12cubed * (x[7] - x[1]) + m3 / r23cubed * (x[13] - x[7]);
		dxdt[11] = -m1 / r12cubed * (x[8] - x[2]) + m3 / r23cubed * (x[14] - x[8]);
		dxdt[12] = x[15];
		dxdt[13] = x[16];
		dxdt[14] = x[17];
		dxdt[15] = -m1 / r13cubed * (x[12] - x[0]) - m2 / r23cubed * (x[12] - x[6]);
		dxdt[16] = -m1 / r13cubed * (x[13] - x[1]) - m2 / r23cubed * (x[13] - x[7]);
		dxdt[17] = -m1 / r13cubed * (x[14] - x[2]) - m2 / r23cubed * (x[14] - x[8]);
		return dxdt;
	}

	/** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
	 *  This method is executed by the propagator at each integration step.
	 * @param t Time.
	 * @param y Data array.
	 */
	public void print(double t, double[] y)
	{

		// handle the first variable for plotting - this is a little mystery but it works
		boolean first = true;
		if (t == 0.0) first = false;

		// print to the screen for warm fuzzy feeling
		System.out.println(t + " " + y[0] + " " + y[1] + " " + y[2] );

		// add data point to the plot
		traj_plot.topPlot.addPoint(0, y[0], y[1], first);
		traj_plot.bottomPlot.addPoint(0, y[6], y[7], first);

	}

	/** Runs the example.
	 * @param args Arguments.
	 */
	public static void main(String[] args)
	{

		// create an RungeKutta8 integrator with step-size of 0.1
		RungeKutta8 rk8 = new RungeKutta8(0.01);

		// create an instance
		threebody_works tb = new threebody_works();

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
		double tf = 100.0;

		// set the initial time to zero
		double t0 = 0.0;

		// integrate the equations
		rk8.integrate(t0, x0, tf, tb, tb, true);

		// make the plot visible
		tb.traj_plot.setVisible(true);
	}
}
