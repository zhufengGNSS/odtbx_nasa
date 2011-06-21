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
import jat.plot.*;
import jat.alg.integrators.*;

/**
 * <P>
 * Threebody dynamics simulation
 * @author Tobias Berthold
 * @version 1.0
 */

public class ThreeBodyExample implements Printable
{

	TwoPlots traj_plot = new TwoPlots();

	/** Creates a new instance of threebody */
	public ThreeBodyExample()
	{
		// set up the trajectory plot
		traj_plot.setTitle("Three body Simulation");
		traj_plot.topPlot.setXLabel("x");
		traj_plot.topPlot.setYLabel("y");
		traj_plot.bottomPlot.setXLabel("x");
		traj_plot.bottomPlot.setYLabel("y");

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
		//System.out.println(t + " " + y[0] + " " + y[1] + " " + y[2] );

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
		ThreeBodyExample tbe = new ThreeBodyExample();
		
		// Threebody dynamics
		ThreeBody tb= new ThreeBody(1.,30.,30.,30.);		

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
		rk8.integrate(t0, x0, tf, tb, tbe, true);

		// make the plot visible
		tbe.traj_plot.setVisible(true);
	}
}
