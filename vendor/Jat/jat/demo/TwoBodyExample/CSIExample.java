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

package jat.demo.TwoBodyExample;

import jat.alg.*;
import jat.alg.integrators.*;
import jat.cm.*;
import jat.cm.eom.*;
import jat.matvec.data.VectorN;
import jat.plot.*;
import java.io.*;

/** This demo computes the trajectory of a spacecraft with a constant specific 
 * impulse propulsion system with the thrust pointing in the velocity direction. 
 * The trajectory is plotted.  
 * @author Tobias Berthold
 *
 */
public class CSIExample
{
	static SinglePlot traj_plot;

	public CSIExample()
	{
		// set up the trajectory plot
		traj_plot.setTitle("CSI Trajectory in a Two Body Force Field");
		traj_plot.plot.setXLabel("x (km)");
		traj_plot.plot.setYLabel("y (km)");
	}

	public class Initial implements Printable
	{
		public void print(double t, double[] y)
		{
			boolean first = true;
			if (t == 0.0)
				first = false;

			// print to the screen for warm fuzzy feeling
			//System.out.println(t + " " + y[0] + " " + y[1] + " " + first);

			// add data point to the plot
			traj_plot.plot.addPoint(0, y[0], y[1], first);
		}
	}

	public class CSI implements Printable, VectorTimeFunction
	{
		public void print(double t, double[] y)
		{
			boolean first = true;
			if (t == 0.0)
				first = false;

			// print to the screen for warm fuzzy feeling
			//			System.out.println(t + " " + y[0] + " " + y[1] + " " + y[6]);

			// add data point to the plot
			traj_plot.plot.addPoint(1, y[0], y[1], first);
		}

		public VectorN evaluate(VectorN x, double t)
		{
			VectorN u = new VectorN(3);
			// Constant thrust in velocity direction 
			u.x[0] = x.x[3];
			u.x[1] = x.x[4];
			u.x[2] = 0.;
			// Constant thrust in x-direction
			//			u.x[0]=1.
			//			u.x[1]=0.;
			// Constant thrust in radial direction
			//			u.x[0] = x.x[0];
			//			u.x[1] = x.x[1];
			//System.out.println(tf.format(t)+"  "+df.format(ux.get_value(t))+" "+uy.get_value(t));
			return u.unitVector();
		}
	}

	/** Runs the example.
	 * @param args Arguments.
	 */
	public static void main(String[] args) throws IOException
	{
		System.out.println("Continuous thrust trajectory example");
		traj_plot = new SinglePlot();

		// create instances of the classes
		CSIExample example = new CSIExample();
		Initial in = example.new Initial();
		CSI csi = example.new CSI();

		// create Two Body orbits (elliptical)
		KeplerElements k_initial = new KeplerElements(10000.0, 0., 0.0, 0.0, 0.0, 0.0);
		TwoBody initialorbit = new TwoBody(cm.mu, k_initial);
		TwoBodyCSI trajectory = new TwoBodyCSI(k_initial, cm.mu, 3000., csi, 0.05, 10000);

		// create an RungeKuttaFehlberg78 integrator
		RungeKutta8 rk8 = new RungeKutta8();
		rk8.setStepSize(1.e2);

		// initialize the state variables and times
		double[] x0initial = initialorbit.randv();
		double[] x0c = trajectory.rvm.getArray();
		double t0 = 0.0;
		double tf = initialorbit.period();
		double tf_CSI = 170000.; // final time for transfer trajectory

		// make the plot visible
		double plotsize = 3.e4;
		traj_plot.plot.setXRange(-plotsize, plotsize);
		traj_plot.plot.setYRange(-plotsize, plotsize);
		traj_plot.setVisible(true);

		// integrate the equations
		rk8.integrate(t0, x0initial, tf, initialorbit, in, true);
		double[] finalstate = rk8.integrate(0., x0c, tf_CSI, trajectory, csi, true);

	}
}
