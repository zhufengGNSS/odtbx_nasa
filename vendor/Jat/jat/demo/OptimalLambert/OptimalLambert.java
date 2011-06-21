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

package jat.demo.OptimalLambert;

import jat.alg.opt.*;
import jat.cm.*;
import jat.matvec.data.VectorN;
import jat.plot.*;
import java.io.*;

/**
 * @author Tobias Berthold
 *
 */
public class OptimalLambert
{
	static SinglePlot plot;
	static PerformanceIndex J;
	static Lambert earth_mars;
	static double muforthisproblem;
	
	public OptimalLambert()
	{}

	/** Runs the example.
	 * @param args Arguments.
	 */

	public static void main(String[] args) throws IOException
	{
		System.out.println("Optimization Lambert targeting example");
		muforthisproblem = Constants.GM_Sun/1.e9;

		// Create instances of the classes
		earth_mars=new Lambert(muforthisproblem);
		J = new PerformanceIndex(earth_mars);

		// Find the optimal parameters for transfer
		double[] x_guess = new double[3]; // initial guess of search vector
		double[] x_search={0,0,0}; //  search vector x=[ tf theta0 theta1 ]
		double tof=190. * 86400.0;
		x_guess[0] = tof;
		x_guess[1] = -10.;
		x_guess[2] = 120.;

		// For a new optimization problem: test the performance index calculation,
		// then create a grid to find good initial guesses, then run the optimization
		int mode = 4;

		if (mode == 0)
		{
			print_result(x_guess);
			plot_trajectories(x_guess);
		}

		if (mode == 1)
		{
			J.initgrid(.8*tof,1.2*tof,0.1*tof);
			J.create_combinations();
		}
		if (mode == 2)
		{
			J.initgrid(.7*tof,1.3*tof,0.1*tof);
			J.print_header();
			J.create_grid();
		}
		if (mode == 3)
		{
			GradientSearch OST = new GradientSearch(J, x_guess);
			OST.err_ods = 1.e-4;
			OST.err_grad = 1.e-1;
			OST.eps_FD = 1.e-8; //-8;
			x_search = OST.find_min_gradient();
			print_result(x_search);
			plot_trajectories(x_search);
		}
		if (mode == 4)
		{
			DFP OST = new DFP(J, x_guess);
			OST.err_ods = 1.e-4;
			OST.err_dfp = 1.e-6;
			OST.eps_CD = 1.e-4;
			OST.max_it = 100;
			x_search = OST.find_min_DFP();
			print_result(x_search);
			plot_trajectories(x_search);
		}

		System.out.println("Done");

	}

	static void plot_trajectories(double[] x)
	{
		// Plot the trajectory
		// Set up the trajectory plot
		plot = new SinglePlot();
		plot.setTitle("Optimal Lambert targeting");
		plot.plot.setXLabel("x (km)");
		plot.plot.setYLabel("y (km)");
		double plotsize = 300000000;
		plot.plot.setXRange(-plotsize, plotsize);
		plot.plot.setYRange(-plotsize, plotsize);
		plot.plot.setMarksStyle("dots");
		plot.plot.addPoint(2, 0., 0., false);
		plot.plot.setMarksStyle("none");
		plot.setVisible(true);	// Make the plot visible

		// Create two body orbits (elliptical)
		Orbit orbit_initial=new Orbit(muforthisproblem,cm.earth_moon_elements,0);
		Orbit orbit_final=new Orbit(muforthisproblem,cm.mars_elements,1);
		// Plot the initial and the final orbit
		orbit_initial.plot_orbit(plot);
		orbit_final.plot_orbit(plot);

		orbit_initial.twobody.setTa(x[1]);
		VectorN dv0 = earth_mars.deltav0;
		VectorN r0 = orbit_initial.twobody.getR();
		VectorN v0 = orbit_initial.twobody.getV();
		// apply the first delta-v
		v0=v0.plus(dv0);
		Orbit orbit_sc=new Orbit(muforthisproblem, r0,v0,2);
		orbit_sc.plotperiod=x[0];
		orbit_sc.plot_orbit(plot);
		plot.plot.setMarksStyle("dots", 3);
		plot.plot.addPoint(3, orbit_initial.twobody.getR().x[0], orbit_initial.twobody.getR().x[1], false);
		plot.plot.addPoint(3, orbit_sc.twobody.getR().x[0], orbit_sc.twobody.getR().x[1], false);

		System.out.println("tof = " + earth_mars.tof);
		//TwoBody spacecraft = new TwoBody(muforthisproblem, r0, v0);
		//spacecraft.print("spacecraft orbit");
		//spacecraft.propagate(0.0, tof, this, true);

	}
	
	static void print_result(double[] x_search)
	{
		double result = J.evaluate(x_search);
		System.out.println("Result");
		System.out.println("J= "+result);
		new VectorN(x_search).print();

	}

}
