package jat.demo.Rendezvous;

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
 * File Created on Aug 25, 2003
 */
 
import jat.cm.*;
import jat.matvec.data.*;
import jat.alg.integrators.*;
//import jat.plot.*;
 
public class CW_Rendezvous {

	public static void main(String[] args) {
		

		// set up the target orbit
		double a = 6765500.0;
		double mu = Constants.GM_Earth;
		TwoBody tgt = new TwoBody( mu, a, 0.0, 0.0, 0.0, 0.0, 0.0);

		// create the CW model		
		ClohessyWiltshire cw = new ClohessyWiltshire(tgt);

		// set up the trajectory plot
		cw.traj_plot.setTitle("Two Body Trajectory");
		cw.traj_plot.plot.setXLabel("y (km)");
		cw.traj_plot.plot.setYLabel("x (km)");

		// create an RungeKutta8 integrator with step-size of 10 sec
		RungeKutta8 rk8 = new RungeKutta8(1.0);
		
		// set up a LinePrinter for printing results to a file
//		LinePrinter lp = new LinePrinter("c:\\temp\\", "sts-iss.dat");


		// set up the initial conditions
		VectorN dr0 = new VectorN(0.0, -15000.0, 0.0);
		VectorN dv0 = new VectorN(3);
		
		// set up the times
		double tof = 5537.0;
		double tcoast = 0.0;
		double t0 = 0.0;
		double tf = tcoast + tof;


		// initialize the variables
		double[] x0 = new double[6];
		x0[0] = dr0.x[0];
		x0[1] = dr0.x[1];
		x0[2] = dr0.x[2];
		x0[3] = dv0.x[0];
		x0[4] = dv0.x[1];
		x0[5] = dv0.x[2];

		// integrate the equations to tcoast
		double[] rv = rk8.integrate(t0, x0, tcoast, cw, cw, true);
				
		VectorN drminus = new VectorN(rv[0], rv[1], rv[2]);
		VectorN dvminus = new VectorN(rv[3], rv[4], rv[5]);
		
		// compute the delta-v's		
		double dvtot = cw.twoImpulseRendezvous(tof, tof, drminus, dvminus);
		System.out.println("dvtot = " + dvtot);
		cw.deltav1.print("delta-v1");
		cw.deltav2.print("delta-v2");
		
		// add in the effects of the first deltav
		VectorN dv0plus = dvminus.plus(cw.deltav1);
		x0[0] = drminus.x[0];
		x0[1] = drminus.x[1];
		x0[2] = drminus.x[2];
		x0[3] = dv0plus.x[0];
		x0[4] = dv0plus.x[1];
		x0[5] = dv0plus.x[2];

		// integrate the equations from tcoast to tf       
		rk8.integrate(tcoast, x0, tf, cw, cw, true);

		// make the plot visible
		cw.traj_plot.setVisible(true);

		// close the LinePrinter
//		lp.close();

	}
		

}
