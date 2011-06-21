/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2004 United States Government as represented by the
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
 * Created on Jun 22, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package jat.demo.OST;

import java.text.*;

import jat.alg.ScalarfromArrayFunction;
import jat.alg.integrators.*;
import jat.cm.eom.TwoBodyCSI;
import jat.matvec.data.VectorN;

/**
 * @author Tobias Berthold
 *
 */
public class PerformanceIndex implements ScalarfromArrayFunction
{

	RungeKutta8 rk;
	//RungeKuttaFehlberg78 rk;
	double[] x0;
	double tf_CSI;
	double[] x_target;
	TwoBodyCSI trajectory;
	ThrustDirection thrust;
	DecimalFormat Jf; // formatting  
	NumberFormat itf;
	double tf_var;

	public PerformanceIndex( double[] x0, double tf_CSI, double[] x_target, TwoBodyCSI trajectory, ThrustDirection thrust)
	{
		this.x0 = x0;
		this.tf_CSI = tf_CSI;
		this.x_target = x_target;
		this.trajectory = trajectory;
		this.thrust = thrust;
		Jf = (DecimalFormat)NumberFormat.getInstance();
		Jf.applyPattern(" ###.##;-###.##");
		Jf.setMinimumFractionDigits(2);
		Jf.setMinimumIntegerDigits(4);
		// Integers
		itf = NumberFormat.getInstance();
		itf.setMinimumIntegerDigits(6);
		
		// create an RungeKuttaFehlberg78 integrator
		rk = new RungeKutta8();
		//rk=new RungeKuttaFehlberg78();
		//rk.setNonAdaptive();
		rk.setStepSize(1.e2);
		//rk.setAdaptive();
		
		

	}

	// Performance index to be minimized is the penalty function |x_target - x_final|   
	public double evaluate(double[] x_search)
	{
		thrust.set_x_search(x_search);
		tf_CSI=x_search[0];
		double[] finalstate = rk.integrate(0., x0, tf_CSI, trajectory, thrust, true);
		VectorN target = new VectorN(x_target);
		VectorN actual = new VectorN(finalstate, 6);
		VectorN diff = actual.minus(target);
		double result;
		result = Math.abs(diff.x[0]/1000.) + Math.abs(diff.x[1]/1000.) +  Math.abs(diff.x[3]) +  Math.abs(diff.x[4]);
		//result =diff.mag();
		return result;
	}

	void create_grid()
	{
		double[] x_search = new double[2];
		double result;

		//		System.out.println(df.format(ux)+"  "+df.format(uy)+"  "+ " " + miss);
		for (x_search[1] = Math.toRadians(75.); x_search[1] < Math.toRadians(145.); x_search[1] += Math.toRadians(1.))
		{
			//System.out.print(Jf.format(x_search[0])+" ");
			System.out.print(" " + Jf.format(Math.toDegrees(x_search[1])));
			for (x_search[0] = 3900.; x_search[0] < 4100; x_search[0] += 20)
			{
				result = evaluate(x_search);
				System.out.print(" " + Jf.format(result));
			}
			System.out.println();
		}
	}

	/*
	static double compute_J(
		RungeKutta8 rk8,
		double[] x0c,
		double tf_CSI,
		double[] x_target,
		TwoBodyCSI trajectory,
		CSI csi,
		double ux,
		double uy)
	{
		csi.set_thrust(ux, uy, 0.);
		double[] finalstate = rk8.integrate(0., x0c, tf_CSI, trajectory, csi, true);
		VectorN target = new VectorN(x_target);
		VectorN arrive = new VectorN(finalstate, 6);
		VectorN diff = arrive.minus(target);
		//new VectorN(arrive).print("arrive");
		//new VectorN(target).print("target");
		//new VectorN(diff).print("diff");
		//System.out.println("final mass: " + finalstate[6]);
	
		// Find a matching final orbit
		//		TwoBody finalorbit = new TwoBody(cm.mu,finalstate);
		//		double vel=finalorbit.getV().mag();
		//		if(vel>5. &&vel<6.)
		//		{
		//			System.out.println();
		//			System.out.println("Thrust "+uf.format(ux)+"  "+uf.format(uy));
		//			System.out.println("Velocity "+finalorbit.getV().mag());
		//			finalorbit.getV().print("Velocity");
		//		}
		//		finalorbit.print("final orbit");
	
		double result;
		result = Math.abs(diff.x[0]) + Math.abs(diff.x[1]) + 1000. * Math.abs(diff.x[3]) + 1000. * Math.abs(diff.x[4]);
		//return arrive.minus(target).mag();
		return result;
	}
	*/

}
