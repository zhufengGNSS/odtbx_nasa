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
 */

package jat.test.alg.integrators;

import jat.alg.integrators.Derivatives;
import jat.alg.integrators.DormandPrince;
import jat.matvec.data.VectorN;

public class DPTest1 implements Derivatives
{

	public double[] derivs(double t, double[] x)
	{
		double[] out = new double[6];
		out[3] = -x[0];
		out[4] = -x[1];
		out[5] = -x[2];
		return out;
	}
	
	
	public static void main(String[] args)
	{
		System.out.println("Dormand Prince integrator test");
		DormandPrince dp = new DormandPrince();
		
		DPTest1 dpt1=new DPTest1();
		// double[] r1 = new double[4];
		// double[] r2 = new double[4];
		// double[] r1prime = new double[4];
		// double[] r2prime = new double[4];
		// r1prime[1] = 10.;
		// r1prime[2] = 0.;
		// double time1 = 0.;
		// double time2 = 2*Math.PI/12;
		// dp.update(r1, r1prime, time1, r2, r2prime, time2);
		// new VectorN(r2).print("r2");
		// new VectorN(r2prime).print("r2prime");
		// double[] r1 = new double[3];
		// double[] r2 = new double[3];
		// double[] r1prime = new double[3];
		// double[] r2prime = new double[3];
		// r1prime[0] = 10.;
		// r1prime[1] = 1.;
		// r1prime[2] = 1.;
		// double time1 = 0.;
		// double time2 = 2*Math.PI/12;
		// dp.integrate(r1, r1prime, time1, r2, r2prime, time2);
		// new VectorN(r2).print("r2");
		// new VectorN(r2prime).print("r2prime");
		double[] state1 = new double[6];
		state1[3] = 10.;
		state1[4] = 1.;
		state1[5] = 2.;
		double time1 = 0.;
		double time2 = 0 * Math.PI / 12;
		double[] state2 = dp.integrate(time1, state1, time2, dpt1);
		new VectorN(state2).print("state2");	}
}
