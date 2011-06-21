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

package jat.alg.opt.test;

import jat.alg.*;
import jat.alg.opt.*;

public class InnerClassExample
{
	public class Function implements ScalarfromArrayFunction
	{
		public double evaluate(double[] x)
		{
			return (0.1 * x[0] * x[0] - 2. * x[0] + 10.);
		}
	}

	static void printmin(double[] x)
	{
		System.out.print("The min is at  ");
		for (int i = 0; i < x.length; i++)
			System.out.print("x" + i + "= " + x[i] + "  ");

		System.out.println();
		System.out.println();
	}

	public static void main(String argv[])
	{
		double[] x; // result
		double[] x1 = new double[1]; // initial guess
		double[] d1 = new double[1]; // search direction

		//System.out.println("Unconstrained Parameter Optimization Test");

		// create instances of the classes
		InnerClassExample ot = new InnerClassExample();

		x1[0] = 0.;
		d1[0] = 1.;
		//opt.G=ot.new Function();
		x = LineSearch.ods(ot.new Function(), x1, d1, 1.e-4);
		printmin(x);

	}
}
