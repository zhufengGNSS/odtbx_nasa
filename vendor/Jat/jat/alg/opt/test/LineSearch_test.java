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

import java.util.Date;
import jat.alg.opt.*;
import jat.alg.opt.test.functions.*;

/**
 * Test the line search with a given search direction.
 * @author Tobias Berthold
 */
public class LineSearch_test
{
	private static void copy(double[] from, double[] to)
	{
		int l = from.length;
		System.arraycopy(from, 0, to, 0, l);
	}

	static void printmin(double[] xg, double[] d, double[] x)
	{
		System.out.print("Initial Guess:   ");
		for (int i = 0; i < xg.length; i++)
			System.out.print("xg" + i + "= " + xg[i] + "  ");
		System.out.println();
		System.out.print("Search Direction: ");
		for (int i = 0; i < d.length; i++)
			System.out.print("d" + i + "= " + d[i] + "  ");
		System.out.println();
		System.out.print("The min is at:    ");
		for (int i = 0; i < x.length; i++)
			System.out.print("x" + i + "= " + x[i] + "  ");
		System.out.println();
		System.out.println();
	}

	static void example1()
	{
		// Function 1
		System.out.println("Function 1");
		double[] xg = { 0. };
		double[] d = { 1. };
		double[] x = new double[1];
		copy(xg, x);

		Function1 G = new Function1();
		LineSearch.ods(G, x, d, 1.e-4);
		printmin(xg, d, x);
	}

	static void example2()
	{
		// Function 2
		double[] xg = { 0. };
		double[] d = { 1. };
		double[] x = new double[1];
		copy(xg, x);

		System.out.println("Function 2, tolerance = 1.e-2");
		LineSearch.ods(new Function2(), x, d, 1.e-2);
		printmin(xg, d, x);

		System.out.println("Function 2, tolerance = 1.e-4");
		copy(xg, x);
		LineSearch.ods(new Function2(), x, d, 1.e-4);
		printmin(xg, d, x);

		System.out.println("Function 2, tolerance = 1.e-8");
		copy(xg, x);
		LineSearch.ods(new Function2(), x, d, 1.e-8);
		printmin(xg, d, x);
	}

	static void example3()
	{
		// Function 3
		double[] xg = { 1., 1. };
		double[] d = { 2., 0.5 };
		double[] x = new double[2];
		copy(xg, x);

		System.out.println("Function 3");
		LineSearch.ods(new Function3(), x, d, 1.e-6);
		printmin(xg, d, x);

		d[0] = 1.;
		d[1] = 2.;
		copy(xg, x);
		LineSearch.ods(new Function3(), x, d, 1.e-6);
		printmin(xg, d, x);
	}

	public static void main(String argv[])
	{
		double epss;

		//System.out.println("Unconstrained Parameter Optimization Test");

		Date start_time = new Date();
		
		example1();
		example2();
		example3();
		
		Date stop_time = new Date();
		double etime = (stop_time.getTime() - start_time.getTime()) / 1000.;
		System.out.println("Elapsed Time = " + etime + " seconds");

	}
}
