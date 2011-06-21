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

package jat.alg.opt;

import jat.alg.*;

/**
 * Line Search
 * 
 * typically used in iterative method, therefore the parameter x (search vector) is
 * being modified in the function for better efficiency
 * 
 * @author Tobias Berthold
 *
 */
public class LineSearch
{
	static public int n; // Dimension of search vector
	static public double err_ods = 1.e-6d; // error tolerance for linesearch
	static public double[] xnew; //
	static ScalarfromArrayFunction G; //
	static int status;

	/*
	double G(double[] x)
	{
		return G.evaluate(x);
	}
	*/

	static double get_new_G(ScalarfromArrayFunction G, double[] x, double alpha, double[] d)
	{
		int i;

		for (i = 0; i < n; i++)
			xnew[i] = x[i] + alpha * d[i];

		return G.evaluate(xnew);
	}

	// CAUTION: the parameter x is changed in the function!!!
	/** Performs a linesearch (one-dimensional search) on the function G(x) in the
	 * search direction d.
	 * @return x : found minimum
	 */
	public static double[] ods(ScalarfromArrayFunction G, double[] x, double[] d, double err_ods)
	{
		// The numbers alpha1, alpha2, alpha3, and alpha4 determine points on the
		// line of search direction such that xnew = x + alpha * d .
		// The minimum is found when alpha4 is determined, and
		// xmin = x + alpha4 * d
		double a1 = 0., a2 = 1., a3 = 0., a4 = 0.;
		double c1 = 5.0;
		double g1, g2, g3;
		double g4, g4_bar, ge;
		boolean more_iter = true;
		int it = 0;

		status = 0;
		n = x.length;
		xnew = new double[n];

		g1 = get_new_G(G, x, a1, d);
		g2 = get_new_G(G, x, a2, d);
		g3 = Double.NEGATIVE_INFINITY;

		// Bracket the minimium: Find an up-down-up sequence in a1, a2, and a3
		while (g2 >= g1)
		{
			a3 = a2;
			g3 = g2;
			a2 = a3 / c1;
			if (Math.abs(a2) < 1.e-12)
			{
				//System.out.println("alpha3 too small");
				System.out.println("Function always increases in search direction");
				status = 1;
				more_iter = false;
				break;
			}
			g2 = get_new_G(G, x, a2, d);
		}
		while (g3 <= g2 && more_iter)
		{
			a3 = a1 + c1 * (a2 - a1);
			if (Math.abs(a3) > 1.e6)
			{
				//System.out.println("alpha3 too large");
				System.out.println("Function always decreases in search direction");
				status = 2;
				more_iter = false;
				break;
			}
			g3 = get_new_G(G, x, a3, d);
			if (g3 < g2)
			{
				a1 = a2;
				g1 = g2;
				a2 = a3;
				g2 = g3;
			}
		}

		// Minimum is bracketed, quad fit
		while (more_iter)
		{
			// Quadratic Fit
			double P = ((g3 - g1) * (a2 - a1) - (g2 - g1) * (a3 - a1)) / ((a3 - a1) * (a3 - a1) * (a2 - a1) - (a2 - a1) * (a2 - a1) * (a3 - a1));
			double Q = ((g2 - g1) - (a2 * a2 - a1 * a1) * P) / (a2 - a1);
			double R = g1 - a1 * a1 * P - a1 * Q;
			a4 = -Q / (2. * P);
			g4_bar = P * a4 * a4 + Q * a4 + R;
			g4 = get_new_G(G, x, a4, d);

			//System.out.println(g4 - g4_bar);

			if (Math.abs(g4) < err_ods)
				//ge = err_ods;
				ge = 1.;
			else
				ge = g4;

			if (Math.abs((g4 - g4_bar) / ge) <= err_ods)
				more_iter = false;
			if (more_iter)
			{
				//g4 = err_ods + 1.;
				//g4_bar = 0.;
				//ge = 1.;

				// find the 3 new bracketting values among the 4
				if (a4 > a2)
				{
					if (g4 > g2)
					{
						g3 = g4;
						a3 = a4;
					} else
					{
						g1 = g2;
						g2 = g4;
						a1 = a2;
						a2 = a4;
					}
				} else
				{
					if (g4 > g2)
					{
						g1 = g4;
						a1 = a4;
					} else
					{
						g3 = g2;
						g2 = g4;
						a3 = a2;
						a2 = a4;
					}
				}
			}
			it++;
			if (it > 200)
			{
				more_iter = false;
				status=3;
				System.out.println("Too many iterations in quadratic fit");
			}
		}
		// Enter minimum into x
		int i;
		for (i = 0; i < n; i++)
			x[i] = x[i] + a4 * d[i];
		return x;
	}
}

/*
// compute temporary result for printing
double[] xtmp = new double[n];
int i;
for (i = 0; i < n; i++)
	xtmp[i] = x[i] + a4 * d[i];
//			System.out.println(xtmp[0]+" "+xtmp[1]);
*/
