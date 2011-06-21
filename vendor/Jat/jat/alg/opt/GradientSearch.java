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
 *
 */

package jat.alg.opt;

import jat.alg.*;

/**
 * @author Tobias Berthold
 * Date        :   6-21-2004
 * @version 1.0
 */

public class GradientSearch extends optimize
{
	public double err_ods = 1.e-4; // error tolerance for linesearch
	public double err_grad = 1.e-6; // error tolerance for gradient search
	public double eps_FD = 1.e-8; // Perturbation for forward difference 
	public int max_it = 50; // maximum iterations 

	public GradientSearch(ScalarfromArrayFunction G, double[] x_init)
	{
		super(G, x_init);
	}

	/**  Set up the optimize class for finding the minimum of a function 
	 * by gradient search method
	 * @param G The function the minimum of which is to be found
	 * @param x_init Initial guess vector
	 * @param err_ods Error tolerance for linesearch
	 */
	public GradientSearch(ScalarfromArrayFunction G, double[] x_init, double err_ods)
	{
		super(G, x_init, err_ods);
	}

	/**  Set up the optimize class for finding the minimum of a function 
	 * by gradient search method
	 * @param G The function the minimum of which is to be found
	 * @param x_init Initial guess vector
	 * @param err_ods Error tolerance for linesearch
	 * @param max_it maximum iterations
	 */
	public GradientSearch(ScalarfromArrayFunction G, double[] x_init, double err_ods, int max_it)
	{
		super(G, x_init, err_ods);
		this.max_it = max_it;
		this.err_ods = err_ods;
	}

	/** Find the minimum of a function by gradient search method
	 * @param G The function the minimum of which is to be found
	 * @param x_init Initial guess vector
	 * @param epss Error tolerance for linesearch
	 * @param max_it maximum iterations
	 * @param eps_FD Perturbation factor for finite difference
	 * @param err_grad error tolerance for gradient search
	 */
	public GradientSearch(ScalarfromArrayFunction G, double[] x_init, double err_ods, int max_it, double eps_FD, double err_grad)
	{
		super(G, x_init, err_ods);

		this.max_it = max_it;
		this.eps_FD = eps_FD;
		this.err_grad = err_grad;
	}

	//	public void set_initial_guess(double[] x)
	//	{
	//		//        for(int i=0;i<x.length;i++)
	//		//            this.x[i]=x[i];
	//		n = x.length - 1;
	//	}

	private void copy(double[] from, double[] to)
	{
		int l = from.length;
		System.arraycopy(from, 0, to, 0, l);
	}

	/** Find the minimum by getting the search direction with the gradient method
	 * 
	 */
	public double[] find_min_gradient()
	{
		int status=0;
		int i, it = 1;
		double norm = 0.;
		boolean more_iter = true;

		// Copy initial guess to x
		copy(x_init, this.x);
		print_header();
		while (more_iter)
		{
			d = NumDerivs.G_x_forward(G, x, eps_FD);
			norm = norm(d);
			print_line(it, x, G.evaluate(x), d, norm);
			if (norm < err_grad)
			{
				more_iter = false;
			} else
			{
				for (i = 0; i < n; i++)
					d[i] = -d[i];
				LineSearch.ods(G, x, d, err_ods);
				it++;
				if (it > max_it)
				{
					more_iter = false;
					status=1;
				}
				if (LineSearch.status >0)
				{
					System.out.println("Linesearch failed, status: "+LineSearch.status);				
					// x might still be minimum, check
					more_iter = false;			
					status=2;
				}
			}
		}
		
		// Conclusion		
		if(status==0)
				System.out.println("Convergence:");
		if(status==1)
				System.out.println("Maximum number of iterations reached");
		if(status==2)
				System.out.println("Linesearch failed");
		for (i = 0; i < x.length; i++)
			System.out.print("x" + i + "= " + x[i] + "  ");
		System.out.println("");
		System.out.println("|Gx|= "+norm);
		
		return x;
	}
}
