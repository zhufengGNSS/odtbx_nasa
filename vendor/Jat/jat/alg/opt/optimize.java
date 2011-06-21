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
import java.text.*;

/**
 * Base class for unconstrained parameter optimization with 
 * iterative methods: Gradient search, DFP, etc. 
 * 
 *
 * @author Tobias Berthold
 * Date        :   6-10-2002
 */

public class optimize
{
	public int n; // Dimension of search vector
	public double[] x_init; // initial guess
	public double[] d_init; // initial search direction    
	public double[] x; // solution
	public double[] d; // search direction
	DecimalFormat df;
	NumberFormat itf;
	ScalarfromArrayFunction G; // 


	public optimize(ScalarfromArrayFunction G, double[] x_init)
	{
		this.G = G;
		this.x_init = x_init;
		set_print_formats();
		//
		n = x_init.length;
		x = new double[n];
		d = new double[n];
		array_copy(x_init, this.x);
	}

	/**
	 * @param G
	 * @param x_init
	 * @param err_ods
	 */
	public optimize(ScalarfromArrayFunction G, double[] x_init, double err_ods)
	{
		this.G = G;
		this.x_init = x_init;
		//set_print_formats();
		//
		n = x_init.length;
		x = new double[n];
		d = new double[n];
		array_copy(x_init, this.x);
	}
	
	private void array_copy(double[] from, double[] to)
	{
		int l = from.length;
		System.arraycopy(from, 0, to, 0, l);
	}
	
	public void set_print_formats()
	{
		df = (DecimalFormat)NumberFormat.getInstance();
		df.applyPattern("  ###.########; -###.#######");
		//df.applyPattern("####0.000000; -###0.000000");
		df.setMinimumFractionDigits(8);
		df.setMinimumIntegerDigits(4);
		itf = NumberFormat.getInstance();
		itf.setMinimumIntegerDigits(4);
	}

	void print_header()
	{
		int i;
		String sp = "                ";

		System.out.print(" It  ");
		for (i = 0; i < x.length; i++)
			System.out.print("|      x("+i+")    ");
		System.out.print("|      G       ");
		for (i = 0; i < x.length; i++)
			System.out.print("|    Gx("+i+")     ");
		System.out.print("       |Gx|");
		System.out.println();
	}

	void print_line(int it, double[] x, double G, double[] d, double norm)
	{
		int i;
		
		System.out.print(itf.format(it));
		for (i = 0; i < x.length; i++)
			System.out.print(df.format(x[i]));
		System.out.print(df.format(G));
		for (i = 0; i < x.length; i++)
			System.out.print(df.format(-d[i]));
		System.out.print(df.format(norm));
		System.out.println();
	}

	double norm(double[] vec) // Vector norm
	{
		int i;
		double norm = 0;
		for (i = 0; i < n; i++)
			norm += vec[i] * vec[i];
		return Math.sqrt(norm);
	}

}
