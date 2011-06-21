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

import java.text.*;

import jat.alg.ScalarfromArrayFunction;
import jat.cm.*;
import jat.matvec.data.VectorN;

/**
 * @author Tobias Berthold
 *
 */
public class PerformanceIndex implements ScalarfromArrayFunction
{
	Lambert earth_mars;
	DecimalFormat Jf; // formatting  
	NumberFormat itf;
	double tf_var;
	double alow, ahigh, astep;

	public PerformanceIndex( Lambert earth_mars)
	{
		this.earth_mars = earth_mars;
		Jf = (DecimalFormat)NumberFormat.getInstance();
		Jf.applyPattern(" ###.##;-###.##");
		Jf.setMinimumFractionDigits(2);
		Jf.setMinimumIntegerDigits(4);
		// Integers
		itf = NumberFormat.getInstance();
		itf.setMinimumIntegerDigits(6);	
	}

	// Performance index to be minimized is the penalty function: delta-v   
	public double evaluate(double[] x_search)
	{
		// compute parameters for current guess
		double muforthisproblem = Constants.GM_Sun/1.e9;
		TwoBody elem0 = new TwoBody(muforthisproblem, cm.earth_moon_elements );
		TwoBody elemf = new TwoBody(muforthisproblem,cm.mars_elements);
		elem0.setTa(x_search[1]);
		elemf.setTa(x_search[2]);
		
		VectorN r0 = elem0.getR();
		VectorN v0 = elem0.getV();
		VectorN rf = elemf.getR();
		VectorN vf = elemf.getV();
		double dt=x_search[0];
		
		return earth_mars.compute(r0, v0, rf, vf, dt);
	}

	public void initgrid(double d, double e, double f)
	{
		alow=d;
		ahigh=e; 
		astep=f;		
	}
	
	void print_header()
	{
		
		System.out.print("x1 / x0 ");
		for (double x = alow; x < ahigh; x += astep)
		{
			System.out.print(" " + Jf.format(x));
		}
		System.out.println();
	}
	
	void create_grid()
	{
		double[] x_search = new double[3];
		double result;

		x_search[2]=120.;
		for (x_search[1] = -100.; x_search[1] <= -80; x_search[1] += 10.)
		{
			//System.out.print(Jf.format(x_search[0])+" ");
			System.out.print(" " + Jf.format(x_search[1]));
			for (x_search[0] = alow; x_search[0] < ahigh; x_search[0] += astep)
			{
				result = evaluate(x_search);
				System.out.print(" " + Jf.format(result));
			}
			System.out.println();
		}
	}
	
	void create_combinations()
	{
		double a,b;
		double[] x_search = new double[3];
		double result=0.;

		//for (x_search[1] = -99.; x_search[1] < -80; x_search[1] += 1.)
		for (b = -90.; b < -70; b += 10.)
		{
			//for (x_search[0] = alow; x_search[0] < ahigh; x_search[0] += astep)
			for (a = alow; a < ahigh;  a+= astep)
			{
				x_search[0]=a;
				x_search[1]=b;
				x_search[2]=130.;
				System.out.print("" + Jf.format(x_search[0]));
				System.out.print(" " + Jf.format(x_search[1]));
				System.out.print(" " + Jf.format(x_search[2]));
				result = evaluate(x_search);
				System.out.println(" " + Jf.format(result));
			}
		}
	}
	
}
