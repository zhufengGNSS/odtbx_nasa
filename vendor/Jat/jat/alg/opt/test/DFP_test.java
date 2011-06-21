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
 */

package jat.alg.opt.test;

import jat.alg.opt.*;
import jat.alg.opt.test.functions.*;

/**
 * Davidon-Fletcher-Powell variable metric method
 * @author Tobias Berthold
 *
 */
public class DFP_test
{

	public static void main(String argv[])
	{
		double[] x_init=new double[2];
		
		System.out.println("Rosenbrock function, Numerical derivs, DFP");

		// create instances of the classes
		DFP_test dt = new DFP_test();
		x_init[0] = -1.2;
		x_init[1] = 1.;
		DFP dfp = new DFP(new Rosenbrock(), x_init);
		dfp.err_ods=1.e-6;
		dfp.err_dfp=1.e-6;
		dfp.eps_CD=1.e-5;
		dfp.max_it=50;
		double[] x=dfp.find_min_DFP();
		
	}
}
