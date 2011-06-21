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

import jat.alg.ScalarfromArrayFunction;
import jat.matvec.data.*;

/**
 * @author Tobias Berthold
 *
 */
public class NumDerivs
{

	public NumDerivs()
	{
	}

	// Forward difference numerical derivative
	static double[] G_x_forward(ScalarfromArrayFunction G, double[] x,double eps_FD) 
	{
		int i;
		int n=x.length;
		double[] gx = new double[n];
		double[] xf = new double[n];
		double dx,GM,GF;

		//double G1 = G.evaluate(x); // used several times, compute once
		GM = G.evaluate(x); // used several times, compute once
		copy(x, xf);
		for (i = 0; i < n; i++)
		{
			// Perturbation step
			if (Math.abs(x[i]) < 1.)
				dx = eps_FD;
			else
				dx = Math.abs(eps_FD * x[i]);
			// Finite difference
			xf[i] = x[i] + dx;
			GF=G.evaluate(xf);
			gx[i] = ( GF- GM) / dx;
			xf[i] = x[i];
		}
		return gx;
	}

	/** Central difference numerical derivative
	 * @param G
	 * @param x
	 * @param eps_CD
	 * @return
	 */
	static double[] G_x_central(ScalarfromArrayFunction G, double[] x,double eps_CD) 
	{
		int i;
		int n=x.length;
		double GB, GF, dx;
		double[] gx = new double[n];
		double[] xc = new double[n];

		copy(x, xc);
		for (i = 0; i < n; i++)
		{
			// Perturbation step
			if (Math.abs(x[i]) < 1.)
				dx = eps_CD;
			else
				dx = Math.abs(eps_CD * x[i]);
			// Central difference
			xc[i] = x[i] + dx;
			GF = (G.evaluate(xc));
			xc[i] = x[i] - dx;
			GB = (G.evaluate(xc));
			xc[i] = x[i];
			gx[i] = 0.5 * (GF - GB) / dx;
		}
		return gx;
	}

	private static void copy(double[] from, double[] to)
	{
		//int l = from.length;
		System.arraycopy(from, 0, to, 0, from.length);
	}

	private void copy(VectorN from, VectorN to)
	{
		System.arraycopy(from.getArray(), 0, to.getArray(), 0, from.length);
	}
}
