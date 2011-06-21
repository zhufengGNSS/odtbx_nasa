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
import jat.matvec.data.VectorN;

/**
 * Davidon-Fletcher-Powell variable metric method
 * @author Tobias Berthold
 *
 */
public class DFP extends optimize
{
	public double err_ods = 1.e-4; // Error tolerance for linesearch
	public double err_dfp = 1.e-6; // Error tolerance for search
	public double eps_CD = 1.e-4; // Perturbation for central difference 
	public int max_it = 50; // maximum iterations 

	public DFP(ScalarfromArrayFunction G, double[] x_init)
	{
		super(G, x_init);
	}

	void DFP_update(VectorN dx, VectorN dg, Matrix H)
	{
		int n = dx.length;
		VectorN H_dg = new VectorN(n);
		double dxT_dg = 0., dgT_H_dg = 0.;
		Matrix dx_dxT, dH, dH1, dH2, H_dg_dgT_H;

		H_dg = H.times(dg);
		dxT_dg = dx.dotProduct(dg);
		dgT_H_dg = dg.dotProduct(H_dg);
		dx_dxT = dx.outerProduct(dx);
		H_dg_dgT_H = H_dg.outerProduct(H_dg);

		dH1 = dx_dxT.ebeDivide(dxT_dg);
		dH2 = H_dg_dgT_H.ebeDivide(dgT_H_dg);
		dH = dH1.minus(dH2);
		H.setMatrix(0, 0, H.plus(dH));
		return;
	}

	public double[] find_min_DFP()
	{
		Matrix H = new Matrix(n); // Set H to identity matrix
		VectorN x, xn, dx, gx, gxn, dgx;
		double[] dummy;
		int i, it = 0;
		double norm = 0.;
		boolean more_iter = true;
		int status = 0;

		// Copy initial guess to x
		x = new VectorN(x_init);
		xn = new VectorN(x_init);
		gxn = new VectorN(x_init);
		//copy(x_init, this.x);
		print_header();

		gx = new VectorN(NumDerivs.G_x_central(G, x.getArray(), eps_CD));
		while (more_iter)
		{
			norm = norm(gx.getArray());
			if (norm < err_dfp)
			{
				more_iter = false;
				print_line(it, x.getArray(), G.evaluate(x.getArray()), gx.getArray(), norm);
			} else
			{
				// Step 4
				dx = H.times(gx);
				dx = dx.times(-1.);
				//H.print("H");
				copy(x.getArray(), xn.getArray());
				dummy = LineSearch.ods(G, xn.getArray(), dx.getArray(), err_ods);
				print_line(it, x.getArray(), G.evaluate(x.getArray()), gx.getArray(), norm);
				// Step 5
				dx = xn.minus(x);
				// Step 6
				gxn = new VectorN(NumDerivs.G_x_central(G, xn.getArray(), eps_CD));
				dgx = gxn.minus(gx);
				// Step 7
				DFP_update(dx, dgx, H);
				// Step 8
				copy(xn, x);
				copy(gxn, gx);
				it++;
				if (it > max_it)
				{
					more_iter = false;
					status = 1;
				}
				if (LineSearch.status > 0)
				{
					System.out.println("Linesearch failed, status: " + LineSearch.status);
					// x might still be minimum, check
					more_iter = false;
					status = 2;
				}
			}
		}

		// Conclusion		
		if (status == 0)
			System.out.println("Convergence:");
		if (status == 1)
			System.out.println("Maximum number of iterations reached");
		if (status == 2)
			System.out.println("Linesearch failed");
		for (i = 0; i < x.length; i++)
			System.out.print("x" + i + "= " + x.x[i] + "  ");
		System.out.println("");
		System.out.println("|Gx|= " + norm);

		return x.getArray();

	}

	private void copy(double[] from, double[] to)
	{
		System.arraycopy(from, 0, to, 0, from.length);
	}

	private void copy(VectorN from, VectorN to)
	{
		System.arraycopy(from.getArray(), 0, to.getArray(), 0, from.length);
	}

}
