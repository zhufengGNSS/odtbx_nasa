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
import java.text.*;

/**
 * @author Tobias Berthold
 *
 */
public class VarMetric extends optimize
{
	public double err_ods = 1.e-4; // Error tolerance for linesearch
	public double err_vm = 1.e-6; // Error tolerance for search
	public double eps_CD = 1.e-4; // Perturbation for central difference
	public int max_it = 50; // maximum iterations
	public double gamma = 1.;
	DecimalFormat df;
	NumberFormat itf;
	boolean print_iter = false;
	public int total_it;
	public double norm = 0.;
	public boolean ods_fail = false;

	public VarMetric(ScalarfromArrayFunction G, double[] x_init)
	{
		super(G, x_init);
	}

	void H_update(double gamma, VectorN dx, VectorN dg, Matrix H)
	{
		int n = dx.length;
		VectorN H_dg = new VectorN(n);
		double dxT_dg = 0., dgT_H_dg = 0.;
		double fac1;
		Matrix DFP, D1, D2, BFSG, B1, B2, B3, VM, VM1;
		Matrix dx_dxT, H_dg_dgT_H, dx_dgT, dx_dgT_H, H_dg_dxT;
		Matrix I = Matrix.identity(n, n);

		// Scalars
		H_dg = H.times(dg);
		dxT_dg = dx.dotProduct(dg);
		dgT_H_dg = dg.dotProduct(H_dg);
		fac1 = gamma * dgT_H_dg / dxT_dg;

		// Matrices
		dx_dxT = dx.outerProduct(dx);
		dx_dgT = dx.outerProduct(dg);
		H_dg_dgT_H = H_dg.outerProduct(H_dg);
		dx_dgT_H = dx_dgT.times(H);
		H_dg_dxT = dx_dgT_H.transpose();

		// Matrix terms
		D1 = dx_dxT.ebeDivide(dxT_dg);
		D2 = H_dg_dgT_H.ebeDivide(dgT_H_dg);
		B1 = dx_dgT.ebeDivide(dxT_dg);
		B2 = I.minus(B1);
		B3 = B2.transpose();

		// Update formulas
		DFP = D1.minus(D2).plus(H);
		BFSG = B2.times(H).times(B3).plus(D1);
		VM1 = DFP.times(fac1);
		VM = BFSG.plus(VM1).ebeDivide(1. + fac1);
		//H.setMatrix(0, 0, DFP);
		//H.setMatrix(0, 0, BFSG);
		H.setMatrix(0, 0, VM);
		return;
	}

	public double[] find_min_VarMetric()
	{
		Matrix H = new Matrix(n); // Set H to identity matrix
		VectorN x, xn, dx, gx, gxn, dgx;
		double[] dummy;

		int i, it = 0;
		boolean more_iter = true;
		ods_fail = false;

		// Copy initial guess to x
		x = new VectorN(x_init);
		xn = new VectorN(x_init);
		gxn = new VectorN(x_init);
		//copy(x_init, this.x);
		if (print_iter)
			print_header();

		gx = new VectorN(NumDerivs.G_x_central(G, x.getArray(), eps_CD));
		while (more_iter)
		{
			norm = norm(gx.getArray());
			if (norm < err_vm)
			{
				more_iter = false;
				if (print_iter)
					print_line(it, x.getArray(), G.evaluate(x.getArray()), gx.getArray(), norm);
			} else
			{
				// Step 4
				dx = H.times(gx);
				dx = dx.times(-1.);
				//H.print("H");
				copy(x.getArray(), xn.getArray());
				dummy = LineSearch.ods(G, xn.getArray(), dx.getArray(), err_ods);
				if (print_iter)
					print_line(it, x.getArray(), G.evaluate(x.getArray()), gx.getArray(), norm);
				// Step 5
				dx = xn.minus(x);
				// Step 6
				gxn = new VectorN(NumDerivs.G_x_central(G, xn.getArray(), eps_CD));
				dgx = gxn.minus(gx);
				// Step 7
				//H_update(-.302301, dx, dgx, H);
				H_update(gamma, dx, dgx, H);
				// Step 8
				copy(xn, x);
				copy(gxn, gx);
				it++;
				if (it > max_it)
					more_iter = false;
				if (LineSearch.status == 1)
				{
					// x might still be minimum, check
					more_iter = false;
					ods_fail = true;
				}
			}
		}
		total_it = it;
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
//H.setMatrix(0, 0, H.plus(DFP));
