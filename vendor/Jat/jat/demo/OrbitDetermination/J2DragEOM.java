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
 * 
 * File Created on May 8, 2003
 */
package jat.demo.OrbitDetermination;

import jat.alg.integrators.Derivatives;
import jat.matvec.data.*;

/**
* The J2DragEOM.java Class provides the equations of motion (derivatives)
* for an orbit including J2 and the effects of exponential drag.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class J2DragEOM implements Derivatives {

	private static double re = 6378136.3; // radius of earth in meters
	private static double h_0 = 920000.0; // atmosphere model parameter
	private static double rho_0 = 4.36E-14; // atmosphere model parameter
	private static double gamma_0 = 5.381E-06; // atmosphere model parameter
	private static double omega_e = 7.2921157746E-05; // earth rotation rate

	/**
	 * Compute the time derivatives
	 * @param t time
	 * @param y double [] containing current state.
	 * @return double[] containing the time derivatives.
	 */
	public double[] derivs(double t, double[] y) {

		double out[] = new double[156];

		// store elements of incoming state in more familiar looking variables
		double xx = y[0];
		double yy = y[1];
		double zz = y[2];
		double vx = y[3];
		double vy = y[4];
		double vz = y[5];
		double mu = y[6];
		double j2 = y[7];
		double beta = y[8];

		// compute derived variables
		VectorN r = new VectorN(xx, yy, zz);
		double rmag = r.mag();
		double rcubed = rmag * rmag * rmag;
		double rsq = rmag * rmag;
		double re_r = re / rmag;
		re_r = re_r * re_r;
		double zsq_rsq = (5.0 * zz * zz / rsq) - 1.0;

		// compute accelerations due to gravity
		double ll = -1.0 * (mu * xx / rcubed) * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
		double mm = -1.0 * (mu * yy / rcubed) * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
		double nn = -1.0 * (mu * zz / rcubed) * (1.0 - 1.5 * re_r * j2 * (zsq_rsq - 2.0));

		// compute accelerations due to drag
		double r0 = re + h_0;
		double exp = -gamma_0 * (rmag - r0);
		double rho = rho_0 * Math.exp(exp);

		double aa = vx + omega_e * yy;
		double bb = vy - omega_e * xx;
		double cc = vz;
		double qq = Math.sqrt(aa * aa + bb * bb + cc * cc);
		double brq = beta * rho * qq;

		// compute state derivatives

		out[0] = y[3];
		out[1] = y[4];
		out[2] = y[5];
		out[3] = ll - brq * aa;
		out[4] = mm - brq * bb;
		out[5] = nn - brq * cc;

		for (int i = 6; i < 12; i++) {
			out[i] = 0.0;
		}

		// compute A matrix

		Matrix a = new Matrix(12, 12); // creates a 12x12 matrix with all zeros

		a.A[0][3] = 1.0;
		a.A[1][4] = 1.0;
		a.A[2][5] = 1.0;

		double r5 = rsq * rcubed;
		double mur5 = mu / r5;
		double mur3 = mu / rcubed;
		double sz2r2 = 7.0 * zz * zz / rsq;

		double muxyr5 = mu * xx * yy / r5;
		double muxzr5 = mu * xx * zz / r5;
		double muyzr5 = mu * yy * zz / r5;

		double bracket1 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 1.0);
		double bracket3 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 3.0);
		double bracket5 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 5.0);
		double bracket2 = 1.5 * re_r * (5.0 * zz * zz / rsq - 1.0);

		double dldx = ll / xx + mur5 * xx * xx * bracket1;
		double dldy = muxyr5 * bracket1;
		double dldz = muxzr5 * bracket3;
		double dldj2 = mur3 * xx * bracket2;

		double dmdx = dldy;
		double dmdy = mm / yy + mur5 * yy * yy * bracket1;
		double dmdz = muyzr5 * bracket3;
		double dmdj2 = mur3 * yy * bracket2;

		double dndx = muxzr5 * bracket3;
		double dndy = dmdz;
		double dndz = nn / zz + mur5 * zz * zz * bracket5;
		double dndj2 = mur3 * zz * (1.5 * re_r * (5.0 * zz * zz / rsq - 3.0));

		double bp = beta * rho;
		double aq = aa / qq;
		double bq = bb / qq;
		double cq = cc / qq;
		double qwe = qq * omega_e;
		double xr = xx / rmag;
		double yr = yy / rmag;
		double zr = zz / rmag;
		double awe = aa * omega_e;
		double bwe = bb * omega_e;

		double bpaq = bp * aq;
		double bpbq = bp * bq;
		double bpcq = bp * cq;

		double paren1 = bwe / qq + qq * gamma_0 * xr;
		double paren2 = qq * gamma_0 * yr - awe / qq;
		double paren3 = brq * gamma_0 * zr;
		double rhoq = rho * qq;

		a.A[3][0] = dldx + bp * aa * paren1;
		a.A[3][1] = dldy + bp * (qq * gamma_0 * aa * yr - qwe - awe * aq);
		a.A[3][2] = dldz + paren3 * aa;
		a.A[3][3] = -bp * (aa * aq + qq);
		a.A[3][4] = -bpbq * aa;
		a.A[3][5] = -bpcq * aa;
		a.A[3][6] = ll / mu;
		a.A[3][7] = dldj2;
		a.A[3][8] = -rhoq * aa;

		a.A[4][0] = dmdx + bp * (qwe + bwe * bq - qq * gamma_0 * bb * xr);
		a.A[4][1] = dmdy + bp * bb * paren2;
		a.A[4][2] = dmdz - paren3 * bb;
		a.A[4][3] = -bpaq * bb;
		a.A[4][4] = -bp * (bb * bq + qq);
		a.A[4][5] = -bpcq * bb;
		a.A[4][6] = mm / mu;
		a.A[4][7] = dmdj2;
		a.A[4][8] = -rhoq * bb;

		a.A[5][0] = dndx + bp * cc * paren1;
		a.A[5][1] = dndy + bp * cc * paren2;
		a.A[5][2] = dndz + paren3 * cc;
		a.A[5][3] = -bpaq * cc;
		a.A[5][4] = -bpbq * cc;
		a.A[5][5] = -bp * (cc * cq + qq);
		a.A[5][6] = nn / mu;
		a.A[5][7] = dndj2;
		a.A[5][8] = -rhoq * cc;

		//	   a.print("A matrix");

		// compute phi derivatives

		Matrix phi = this.phi(y);

		Matrix dphi = a.times(phi);

		//	   dphi.print("dphi");

		// put phi derivatives into output array
		int k = 12;
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 12; j++) {
				out[k] = dphi.A[i][j];
				k = k + 1;
			}
		}

		return out;
	}

	private Matrix phi(double[] in) {
		Matrix out = new Matrix(12, 12);
		int k = 12;
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 12; j++) {
				out.A[i][j] = in[k];
				k = k + 1;
			}
		}
		return out;
	}

}
