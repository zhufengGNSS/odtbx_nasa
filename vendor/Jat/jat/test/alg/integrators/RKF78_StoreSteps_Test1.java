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
 */

package jat.test.alg.integrators;

import jat.alg.integrators.Derivatives;
import jat.alg.integrators.RKF78_StoreSteps;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.matvec.io.data.fileTools.MatrixString;
import jat.util.PrintMatrix;

/**
 * Test the Runge Kutta integrator by comparing it's output of a harmonic oscillator with the sine and cosine functions
 * @author Tobias Berthold
 */
public class RKF78_StoreSteps_Test1 implements Derivatives
{
	public double[] derivs(double t, double[] x)
	{
		double[] out = new double[2];
		out[0] = x[1];
		out[1] = -x[0];
		return out;
	}

	public static void main(String[] args)
	{
		double time1 = 0.;
		double time2 = 24 * Math.PI / 12;
		double h1 = 0.01;
		double hmin = 0.0;
		double[] state2;
		RKF78_StoreSteps_Test1 rt = new RKF78_StoreSteps_Test1();
		RKF78_StoreSteps r = new RKF78_StoreSteps();
		//
		r.eps = 1.e-4;
		double[] state1 = new double[2];
		state1[0] = 0.;
		state1[1] = 1.;
		// Integrate
		state2 = r.integrate(time1, state1, time2, rt);
		new VectorN(state2).print("state");
		Matrix M = new Matrix(r.steplist);
		PrintMatrix PM = new PrintMatrix(M);
		String Titles[] =
		{ "t", "x", "x'", "retried steps" };
		PM.titles = Titles;
		PM.setMinFracDig(8);
		PM.print();
		System.out.println("good " + r.nok + " bad " + r.nbad + " steps");
		// N.toFile("outrk45.txt");
		// again to verify that the same integrator can be reused
		state2 = r.integrate(time1, state1, time2, rt);
		new VectorN(state2).print("state");
		Matrix N = new Matrix(r.steplist);
		PrintMatrix PN = new PrintMatrix(N);
		PN.setMinFracDig(8);
		PN.print();
		System.out.println("good " + r.nok + " bad " + r.nbad + " steps");
	}
}
