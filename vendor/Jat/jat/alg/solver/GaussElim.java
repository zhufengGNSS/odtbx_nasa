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

package jat.alg.solver;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

class GaussElim
{
	public static Matrix C; // Augmented matrix  C= [ A | y ]

	static double DABS(double d)
	{
		return Math.abs(d);
	}

	/**
	 *  Solve the linear system y = Ax by Gaussian Elimination
	 * as described in "Numerical Analysis" by Richard Burden and J. Douglas Faires
	 * @param A Matrix
	 * @param y vector
	 * @return solution x
	 *
	 * @author Tobias Berthold
	 * Date        :   6-16-2004
	 * @version 1.0
	 */
	static VectorN solve(Matrix A, VectorN y)
	{
		//		double R, S, Z, V, EX, AER, RATIO, sum;
		//		int J, K, L, KP, LP1, i, j;
		int rows, cols;
		boolean OK = true;
		int ICHG, M, I, IP, NN, K, KK, J, JJ;
		double XM, SUM, temp;
		double ZERO = 1.0E-20;
		int ip;

		rows = A.m;
		cols = A.n;
		//		int n = rows - 1;
		int N = A.n;

		C = new Matrix(rows, cols + 1);
		C.setMatrix(0, 0, A);
		C.setColumn(cols, y);
		//C.print("Augmented matrix C:");
		VectorN X = new VectorN(rows);

		/* STEP 1 */
		NN = N - 1;
		M = N + 1;
		ICHG = 0;
		I = 0;
		while ((OK) && (I < NN))
		{
			// STEP 2 : Find pivot 
			/* use IP instead of p */
			IP = I;
			while ((DABS(C.A[IP][I]) <= ZERO) && (IP <= N))
				IP++;
			if (IP == M)
				OK = false;
			else
			{
				// STEP 3 : Interchange rows
				if (IP != I)
				{
					for (JJ = 1; JJ < M; JJ++)
					{
						temp = C.A[I][JJ];
						C.A[I][JJ] = C.A[IP][JJ];
						C.A[IP][JJ] = temp;
					}
					ICHG = ICHG + 1;
				}
				/* STEP 4 */
				JJ = I + 1;
				for (J = JJ; J < N; J++)
				{
					/* STEP 5 */
					/* use XM in place of m(J,I) */
					XM = C.A[J][I] / C.A[I][I];
					/* STEP 6 */
					for (K = JJ; K < M; K++)
						C.A[J][K] = C.A[J][K] - XM * C.A[I][K];
					/*  Multiplier XM could be saved in C.A[J,I].  */
					C.A[J][I] = 0.0;
				}
			}
			I = I + 1;
		}
		//C.print("C before back subs in solve");

		if (OK)
		{
			/* STEP 7 */
			if (DABS(C.A[N - 1][N - 1]) <= ZERO)
				OK = false;
			else
			{
				/* STEP 8 */
				/* start backward substitution */
				X.x[N - 1] = C.A[N - 1][M - 1] / C.A[N - 1][N - 1];
				/* STEP 9 */
				for (K = 1; K <= NN; K++)
				{
					I = NN - K + 1;
					JJ = I + 1;
					SUM = 0.0;
					for (KK = JJ; KK <= N; KK++)
						SUM = SUM - C.A[I - 1][KK - 1] * X.x[KK - 1];
					X.x[I - 1] = (C.A[I - 1][M - 1] + SUM) / C.A[I - 1][I - 1];
				}
				/* STEP 10 */
				/* procedure completed successfully */
				//OUTPUT(N, M, ICHG, X, A);
			}
		}
		if (!OK)
			System.out.println("System has no unique solution\n");

		return X;
	}

}
