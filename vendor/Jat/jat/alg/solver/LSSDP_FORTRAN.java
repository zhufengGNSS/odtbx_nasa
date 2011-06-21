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

package jat.alg.solver;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

class LSSDP_FORTRAN
{

	static double DABS(double d)
	{
		return Math.abs(d);
	}

	// Solve y=Ax by Gaussian Elimination
	// C= [A | y]

	static void LSSDP_solve(Matrix C, VectorN X)
	{
		// IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		//     SOLUTION OF A LINEAR SYSTEM BY GAUSS ELIMINATION
		// NR=number of rows
		// NPM=number of columns
		//  DIMENSION C(NR,NPM),X(NR)
		//  EQUIVALENCE (R,S),(RATIO,V)
		double R, S, Z, V, EX, AER, RATIO;
		int I, II, IIP1, J, K, L, KP, LP1;
		int NR, NPM;

		System.out.println("" + C.m);
		NR = C.m - 1;
		NPM = C.n - 1;
		// NORM OF MATRIX C
		R = 0.0;
		for (J = 1; J <= NR; J++)
		{
			Z = 0.0;
			for (K = 1; K <= NR; K++)
			{
				V = C.A[J][K];
				Z = Z + DABS(V);
			}
			if (R < Z)
				R = Z;
		}
		EX = (1.0e-30) * R;
		// TRIANGULARIZATION
		for (L = 1; L <= NR; L++)
		{
			Z = 0.0;
			KP = 0;
			// FIND ELEMENT FOR ROW PIVOT
			for (K = L; K <= NR; K++)
			{
				AER = C.A[K][L];
				AER = DABS(AER);
				if (Z < AER)
				{
					Z = AER;
					KP = K;
				}
			}
			if (L < KP)
			{
				// INTERCHANGE ROWS
				for (J = L; J <= NPM; J++)
				{
					S = C.A[L][J];
					C.A[L][J] = C.A[KP][J];
					C.A[KP][J] = S;
				}
			}
			// TEST FOR A SINGULAR MATRIX
			AER = C.A[L][L];
			AER = DABS(AER);
			if (AER <= EX)
			{
				System.out.println(" MATRIX SINGULAR IN SUBROUTINE LSSSP");
				System.exit(0);
			}
			if (!(L < NR))
			{
				break;
			}
			LP1 = L + 1;
			for (K = LP1; K <= NR; K++)
			{
				AER = C.A[K][L];
				AER = DABS(AER);
				if (AER == 0)
				{
					RATIO = C.A[K][L] / C.A[L][L];
					for (J = LP1; J <= NPM; J++)
					{
						C.A[K][J] = C.A[K][J] - RATIO * C.A[L][J];
					}
				}
			}
		} //34
		C.print(" before back subst");
		// BACK SUBSTITUTION
		for (I = 1; I <= NR; I++)
		{
			S = 0.0;
			II = NPM - I;
			if (II < NR)
			{
				IIP1 = II + 1;
				for (K = IIP1; K <= NR; K++)
				{
					S = S + C.A[II][K] * X.x[K];
				}
			}
			RATIO = C.A[II][NPM];
			X.x[II] = (RATIO - S) / C.A[II][II];
		}
	}

	public static void main(String argv[])
	{
		double[][] C_data = { 
				{ 0, 0, 0, 0, 0 }, 
				{ 0, 1, 0, 1, 1}, 
				{ 0, 0, 3, 0, 2 }, 
				{ 0, 1, 0, 1, 3 }
		};

		System.out.println("Gaussian Elimination Test");
		System.out.println("FORTRAN Arrays starting at 1");
		Matrix C = new Matrix(C_data);
		C.print();
		VectorN X = new VectorN(4);
		LSSDP_solve(C, X);
		C.print();
		X.print();
	}
}
/*
// Set A
C.set(1,1,1.);
C.set(1,2,2.);
C.set(2,1,3.);
C.set(2,2,4.);
// set y
C.set(1,3,5.);
C.set(2,3,6.);
*/