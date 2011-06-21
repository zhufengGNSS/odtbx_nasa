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
 
package jat.matvec.data.matrixDecompositions;
import jat.matvec.data.*;

/**
 * <P>
 * The Householder Class provides the solution to the equation y = Hx, where we are solving for x,
 * a double array by using a Householder decomposition.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class Householder
// Householder
{
	public double xhat[];
	public Matrix R;
	public Matrix P;
	int n;
	int m;

	public Householder(Matrix H, double [] y)
	{
		this.n = H.n;
		this.m = y.length;
		H.checkRowDimension(y.length);
		double [] out = new double[n];
		out = this. compute(H, y);
	}

	public Householder(int n_in, int m_in)
	{
		this.n = n_in;
		this.m = m_in;
	}

	public double [] compute( Matrix H, double [] y)
	{
//		long ts = System.currentTimeMillis();

		Matrix A = new Matrix((n+m),(n+1));

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				A.A[i][j] = 0.0;
			}
			A.A[i][n] = 0.0;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				A.A[i+n][j] = H.A[i][j];
			}
			A.A[i+n][n] = y[i];
		}

		this.R = new Matrix(n,n);
		this.xhat = new double[n];
		this.P = new Matrix(n,n);
		Matrix Rtrans = new Matrix(n,n);
		Matrix RTR = new Matrix(n,n);

		double sigma = 0.0;
		double [] u = new double[m+n];

		for (int k = 0; k < n; k++)
		{
			double sum1 = 0.0;
			for (int i = k; i < (m+n); i++)
			{
				sum1 = sum1 + A.A[i][k]*A.A[i][k];
			}
			if (A.A[k][k] < 0.0)
			{
				sigma = -1.0 * Math.sqrt(sum1);
			}
			else
			{
				sigma = Math.sqrt(sum1);
			}

			u[k] = A.A[k][k] + sigma;
			A.A[k][k] = -1.0 * sigma;

			for (int i = (k+1); i < (m+n); i++)
			{
				u[i] = A.A[i][k];
			}

			double beta = 1.0/(sigma*u[k]);

			for (int j = (k+1); j < (n+1); j++)
			{
				double sum2 = 0.0;
				for (int i = k; i < (m+n); i++)
				{
					sum2 = sum2 + u[i]*A.A[i][j];
				}
				double gamma = beta * sum2;
				for (int i = k; i < (m+n); i++)
				{
					A.A[i][j] = A.A[i][j] - gamma * u[i];
				}
			} // Next j

			for (int i = (k+1); i < (m+n); i++)
			{
				A.A[i][k] = 0.0;
			}
		} // Next k

		double sum = 0.0;
		for (int i = 0; i < m; i++)
		{
			sum = sum + A.A[n+i][n]*A.A[n+i][n];
		}


        double [] b = new double[n];
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				this.R.A[i][j] = A.A[i][j];
			}
			b[i] = A.A[i][n];
			this.xhat[i] = 0.0;
		}


		for (int i = (n-1); i > -1; i--)
		{
			double rijxj = 0.0;
			for (int j = (i+1); j < n; j++)
			{
				rijxj = rijxj + R.A[i][j]*this.xhat[j];
			}

			this.xhat[i] = (b[i] - rijxj)/R.A[i][i];
		}

		Rtrans = this.R.transpose();
		RTR = Rtrans.times(this.R);

		return this.xhat;

	}

}