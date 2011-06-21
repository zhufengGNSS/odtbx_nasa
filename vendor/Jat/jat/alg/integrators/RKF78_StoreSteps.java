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
package jat.alg.integrators;

import java.util.ArrayList;

/**
 * This class contains a set of functions (rkck(), rkqs(), and odeint()), which are each detailed in their individual
 * function headers. These functions are intended to work together to create an adaptive stepsize numerical integrator
 * routine.
 * @author Dave Gaylor, Brent Barbee, Tobias Berthold
 * @version 1.0
 */
public class RKF78_StoreSteps
{
	// Variables used as inputs for user options
	public double eps = 1.e-10; // desired accuracy
	public double dxsav = 0.001; // Minimum difference between indep. var. values to store steps
	public boolean store_steps = true;
	public int MAXSTEP = 10000;
	public double stepSize_ = 0.01;
	public double minStepSize_ = 1.2E-10;
	// Variables used as output
	public double hdid; // the step size that was actually applied
	public double hnext; // estimated next step size
	public int nok; // Total steps
	public int nbad;// Total step adjustments
	public double[][] steplist; // Array storing intermediate steps
	// Variables internally used by several methods
	public double x; // independent variable of integration, used by rkqs and odeint
	private int ntry;
	// Constants
	static double TINY = 1.0e-30;
	// Brent:
	// static double SAFETY = 0.9;
	// static double PGROW = -0.2;
	// static double PSHRNK = -0.25;
	// static double ERRCON = 1.89e-4; // This value is equal to (5/SAFETY raised to the power (1/PGROW)
	// Dave:
	private static final double SAFETY = 0.9;
	private static final double PGROW = -1.0 / 8.0;
	private static final double PSHRNK = -1.0 / 7.0;
	private static final double ERRCON = 2.56578451395034701E-8;

	/**
	 * Propagate a state vector
	 * @param start initial value of independent variable
	 * @param x0 initial value of state vector
	 * @param end final value of independent variable
	 * @param dv class implementing the Derivatives interface
	 * @return state vector at the end of the interval
	 */
	public double[] integrate(double start, double[] x0, double end, Derivatives dv)
	{
		double[] y = x0.clone();
		double h = stepSize_;
		double hmin = minStepSize_;
		odeint(y, start, end, h, hmin, dv);
		return y;
	}

	/**
	 * Driver program for Runge Kutta numerical integrator that employs adaptive
	 * stepsize control and has 8th order error estimation
	 * @param ystart Starting dependent variable vector
	 * @param x1 Starting value of the independent variable
	 * @param x2 Ending value of the independent variable
	 * @param h1 Guessed first stepsize
	 * @param hmin Minimum allowed stepsize (may be zero)
	 * @param dv user-supplied routine that computes the derivatives
	 */
	public void odeint(double ystart[], double x1, double x2, double h1, double hmin, Derivatives dv)
	{
		int nvar = ystart.length;
		int nstp, i;
		double xsav = 0., h;
		double[] yscal = new double[nvar];
		double[] y = new double[nvar];
		double[] dydx = new double[nvar];
		ArrayList steps = new ArrayList();
		//
		x = x1;
		// If backward integration required, make stepsize negative
		h = Math.abs(h1);
		if (x2 < x1)
			h = -h;
		// Initialize the counters for the number of good and bad (but corrected) integration steps.
		nok = nbad = 0;
		// Copy initial state to temp array
		y = ystart.clone();
		// Store initial values as first entry in intermediate storage array
		if (store_steps)
			store_step(steps, nvar, x, y, ntry);
		// Main integration loop
		for (nstp = 1; nstp <= MAXSTEP; nstp++)
		{
			dydx = dv.derivs(x, y);
			// Scaling used to monitor accuracy
			for (i = 0; i < nvar; i++)
				yscal[i] = Math.abs(y[i]) + Math.abs(dydx[i] * h) + TINY;
			// If stepsize can overshoot, decrease
			if ((x + h - x2) * (x + h - x1) > 0.0)
				h = x2 - x;
			// Take Runge Kutta step with monitoring of local truncation error
			rkqs(y, dydx, h, yscal, dv);
			// Store intermediate results
			if (store_steps && Math.abs(x - xsav) > Math.abs(dxsav))
			{
				store_step(steps, nvar, x, y, ntry);
				xsav = x;
			}
			// If rkqs did not decrease the step size, it was a good step
			if (hdid == h)
				nok++;
			else
				nbad++;
			// If the independent variable has been advanced to the final value, then the integration is complete.
			if ((x - x2) * (x2 - x1) >= 0.0)
			{
				for (i = 0; i < nvar; i++)
					ystart[i] = y[i];
				// convert list to array
				if (store_steps)
				{
					// store_step(steps, nvar, x, y);
					steplist = steps_to_array(steps, nvar);
				}
				return;
			}
			if (Math.abs(hnext) <= hmin)
				error("Step size too small in odeint");
			h = hnext;
		}
		error("Too many steps in routine odeint");
	}

	/**
	 * Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust stepsize
	 * @param y Dependent variable vector
	 * @param dydx Derivative of the dependent variable vector
	 * @param htry Stepsize to attempt
	 * @param yscal
	 * @param dv user-supplied routine that computes the derivatives
	 */
	public void rkqs(double y[], double dydx[], double htry, double[] yscal, Derivatives dv)
	{
		int n = y.length;
		int i;
		double errmax, h, htemp, xnew, yerr[], ytemp[];
		yerr = new double[n];
		ytemp = new double[n];
		// Set the stepsize to the size that is to be tried first
		h = htry;
		ntry = 0;
		for (;;)
		{
			// Take a Runge-Kutta Cash-Karp integration step
			rkck(y, dydx, x, h, ytemp, yerr, dv);
			// Initialize the maximum error value in the state vector to zero
			errmax = 0.0;
			// Evaluate accuracy
			for (i = 0; i < n; i++)
			{
				// The largest absolute valued element in the error vector is considered the maximum error value
				errmax = Math.max(errmax, Math.abs(yerr[i]) / yscal[i]);
			}
			// Scale the maximum error relative to the required tolerance
			errmax /= eps;
			// If the error is within tolerance then break out of this loop
			if (errmax <= 1.0)
			{
				// Step succeeded. Break and go on to compute the size of the next step
				break;
			}
			// Need more steps
			ntry++;
			// Temporary stepsize value.
			htemp = SAFETY * h * Math.pow(errmax, PSHRNK);
			// Truncation error too large, reduce stepsize (by no more than a factor of 10)
			h = (h >= 0.0 ? Math.max(htemp, 0.1 * h) : Math.min(htemp, 0.1 * h));
			// Compute the next value of the independent variable with the given stepsize
			xnew = x + h;
			// If the independent variable has not been changed, we have a stepsize underflow
			if (xnew == x)
			{
				// Error-out with a message
				error("Stepsize underflow in rkqs().");
				// TODO: exception handling
			}
		}
		// Alter the stepsize, but allow no more than a factor of 5 increase
		if (errmax > ERRCON)
		{
			hnext = SAFETY * h * Math.pow(errmax, PGROW);
		} else
		{
			hnext = 5.0 * h;
		}
		// Update the independent variable
		hdid = h;
		x += hdid;
		// Update the state vector
		for (i = 0; i < n; i++)
		{
			y[i] = ytemp[i];
		}
	}

	private static final double[] a =
	{ 0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 0.5, 5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0 };
	private static final double[][] b = new double[13][12];
	static
	{
		for (int i = 0; i < 13; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				b[i][j] = 0.0;
			}
		}
		b[1][0] = 2.0 / 27.0;
		b[2][0] = 1.0 / 36.0;
		b[2][1] = 1.0 / 12.0;
		b[3][0] = 1.0 / 24.0;
		b[3][2] = 1.0 / 8.0;
		b[4][0] = 5.0 / 12.0;
		b[4][2] = -25.0 / 16.0;
		b[4][3] = 25.0 / 16.0;
		b[5][0] = 1.0 / 20.0;
		b[5][3] = 0.25;
		b[5][4] = 0.2;
		b[6][0] = -25.0 / 108.0;
		b[6][3] = 125.0 / 108.0;
		b[6][4] = -65.0 / 27.0;
		b[6][5] = 125.0 / 54.0;
		b[7][0] = 31.0 / 300.0;
		b[7][4] = 61.0 / 225.0;
		b[7][5] = -2.0 / 9.0;
		b[7][6] = 13.0 / 900.0;
		b[8][0] = 2.0;
		b[8][3] = -53.0 / 6.0;
		b[8][4] = 704.0 / 45.0;
		b[8][5] = -107.0 / 9.0;
		b[8][6] = 67.0 / 90.0;
		b[8][7] = 3.0;
		b[9][0] = -91.0 / 108.0;
		b[9][3] = 23.0 / 108.0;
		b[9][4] = -976.0 / 135.0;
		b[9][5] = 311.0 / 54.0;
		b[9][6] = -19.0 / 60.0;
		b[9][7] = 17.0 / 6.0;
		b[9][8] = -1.0 / 12.0;
		b[10][0] = 2383.0 / 4100.0;
		b[10][3] = -341.0 / 164.0;
		b[10][4] = 4496.0 / 1025.0;
		b[10][5] = -301.0 / 82.0;
		b[10][6] = 2133.0 / 4100.0;
		b[10][7] = 45.0 / 82.0;
		b[10][8] = 45.0 / 164.0;
		b[10][9] = 18.0 / 41.0;
		b[11][0] = 3.0 / 205.0;
		b[11][5] = -6.0 / 41.0;
		b[11][6] = -3.0 / 205.0;
		b[11][7] = -3.0 / 41.0;
		b[11][8] = 3.0 / 41.0;
		b[11][9] = 6.0 / 41.0;
		b[12][0] = -1777.0 / 4100.0;
		b[12][3] = -341.0 / 164.0;
		b[12][4] = 4496.0 / 1025.0;
		b[12][5] = -289.0 / 82.0;
		b[12][6] = 2193.0 / 4100.0;
		b[12][7] = 51.0 / 82.0;
		b[12][8] = 33.0 / 164.0;
		b[12][9] = 12.0 / 41.0;
		b[12][11] = 1.0;
	}
	private static final double[] c =
	{ 41.0 / 840.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 41.0 / 840.0, 0.0, 0.0 };
	private static final double[] chat =
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 0.0, 41.0 / 840.0, 41.0 / 840.0 };

	/**
	 * This function will use the Runge-Kutta method to advance the solution over an interval h and return the
	 * incremented variables
	 * @param y Dependent variable vector
	 * @param dydx Derivative of the dependent variable vector
	 * @param x independent variable
	 * @param h Stepsize
	 * @param yout Incremented dependent variable vector
	 * @param yerr Estimate of the local truncation error in yout using the embedded 7th order method
	 * @param dv user-supplied routine that computes the derivatives
	 */
	protected void rkck(double[] y, double[] dydx, double x, double h, double[] yout, double[] yerr, Derivatives dv)
	{
		int n = y.length;
		double f[][] = new double[13][n];
		double yt[][] = new double[13][n];
		// double y7th[] = new double[n];
		double ytmp[] = new double[n];
		double sum[] = new double[n];
		// System.out.println("step size = "+h);
		double xeval[] = new double[13];
		for (int i = 0; i < 13; i++) // find times for function evals
		{
			xeval[i] = x + a[i] * h;
		}
		// build f matrix
		// f[0] = derivs.derivs(x, y);
		f[0] = dydx;
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * b[1][0] * f[0][i];
		}
		f[1] = dv.derivs(xeval[1], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[2][0] * f[0][i] + b[2][1] * f[1][i]);
		}
		f[2] = dv.derivs(xeval[2], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[3][0] * f[0][i] + b[3][2] * f[2][i]);
		}
		f[3] = dv.derivs(xeval[3], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[4][0] * f[0][i] + b[4][2] * f[2][i] + b[4][3] * f[3][i]);
		}
		f[4] = dv.derivs(xeval[4], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[5][0] * f[0][i] + b[5][3] * f[3][i] + b[5][4] * f[4][i]);
		}
		f[5] = dv.derivs(xeval[5], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[6][0] * f[0][i] + b[6][3] * f[3][i] + b[6][4] * f[4][i] + b[6][5] * f[5][i]);
		}
		f[6] = dv.derivs(xeval[6], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[7][0] * f[0][i] + b[7][4] * f[4][i] + b[7][5] * f[5][i] + b[7][6] * f[6][i]);
		}
		f[7] = dv.derivs(xeval[7], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[8][0] * f[0][i] + b[8][3] * f[3][i] + b[8][4] * f[4][i] + b[8][5] * f[5][i] + b[8][6] * f[6][i] + b[8][7] * f[7][i]);
		}
		f[8] = dv.derivs(xeval[8], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[9][0] * f[0][i] + b[9][3] * f[3][i] + b[9][4] * f[4][i] + b[9][5] * f[5][i] + b[9][6] * f[6][i] + b[9][7] * f[7][i] + b[9][8] * f[8][i]);
		}
		f[9] = dv.derivs(xeval[9], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h
					* (b[10][0] * f[0][i] + b[10][3] * f[3][i] + b[10][4] * f[4][i] + b[10][5] * f[5][i] + b[10][6] * f[6][i] + b[10][7] * f[7][i] + b[10][8] * f[8][i] + b[10][9] * f[9][i]);
		}
		f[10] = dv.derivs(xeval[10], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i] + h * (b[11][0] * f[0][i] + b[11][5] * f[5][i] + b[11][6] * f[6][i] + b[11][7] * f[7][i] + b[11][8] * f[8][i] + b[11][9] * f[9][i]);
		}
		f[11] = dv.derivs(xeval[11], ytmp);
		for (int i = 0; i < n; i++)
		{
			ytmp[i] = y[i]
					+ h
					* (b[12][0] * f[0][i] + b[12][3] * f[3][i] + b[12][4] * f[4][i] + b[12][5] * f[5][i] + b[12][6] * f[6][i] + b[12][7] * f[7][i] + b[12][8] * f[8][i] + b[12][9] * f[9][i] + f[11][i]);
		}
		f[12] = dv.derivs(xeval[12], ytmp);
		// construct solutions
		// yout is the 8th order solution
		for (int i = 0; i < n; i++)
		{
			// y7th[i] = y[i] + h*(c[0]*f[0][i] +c[5]*f[5][i] + c[6]*f[6][i] + c[7]*f[7][i] + c[8]*f[8][i] +
			// c[9]*f[9][i] + c[10]*f[10][i]);
			yout[i] = y[i] + h * (chat[5] * f[5][i] + chat[6] * f[6][i] + chat[7] * f[7][i] + chat[8] * f[8][i] + chat[9] * f[9][i] + chat[11] * f[11][i] + chat[12] * f[12][i]);
			yerr[i] = h * c[0] * (f[11][i] + f[12][i] - f[0][i] - f[10][i]);
		}
	}

	/**
	 * Store one integration step in an ArrayList
	 * @param steps The ArrayList to store the step into
	 * @param nvar dimension of state vector
	 * @param x independent variable 
	 * @param y state vector
	 * @param ntry number of steps to achieve requested accuracy
	 */
	private void store_step(ArrayList steps, int nvar, double x, double[] y, int ntry)
	{
		double[] step = new double[nvar + 2];
		step[0] = x;
		for (int i = 0; i < nvar; i++)
			step[i + 1] = y[i];
		step[nvar + 1] = ntry;
		steps.add(step);
	}

	/**
	 * Convert an ArrayList to a two-dimensional double array
	 * @param steps the ArrayList
	 * @param nvar dimension of state vector
	 * @return
	 */
	private double[][] steps_to_array(ArrayList steps, int nvar)
	{
		int number_of_steps = steps.size();
		double[][] steps_array = new double[number_of_steps][nvar + 2];
		for (int i = 0; i < number_of_steps; i++)
		{
			double[] step = (double[]) steps.get(i);
			for (int j = 0; j < step.length; j++)
			{
				steps_array[i][j] = step[j];
			}
		}
		return steps_array;
	}

	private void error(String msg)
	{
		System.err.println("RKF78_PrintSteps: " + msg);
	}
}
