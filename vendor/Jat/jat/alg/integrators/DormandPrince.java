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

/**
 * Implements a Dormand-Prince 8(7) integrator
 * adapted to JAT from the CODES software by Jim Baer
 * http://home.earthlink.net/~jimbaer1/
 * Only works with acceleration for derivs and 6-element state (r and v)
 * @author Tobias Berthold
 * @version 1.0
 */
public class DormandPrince
{
	// high (6e-9) and medium (1e-5) integration accuracy
	// public double integration_error = 6e-9;
	public double integration_error = 1.e-6;
	
	int counter_limit = 100;
	
	/*
	 * Since the dopri8 integration coefficients are reinitialized for every integration step, I'm going to try making
	 * them class variables; they will be initialized only once.
	 */
	double[][] dopri8_a = new double[14][14];
	double[] dopri8_c = new double[14];
	double[] dopri8_b = new double[14];
	double[] dopri8_bhat = new double[14];

	void initialize_dopri8()
	{
		/* Initialization of the integration coefficients for dopri8 */
		/* Initialize arrays */
		dopri8_a[2][1] = 1.0 / 18.0;
		dopri8_a[3][1] = 1.0 / 48.0;
		dopri8_a[3][2] = 1.0 / 16.0;
		dopri8_a[4][1] = 1.0 / 32.0;
		dopri8_a[4][2] = 0;
		dopri8_a[4][3] = 3.0 / 32.0;
		dopri8_a[5][1] = 5.0 / 16.0;
		dopri8_a[5][2] = 0;
		dopri8_a[5][3] = -75.0 / 64.0;
		dopri8_a[5][4] = 75.0 / 64.0;
		dopri8_a[6][1] = 3.0 / 80.0;
		dopri8_a[6][2] = 0;
		dopri8_a[6][3] = 0;
		dopri8_a[6][4] = 3.0 / 16.0;
		dopri8_a[6][5] = 3.0 / 20.0;
		dopri8_a[7][1] = 29443.841 / 614563.906;
		dopri8_a[7][2] = 0;
		dopri8_a[7][3] = 0;
		dopri8_a[7][4] = 77736.538 / 692538.347;
		dopri8_a[7][5] = -28693.883 / 1125000.000;
		dopri8_a[7][6] = 23124.283 / 1800000.000;
		dopri8_a[8][1] = 16016.141 / 946692.911;
		dopri8_a[8][2] = 0;
		dopri8_a[8][3] = 0;
		dopri8_a[8][4] = 61564.180 / 158732.637;
		dopri8_a[8][5] = 22789.713 / 633445.777;
		dopri8_a[8][6] = 545815.736 / 2771057.229;
		dopri8_a[8][7] = -180193.667 / 1043307.555;
		dopri8_a[9][1] = 39632.708 / 573591.083;
		dopri8_a[9][2] = 0;
		dopri8_a[9][3] = 0;
		dopri8_a[9][4] = -433636.366 / 683701.615;
		dopri8_a[9][5] = -421739.975 / 2616292.301;
		dopri8_a[9][6] = 100302.831 / 723423.059;
		dopri8_a[9][7] = 790204.164 / 839813.087;
		dopri8_a[9][8] = 800635.310 / 3783071.287;
		dopri8_a[10][1] = 246121.993 / 1340847.787;
		dopri8_a[10][2] = 0;
		dopri8_a[10][3] = 0;
		dopri8_a[10][4] = -37695042.795 / 15268766.246;
		dopri8_a[10][5] = -309121.744 / 1061227.803;
		dopri8_a[10][6] = -12992.083 / 490766.935;
		dopri8_a[10][7] = 6005943.493 / 2108947.869;
		dopri8_a[10][8] = 393006.217 / 1396673.457;
		dopri8_a[10][9] = 123872.331 / 1001029.789;
		dopri8_a[11][1] = -1028468.189 / 846180.014;
		dopri8_a[11][2] = 0;
		dopri8_a[11][3] = 0;
		dopri8_a[11][4] = 8478235.783 / 508512.852;
		dopri8_a[11][5] = 1311729.495 / 1432422.823;
		dopri8_a[11][6] = -10304129.995 / 1701304.382;
		dopri8_a[11][7] = -48777925.059 / 3047939.560;
		dopri8_a[11][8] = 15336726.248 / 1032824.649;
		dopri8_a[11][9] = -45442868.181 / 3398467.696;
		dopri8_a[11][10] = 3065993.473 / 597172.653;
		dopri8_a[12][1] = 185892.177 / 718116.043;
		dopri8_a[12][2] = 0;
		dopri8_a[12][3] = 0;
		dopri8_a[12][4] = -3185094.517 / 667107.341;
		dopri8_a[12][5] = -477755.414 / 1098053.517;
		dopri8_a[12][6] = -703635.378 / 230739.211;
		dopri8_a[12][7] = 5731566.787 / 1027545.527;
		dopri8_a[12][8] = 5232866.602 / 850066.563;
		dopri8_a[12][9] = -4093664.535 / 808688.257;
		dopri8_a[12][10] = 3962137.247 / 1805957.418;
		dopri8_a[12][11] = 65686.358 / 487910.083;
		dopri8_a[13][1] = 403863.854 / 491063.109;
		dopri8_a[13][2] = 0;
		dopri8_a[13][3] = 0;
		dopri8_a[13][4] = -5068492.393 / 434740.067;
		dopri8_a[13][5] = -411421.997 / 543043.805;
		dopri8_a[13][6] = 652783.627 / 914296.604;
		dopri8_a[13][7] = 11173962.825 / 925320.556;
		dopri8_a[13][8] = -13158990.841 / 6184727.034;
		dopri8_a[13][9] = 3936647.629 / 1978049.680;
		dopri8_a[13][10] = -160528.059 / 685178.525;
		dopri8_a[13][11] = 248638.103 / 1413531.060;
		dopri8_a[13][12] = 0;
		dopri8_c[1] = 0;
		dopri8_c[2] = 1.0 / 18.0;
		dopri8_c[3] = 1.0 / 12.0;
		dopri8_c[4] = 1.0 / 8.0;
		dopri8_c[5] = 5.0 / 16.0;
		dopri8_c[6] = 3.0 / 8.0;
		dopri8_c[7] = 59.0 / 400.0;
		dopri8_c[8] = 93.0 / 200.0;
		dopri8_c[9] = 5490023.248 / 9719169.821;
		dopri8_c[10] = 13.0 / 20.0;
		dopri8_c[11] = 1201146.811 / 1299019.798;
		dopri8_c[12] = 1;
		dopri8_c[13] = 1;
		dopri8_b[1] = 13451.932 / 455176.623;
		dopri8_b[2] = 0;
		dopri8_b[3] = 0;
		dopri8_b[4] = 0;
		dopri8_b[5] = 0;
		dopri8_b[6] = -808719.846 / 976000.145;
		dopri8_b[7] = 1757004.468 / 5645159.321;
		dopri8_b[8] = 656045.339 / 265891.186;
		dopri8_b[9] = -3867574.721 / 1518517.206;
		dopri8_b[10] = 465885.868 / 322736.535;
		dopri8_b[11] = 53011.238 / 667516.719;
		dopri8_b[12] = 2.0 / 45.0;
		dopri8_b[13] = 0;
		dopri8_bhat[1] = 14005.451 / 335480.064;
		dopri8_bhat[2] = 0;
		dopri8_bhat[3] = 0;
		dopri8_bhat[4] = 0;
		dopri8_bhat[5] = 0;
		dopri8_bhat[6] = -59238.493 / 1068277.825;
		dopri8_bhat[7] = 181606.767 / 758867.731;
		dopri8_bhat[8] = 561292.985 / 797845.732;
		dopri8_bhat[9] = -1041891.430 / 1371343.529;
		dopri8_bhat[10] = 760417.239 / 1151165.299;
		dopri8_bhat[11] = 118820.643 / 751138.087;
		dopri8_bhat[12] = -528747.749 / 2220607.170;
		dopri8_bhat[13] = 0.25;
	}

	/**
	 * Default constructor.
	 */
	public DormandPrince()
	{
		initialize_dopri8();
	}

	public void setIntegration_error(double integration_error)
	{
		this.integration_error = integration_error;
	}

	public double[] integrate(double time1, double state1[], double time2, Derivatives dv)
	{
		double[] state2 = new double[6];
		double r1[] = new double[4];
		double r1prime[] = new double[4];
		double r2[] = new double[4];
		double r2prime[] = new double[4];
		r1[1] = state1[0];
		r1[2] = state1[1];
		r1[3] = state1[2];
		r1prime[1] = state1[3];
		r1prime[2] = state1[4];
		r1prime[3] = state1[5];
		update(r1, r1prime, time1, r2, r2prime, time2,dv);
		state2[0] = r2[1];
		state2[1] = r2[2];
		state2[2] = r2[3];
		state2[3] = r2prime[1];
		state2[4] = r2prime[2];
		state2[5] = r2prime[3];
		return state2;
	}

	public void update(double r1[], double r1prime[], double time1, double r2[], double r2prime[], double time2, Derivatives dv)
	{
		// Procedure to calculate the position and velocity of a body at time2 given the position and velocity
		// of the body at time1. The Cowell method is used, with integration performed using Runge-Kutta
		// Prince-Dormand 8(7). Local error per unit step is controlled by variable step size.
		// Positions and velocities are given in A.U.s and A.U.s/day, while the times are Julian dates.
		// Reference for the Cowell method is Bate, Mueller, and White, "Fundamentals of Astrodynamics", pp 386-390.
		// Reference for RK dopri 8(7) is P.J. Prince and J.R. Dormand, "High-order embedded Runge-Kutta formulae",
		// J. Comp. Appl. Math., vol 7, 1981, p. 67-75.
		// The error control is adapted from the class notes of Math 551 taught at The Pennsylvania State University
		// during Fall 1989 by Prof. Douglas Arnold.
		double[] fatn = new double[4];
		double[] r = new double[4];
		double[] rprime = new double[4];
		double[] r_at_n = new double[4];
		double[] rprime_at_n = new double[4];
		double[] fatnp1 = new double[4];
		double toler = 0, h = 0, n1time = 0, ntime = 0, localerror = 0, ratio = 0;
		int k = 0, istep = 0, counter = 0;
		if (time1 == time2)
		{
			/* Return the state vector unchanged */
			for (k = 1; k <= 3; k++)
			{
				r2[k] = r1[k];
				r2prime[k] = r1prime[k];
			}
			return;
		}
		/* Set local error per unit step tolerance */
		toler = integration_error / Math.abs(time2 - time1);
		/* Determine the maximum number of attempted steps to be permitted */
		/* counter_limit = (int)(Math.abs(time2 - time1)/0.01); */
		/* Initialize variables for first time step */
		/* istep is the step being calculated */
		istep = 1;
		/* Set initial step size */
		h = time2 - time1;
		/* Limit initial h to +/- 20 days */
		//if (Math.abs(h) > 20.0)
			//h = 20.0 * h / Math.abs(h);
		/* n1time is the time at the (n+1)'st step, ntime is the time at the n'th step */
		n1time = time1;
		ntime = time1;
		/* Initialize position, velocity, and acceleration vectors */
		for (k = 1; k <= 3; k++)
		{
			r[k] = r1[k];
			r_at_n[k] = r1[k];
			rprime[k] = r1prime[k];
			rprime_at_n[k] = r1prime[k];
		}
		get_acceleration(ntime, r_at_n, rprime_at_n, dv, fatn);
		/* Begin integration loop */
		integrationloop: while ((((n1time < time2) && (time2 > time1)) || ((n1time > time2) && (time2 < time1))) && (counter < counter_limit))
		{
			/* So long as we've not reached time2 */
			/* Increment time to the (n+1)'st step */
			n1time = ntime + h;
			counter = counter + 1;
			/* if (Math.abs(h) < 0.01) System.out.println("step size = " + h); */
			/* System.out.println("step size = " + h); */
			/* Integrate to get r and rprime at the (n+1)'st time step */
			/* Use Runge-Kutta Dormand-Prince 8(7) */
			localerror = dopri8(ntime, h, fatn, rprime_at_n, r_at_n, rprime, r,dv);
			/*
			 * If a variant trajectory collides with a planet, the error may go to infinite or NaN. If so, end the
			 * integration.
			 */
			if (!(Math.abs(localerror) < 1e+10))
				counter = counter_limit;
			/* Compute ratio of projected new step size to present step size */
			if (Math.abs(localerror / h) > 0.0)
				/* Avoid a division by zero for small step sizes */
				ratio = .8 * Math.pow(Math.abs(toler / (localerror / h)), (1.0 / 7.0));
			else
				/* Error and step size are very small; allow step size to double */
				ratio = 2;
			/* Test to see if step will be accepted; if already at minimum step size, ignore */
			if ((Math.abs(localerror / h) > toler) && (Math.abs(h) > 0.001))
			{
				/* Step has failed */
				/* Determine a new step size */
				if ((ratio < .2) && (istep > 1))
					ratio = .2;
				/*
				 * Note that new step size is limited to one-fifth of the previous one, except on the first step, which
				 * allows faster determination of an initial step size
				 */
				n1time = ntime;
				h = ratio * h;
				/* Enforce minimum step size to prevent ultrasmall steps near collision */
				if (Math.abs(h) < 0.001)
					h = 0.001 * (h / Math.abs(h));
				continue integrationloop;
			}
			/* Step has been accepted, or we're at the minimum step size. In either case, reset the counter. */
			counter = 0;
			/* Determine a new step size */
			if (ratio > 2)
				ratio = 2;
			/* Note that new step size is limited to double the previous one */
			h = ratio * h;
			/* Enforce minimum step size to prevent ultrasmall steps near collision */
			if (Math.abs(h) < 0.001)
				h = 0.001 * (h / Math.abs(h));
			/* Calculate second derivative of r at time n+1 */
			get_acceleration(n1time, r, rprime, dv, fatnp1);
			/* Increment step counter */
			istep = istep + 1;
			/* Check to prevent overrunning endpoint of integration */
			if (((((n1time + h) > time2) && (time2 > time1)) || (((n1time + h) < time2) && (time2 < time1))) && (n1time != time2))
				/* Reduce step size to avoid overrunning endpoint of integration */
				h = time2 - n1time;
			/*
			 * If a variant trajectory collides with a planet, the acceleration may go to infinite or NaN. If so, end
			 * the integration.
			 */
			if (!(Math.abs(fatnp1[1]) < 1e+10) || !(Math.abs(fatnp1[2]) < 1e+10) || !(Math.abs(fatnp1[3]) < 1e+10))
				counter = counter_limit;
			/* Update integration variables */
			ntime = n1time;
			for (k = 1; k <= 3; k++)
			{
				r_at_n[k] = r[k];
				rprime_at_n[k] = rprime[k];
				fatn[k] = fatnp1[k];
			}
		}
		if (!(counter < counter_limit))
		{
			System.out.println("WARNING: update exceeds max steps");
		}
		/*
		 * Integration complete; return new state vector
		 */
		for (k = 1; k <= 3; k++)
		{
			r2[k] = r[k];
			r2prime[k] = rprime[k];
		}
	}

	double dopri8(double ntime, double h, double fatn[], double rprimen[], double rn[], double rprime[], double r[], Derivatives dv)
	{
		/*
		 * Procedure to integrate the equations of motion (Cowell's method) from ntime to ntime+h using the
		 * Runge-Kutta-type Prince-Dormand 8(7) method (see procedure update). Also finds local error in the integration
		 * step, which is returned. Reference: P.J. Prince and J.R. Dormand, "High-order embedded Runge-Kutta formulae",
		 * J. Comp. Appl. Math., vol 7, 1981, p. 67-75. Also, E. Hairer, S.P. Norsett, G. Wanner, "Solving Ordinary
		 * Differential Equations I (Non-stiff problems), Springer-Verlag, 1987, pp. 132-133, 167-171, 193-195, and 260.
		 */
		double[][] dopri8_kprime = new double[14][4];
		double[][] dopri8_k = new double[14][4];
		double[] low_r = new double[4];
		double[] rvariant = new double[4];
		double[] rprimevariant = new double[4];
		double[] kprime_temp = new double[4];
		double localerror = 0;
		int i = 0, j = 0, n = 0;
		/* Begin the integration from ntime to ntime+h */
		/* Stage 1 */
		for (n = 1; n <= 3; n++)
		{
			/* Initialize the Runge-Kutta variables, and the pos/vel vectors. */
			dopri8_kprime[1][n] = fatn[n];
			dopri8_k[1][n] = rprimen[n];
			/* r is the 8'th order position solution; low_r is the 7'th order solution */
			r[n] = rn[n] + h * dopri8_bhat[1] * dopri8_k[1][n];
			low_r[n] = rn[n] + h * dopri8_b[1] * dopri8_k[1][n];
			rprime[n] = rprimen[n] + h * dopri8_bhat[1] * dopri8_kprime[1][n];
		}
		/* Stages 2-13 */
		for (i = 2; i <= 13; i++)
		{
			for (n = 1; n <= 3; n++)
			{
				/*
				 * rvariant and rprimevariant are the intermediate pos/vel used to evaluate the acceleration function at
				 * each stage
				 */
				rvariant[n] = rn[n];
				rprimevariant[n] = rprimen[n];
				dopri8_k[i][n] = rprimen[n];
				for (j = 1; j <= (i - 1); j++)
				{
					/*
					 * Sum the terms in k, rvariant, and rprimevariant needed to evaluate the acceleration function and
					 * determine kprime at this stage
					 */
					rvariant[n] = rvariant[n] + h * (dopri8_a[i][j] * dopri8_k[j][n]);
					rprimevariant[n] = rprimevariant[n] + h * (dopri8_a[i][j] * dopri8_kprime[j][n]);
					dopri8_k[i][n] = dopri8_k[i][n] + h * (dopri8_a[i][j] * dopri8_kprime[j][n]);
				}
			}
			get_acceleration(ntime + dopri8_c[i] * h, rvariant, rprimevariant, dv, kprime_temp);
			for (n = 1; n <= 3; n++)
			{
				dopri8_kprime[i][n] = kprime_temp[n];
				r[n] = r[n] + h * dopri8_bhat[i] * dopri8_k[i][n];
				low_r[n] = low_r[n] + h * dopri8_b[i] * dopri8_k[i][n];
				rprime[n] = rprime[n] + h * dopri8_bhat[i] * dopri8_kprime[i][n];
			}
		}
		/* Compute the local error in position */
		localerror = 0;
		for (n = 1; n <= 3; n++)
			localerror = localerror + Math.pow((r[n] - low_r[n]), 2);
		localerror = Math.sqrt(localerror);
		return localerror;
	}

	void get_acceleration(double jultime, double r[], double rprime[], Derivatives dv, double acceleration[])
	{
		// Procedure to calculate the accelerations on a body for use in the differential equations of motion
		// in Cowell's method.
		//int k = 0;
		/* Calculate contribution to acceleration of object */
		//for (k = 1; k <= 3; k++)
			// acceleration[k] = 1.;
			//acceleration[k] = -r[k];
		double[] state=new double[6];
		state[0] = r[1];
		state[1] = r[2];
		state[2] = r[3];
		state[3] = rprime[1];
		state[4] = rprime[2];
		state[5] = rprime[3];
		double[] deriv=dv.derivs(jultime,state);
		acceleration[1] = deriv[3];
		acceleration[2] = deriv[4];
		acceleration[3] = deriv[5];
	
	}

	public static void main(String[] args)
	{

	}

}
