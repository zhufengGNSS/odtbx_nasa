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
package jat.cm;
import jat.alg.*;
import jat.matvec.data.*;
import jat.alg.integrators.*;

/**
 * <P>
 * The Lambert class provides the means to solve Lambert's problem.
 * Ref: Vallado
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class Lambert implements ScalarFunction
{

	private double s = 0.0;
	private double c = 0.0;
	private double dt = 0.0;
	private double mu = 0.0;
	private boolean aflag = false;
	private boolean bflag = false;
	public boolean debug_print=false;

	/** Contains the computed initial delta-v.
	 */
	public VectorN deltav0;
	/** Contains the computed final delta-v.
	 */
	public VectorN deltavf;

	/** time of flight */
	public double tof;

	
	public void reset()
	{
		s = 0.0;
		c = 0.0;
		aflag = false;
		bflag = false;
		debug_print=false;
	}
	
	
    /** Constructor with mu
     * @param mu mu of the central body
     */
	public Lambert(double mu)
	{
		this.mu=mu;
	}

	private double getalpha(double a)
	{
		double alpha = 2.0 * Math.asin(Math.sqrt(s / (2.0 * a)));
		if (this.aflag)
		{
			alpha = 2.0 * Constants.pi - alpha;
		}
		return alpha;
	}

	private double getbeta(double a)
	{
		double beta = 2.0 * Math.asin(Math.sqrt((s - c) / (2.0 * a)));
		if (this.bflag)
		{
			beta = -1.0 * beta;
		}
		return beta;
	}

	private double getdt(double a, double alpha, double beta)
	{
		double sa = Math.sin(alpha);
		double sb = Math.sin(beta);
		double dt =
			Math.pow(a, 1.5) * (alpha - sa - beta + sb) / Math.sqrt(mu);
		return dt;
	}

	/** Evaluate the delta-t function f
	 * @param a semi-major axis.
	 * @return
	 */
	public double evaluate(double a)
	{
		double alpha = getalpha(a);
		double beta = getbeta(a);
		double out = dt - getdt(a, alpha, beta);
		return out;
	}

	/** Computes the delta-v's required to go from r0,v0 to rf,vf.
         * @return Total delta-v (magnitude) required.
         * @param dt Time of flight
         * @param r0 Initial position vector.
         * @param v0 Initial velocity vector.
         * @param rf Desired final position vector.
         * @param vf Desired final velocity vector.
         */
	public double compute(
		VectorN r0,
		VectorN v0,
		VectorN rf,
		VectorN vf,
		double dt)
	{
		reset();
		double tp = 0.0;

		this.dt = dt;
		double magr0 = r0.mag();
		double magrf = rf.mag();

		VectorN dr = r0.minus(rf);
		this.c = dr.mag();
		this.s = (magr0 + magrf + c) / 2.0;
		double amin = s / 2.0;
		if(debug_print)
			System.out.println("amin = " + amin);

		double dtheta = Math.acos(r0.dotProduct(rf) / (magr0 * magrf));

		//dtheta = 2.0 * Constants.pi - dtheta;

		if(debug_print)
			System.out.println("dtheta = " + dtheta);

		if (dtheta < Constants.pi)
		{
			tp =
				Math.sqrt(2.0 / (mu))
					* (Math.pow(s, 1.5) - Math.pow(s - c, 1.5))
					/ 3.0;
		}
		if (dtheta > Constants.pi)
		{
			tp =
				Math.sqrt(2.0 / (mu))
					* (Math.pow(s, 1.5) + Math.pow(s - c, 1.5))
					/ 3.0;
			this.bflag = true;
		}
		
		if(debug_print)
			System.out.println("tp = " + tp);

		double betam = getbeta(amin);
		double tm = getdt(amin, Constants.pi, betam);

		if(debug_print)
			System.out.println("tm = " + tm);

		if (dtheta == Constants.pi)
		{
			System.out.println(" dtheta = 180.0. Do a Hohmann");
			System.exit(0);
		}

		double ahigh = 1000.0 * amin;
		double npts = 3000.0;
		//this.dt = (2.70-0.89)*86400;
		if(debug_print)
			System.out.println("dt = " + dt);

		if(debug_print)
			System.out.println("************************************************");

		if (this.dt < tp)
		{
			System.out.println(" No elliptical path possible ");
			System.exit(0);
		}

		if (this.dt > tm)
		{
			this.aflag = true;
		}

		double fm = evaluate(amin);
		double ftemp = evaluate(ahigh);

		if ((fm * ftemp) >= 0.0)
		{
			System.out.println(" initial guesses do not bound ");
			System.exit(0);
		}

		ZeroFinder regfalsi = new ZeroFinder(this, 10000, 1.0E-6, 1.0E-15);

		double sma = regfalsi.regulaFalsi(amin, ahigh);

		double alpha = getalpha(sma);
		double beta = getbeta(sma);

		double de = alpha - beta;

		double f = 1.0 - (sma / magr0) * (1.0 - Math.cos(de));
		double g = dt - Math.sqrt(sma * sma * sma / mu) * (de - Math.sin(de));

		VectorN newv0 = new VectorN(3);
		VectorN newvf = new VectorN(3);

		newv0.x[0] = (rf.x[0] - f * r0.x[0]) / g;
		newv0.x[1] = (rf.x[1] - f * r0.x[1]) / g;
		newv0.x[2] = (rf.x[2] - f * r0.x[2]) / g;

		this.deltav0 = newv0.minus(v0);
		if(debug_print)
			this.deltav0.print("deltav-0");

		double dv0 = deltav0.mag();

		double fdot =
			-1.0 * (Math.sqrt(mu * sma) / (magr0 * magrf)) * Math.sin(de);
		double gdot = 1.0 - (sma / magrf) * (1.0 - Math.cos(de));

		newvf.x[0] = fdot * r0.x[0] + gdot * newv0.x[0];
		newvf.x[1] = fdot * r0.x[1] + gdot * newv0.x[1];
		newvf.x[2] = fdot * r0.x[2] + gdot * newv0.x[2];

		this.deltavf = vf.minus(newvf);
		double dvf = deltavf.mag();
		if(debug_print)
			this.deltavf.print("deltav-f");

		double totaldv = dv0 + dvf;

		this.tof = dt;

		if(debug_print)
			System.out.println(
			"dt = "
				+ dt
				+ " dv0 = "
				+ dv0
				+ " dvf = "
				+ dvf
				+ " total dv = "
				+ totaldv
				+ " sma = "
				+ sma);
		return totaldv;
	}

	/** Test case.
	 * @param args arguments (none).
	 */
	public static void main(String args[])
	{
		LinePrinter lp = new LinePrinter();

		TwoBody elem0 = new TwoBody(40000.0, 0.2, 0.0, 0.0, 45.0, 0.0);
		TwoBody elemf = new TwoBody(80000.0, 0.2, 0.0, 0.0, 270.0, 286.0);

		elem0.propagate(0.0, (0.89 * 86400.0), lp, false);
		elemf.propagate(0.0, (2.70 * 86400.0), lp, false);

		elem0.print("SC1");
		elemf.print("SC2");

		VectorN r0 = elem0.getR();
		VectorN v0 = elem0.getV();
		VectorN rf = elemf.getR();
		VectorN vf = elemf.getV();

		Lambert lambert = new Lambert(Constants.mu);
		double totaldv = lambert.compute(r0, v0, rf, vf, (2.70 - 0.89) * 86400);

	}
}
