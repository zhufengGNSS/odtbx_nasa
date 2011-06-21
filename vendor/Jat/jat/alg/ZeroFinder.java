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

package jat.alg;
import jat.math.*;

/**
 * <P>
 * The ZeroFinder Class provides a way to solve a scalar f(x) = 0.
 * These functions have been translated from Numerical Recipes.
 * Currently there are: Regula Falsi, Secant and Ridder's methods.
 * The function f is passed via the ScalarFunction interface.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class ZeroFinder {
	
	/** maximum number of iterations */
	private int maxit = 500;
	
	/** accuracy of the function, how close to f = 0 */
	private double accuracy = 1.0E-12;
	
	/** accuracy of the solution */
	private double dxmin = 1.0E-12;
	
	/** the function to be solved */
	private ScalarFunction func;
	
	/** bad number */
	private double UNUSED = -1.11E+30;
		
	/** error counter */
	private int err = 0;
	
	public double bracket1;
	
	public double bracket2;
	
	public double[] xb1;
	public double[] xb2;

	/** Construct a zero finding problem
	 * @param f ScalarFunction to be solved.
	 */
	public ZeroFinder(ScalarFunction f) {
		this.func = f;
	}

	/** Construct a zero finding problem.
	 * @param f ScalarFunction to be solved.
	 * @param max Maximum number of iterations.
	 * @param acc How close does f(x) have to be to zero?
	 * @param dxm How close to the zero does x have to be?
	 */
	public ZeroFinder(ScalarFunction f, int max, double acc, double dxm) {
		setMaxIterations(max);
		setAccuracy(acc);
		set_dxmin(dxm);
		this.func = f;
	}

	/** Set the maximum number of iterations
	 * @param max Maximum number of iterations.
	 */
	public void setMaxIterations(int max) {
		this.maxit = max;
	}

	/** How close does f(x) have to be to zero?
	 * @param acc How close does f(x) have to be to zero?
	 */
	public void setAccuracy(double acc) {
		this.accuracy = acc;
	}

	/** How close to the zero does x have to be?
	 * @param x How close to the zero does x have to be?
	 */
	public void set_dxmin(double x) {
		this.dxmin = x;
	}
	
	/** Find a pair of brackets inside or outside the interval (x1, x2). From Numerical Recipes, chap 9.1
	 * To access the brackets, get this.bracket1 and this.bracket2
	 * @param x1 lower value
	 * @param x2 upper value
	 * @return boolean true=success, false=failure
	 */			
	public boolean zbrac(double x1, double x2){
		double factor = 1.6;
		int ntry = 50;
		
		boolean out = false;
		
		if (x1 == x2) {
			System.out.println("You must provide ZBRAC with an interval");
			return out;
		}
		
		double f1 = this.func.evaluate(x1);
		double f2 = this.func.evaluate(x2);
		out = true;
		
		for (int i = 0; i < ntry; i++) {
			
			if ((f1*f2) < 0.0) {
				this.bracket1 = x1;
				this.bracket2 = x2;
				return out;
			}
			
			if (Math.abs(f1) < Math.abs(f2)) {
				x1 = x1 + factor*(x1-x2);
				f1 = this.func.evaluate(x1);
			}
			else {
				x2 = x2 + factor*(x2-x1);
				f2 = this.func.evaluate(x2);
			}
		}
		
//		System.out.println("ZBRAC: Failed to find a bracket in "+ntry+" attempts");
		out = false;
		return out;
	}
	
	/** Find pairs of brackets inside the interval (x1, x2). From Numerical Recipes, chap 9.1
	 * To obtain the brackets, get this.xb1[] and this.xb2[]
	 * @param x1 lower value
	 * @param x2 upper value
	 * @param nint number of intervals
	 * @param nrt number of roots
	 * @return number of brackets found
	 */			
	public int zbrak(double x1, double x2, int nint, int nrt){
		int ntry = 50;
				
		if (x1 == x2) {
			System.out.println("You must provide ZBRAC with an interval");
			return 0;
		}
		
		int nbb = 0;
		
		double x = x1;
		double dx = (x2 - x1)/nint;
		
		this.xb1 = new double[nrt];
		this.xb2 = new double[nrt];
		
		double fp = this.func.evaluate(x);
		
		for (int i = 0; i < nint; i++) {			
			x = x + dx;
			double fc = this.func.evaluate(x);
			if ( fc*fp <= 0.0 ) {
				xb1[nbb] = x - dx;
				xb2[nbb] = x;
				nbb = nbb + 1;
				if (nbb == nrt) return nbb;
			}
			fp = fc;
		}
		
//		if (nbb == 0){
//			System.out.println("ZBRAK: Failed to find any brackets");
//		}
		return nbb;
	}
			

	/** Find the solution using RegulaFalsi.
	 * @param x1 lower limit on x.
	 * @param x2 upper limit on x.
	 * @param dxmin
	 * @return
	 */
	public double regulaFalsi(double x1, double x2) {
		double xlow;
		double xhigh;
		double del;
		double out = 0.0;
		double f;

		double fl = this.func.evaluate(x1);
		double fh = this.func.evaluate(x2);

		double test = fl * fh;

		if (test > 0.0) {
			if (fl == 0.0) return x1;
			if (fh == 0.0) return x2;
			err++;
			System.out.println("Root must be bracketed in ZeroFinder "+err);
		}

		if (fl < 0.0) {
			xlow = x1;
			xhigh = x2;
		} else {
			xlow = x2;
			xhigh = x1;
			double temp = fl;
			fl = fh;
			fh = temp;
		}

		double dx = xhigh - xlow;

		for (int i = 1; i < this.maxit; i++) {
			out = xlow + dx * fl / (fl - fh);
			f = this.func.evaluate(out);

			if (f < 0.0) {
				del = xlow - out;
				xlow = out;
				fl = f;
			} else {
				del = xhigh - out;
				xhigh = out;
				fh = f;
			}

			dx = xhigh - xlow;

			if ((Math.abs(del) < this.dxmin)
				|| (Math.abs(f) < this.accuracy)) {
				return out;
			}
		}

		System.out.println(
			" Regula Falsi exceeded " + this.maxit + " iterations ");
		return out;
	}

	/** Find the solution using the secant method.
	 * @param x1 lower limit on x.
	 * @param x2 upper limit on x.
	 * @return
	 */
	public double secant(double x1, double x2) {
		double xlow;
		double xhigh;
		double del;
		double out = 0.0;
		double f;

		double fl = this.func.evaluate(x1);
		double fh = this.func.evaluate(x2);

		if (Math.abs(fl) < Math.abs(fh)) {
			xlow = x2;
			xhigh = x1;
			double temp = fl;
			fl = fh;
			fh = temp;
		} else {
			xlow = x1;
			xhigh = x2;
		}

		for (int i = 1; i < this.maxit; i++) {
			double dx = (xlow - xhigh) * fh / (fh - fl);
			xlow = xhigh;
			fl = fh;
			xhigh = xhigh + dx;
			fh = this.func.evaluate(xhigh);

			if ((Math.abs(dx) < this.dxmin)
				|| (Math.abs(fh) < this.accuracy)) {
				return xhigh;
			}
		}

		System.out.println(
			" Secant Method exceeded " + this.maxit + " iterations ");
		return 0.0;
	}

	/** Find the solution using Ridder's method.
	 * @param x1 lower limit on x.
	 * @param x2 upper limit on x.
	 * @return
	 */
	public double ridder(double x1, double x2) {
		double fl = this.func.evaluate(x1);
		double fh = this.func.evaluate(x2);
		if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
			double xl = x1;
			double xh = x2;
			double ans = UNUSED;
			for (int j = 1; j <= this.maxit; j++) {
				double xm = 0.5 * (xl + xh);
				double fm = this.func.evaluate(xm);
				double s = Math.sqrt(fm * fm - fl * fh);
				if (s == 0.0)
					return ans;
				double xnew =
					xm + (xm - xl) * (MathUtils.sign(1.0, (fl - fh)) * fm / s);
				if (Math.abs(xnew - ans) <= this.dxmin)
					return ans;
				ans = xnew;
				double fnew = this.func.evaluate(ans);
				if (fnew == 0.0)
					return ans;
				if (MathUtils.sign(fm, fnew) != fm) {
					xl = xm;
					fl = fm;
					xh = ans;
					fh = fnew;
				} else if (MathUtils.sign(fl, fnew) != fl) {
					xh = ans;
					fh = fnew;
				} else if (MathUtils.sign(fh, fnew) != fh) {
					xl = ans;
					fl = fnew;
				} else
					System.out.println("never get here.");
				if (Math.abs(xh - xl) <= this.dxmin)
					return ans;
			}
			System.out.println("Ridder exceed maximum iterations");
		} else {
			if (fl == 0.0)
				return x1;
			if (fh == 0.0)
				return x2;
			err++;
			System.out.println("Root must be bracketed in Ridder "+err);
		}
		return 0.0;
	}
	
	/** Find the solution using Fixed Point Iteration method.
	 * @param x0 initial guess for x.
	 * @param dum not used.
	 * @return
	 */
	public double fixedPtIteration(double x0, double dum) {
		double xr = x0;
		
		for (int i = 1; i < this.maxit; i++) {
			double xrold = xr;
			double f = this.func.evaluate(xrold);
			xr = f + xrold;
			double dx = Math.abs(xr - xrold); 
			if ((Math.abs(dx) < this.dxmin)|| (Math.abs(f) < this.accuracy)) {
				return xr;
			}
		}
		System.out.println(
			" Fixed Point Iteration Method exceeded " + this.maxit + " iterations ");
		return 0.0;
	}

}
