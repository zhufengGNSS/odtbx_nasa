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

package jat.alg.integrators;
import jat.cm.*;
import jat.matvec.data.VectorN;

/**
 * <P>
 * RungeKutta8 is a fixed stepsize Runge-Kutta 8th order integrator.
 * Integrator parameter values came from:
 * Fehlberg, Erwin, NASA TR R-287.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class Nystrom4 {

    /** Step size.
     */
    public double step_size;
    
    private LinePrinter _lp = null;;

    /** Default constructor
     */
    public Nystrom4() {
        this.step_size = 1.0;
        this._lp = new LinePrinter();
    }

    /** Construct RungeKutta8 integrator with specified step size.
     * @param s Step size.
     */

    public Nystrom4(double s) {
        this.step_size = s;
        this._lp = new LinePrinter();        
    }

    /** Construct RungeKutta8 integrator with specified step size.
     * @param s Step size.
     */

    public Nystrom4(double s, LinePrinter lp) {
        this.step_size = s;
        this._lp = lp;        
    }    

    /** Set step size
     * @param s     step size
     */

    public void setStepSize(double s) {
        this.step_size = s;
    }

    /** Take a single integration step. Note: it is assumed that any printing
     * will be done by the calling function.
     * @param x     time or independent variable
     * @param y     double[] containing needed inputs (usually the state)
     * @param neqns number of equations to integrate
     * @param derivs   Object containing the Equations of Motion
     * @return      double[] containing the new state
     */
    public double[] step(double t, double[] x, SecondDerivatives derivs) {
        int n = x.length/2;
        double k[][] = new double[3][n];
        double [] y = new double[n];
        double [] ydot = new double[n];

        double out[] = new double[2*n];
       
        // split out state and first derivatives
        for (int i = 0; i < n; i++) {
        	y[i] = x[i];
        	ydot[i] = x[i+n];
        }

        double ytmp[] = new double[n];
        double ydtmp[] = new double[n];

        double h = this.step_size;
        
        double hh = h * 0.5;
        double h6 = h/6.0;
        double h8 = h*h/8.0;

        // build f matrix
        k[0] = derivs.derivs(t, y, ydot);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * ydot[i] + h8*k[0][i];
            ydtmp[i] = ydot[i] + hh*k[0][i];
        }
        k[1] = derivs.derivs(t + hh, ytmp, ydtmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h*ydot[i] + hh*h * k[1][i];
            ydtmp[i] = ydot[i] + h*k[1][i];
        }
        k[2] = derivs.derivs(t + h, ytmp, ydtmp);

        // construct solutions
        for (int i = 0; i < n; i++) {
            out[i] = y[i] + h*ydot[i] + h*h6*(k[0][i] + 2.0 *k[1][i]);
            out[i+n] = ydot[i] + h6*(k[0][i] + 4.0*k[1][i] + k[2][i]);
        }
        
        return out;
    }

    /** Integrate from t0 to tf.
     * @param t0    initial time or independent variable
     * @param x0    double[] containing the initial state.
     * @param tf    final time or independent variable.
     * @param derivs   Object containing the Equations of Motion
     * @param pr    Object containing the print method.
     * @param print_switch  Boolean (true = print, false = don't print).
     * @return      double[] containing the new state
     */
    public double[] integrate(double t0, double[] x0, double tf, SecondDerivatives derivs, Printable pr, boolean print_switch) {
        int neqns = x0.length;

        double dt = this.step_size;
        double t = t0;
        double[] newstate = new double[neqns];
        double[] oldstate = new double[neqns];

        // put initial conditions into the previous state array

        for (int i = 0; i < neqns; i++) {
            oldstate[i] = x0[i];
        }

        if ((t + dt) > tf) {
            dt = tf - t;
        }

        if (print_switch) {
            pr.print(t, oldstate);
        }

        // main integration loop

        while (t < tf) {
            newstate = step(t, oldstate, derivs);
            for (int i = 0; i < neqns; i++) {
                oldstate[i] = newstate[i];
            }

            t = t + dt;

            if (print_switch) {
                pr.print(t, oldstate);
            }

            if ((t + dt) > tf) {
                dt = tf - t;
            }
        }
        return oldstate;
    }

    /** Integrate from t0 to tf. No printing/plotting interface provided.
     * @param t0    initial time or independent variable
     * @param x0    double[] containing the initial state.
     * @param tf    final time or independent variable.
     * @param derivs   Object containing the Equations of Motion
     * @return      double[] containing the new state
     */
    public double[] integrate(double t0, double[] x0, double tf, SecondDerivatives derivs) {
        int neqns = x0.length;

        double dt = this.step_size;
        double t = t0;
        double[] newstate = new double[neqns];
        double[] oldstate = new double[neqns];

        // put initial conditions into the previous state array

        for (int i = 0; i < neqns; i++) {
            oldstate[i] = x0[i];
        }

        if ((t + dt) > tf) {
            dt = tf - t;
        }

        // main integration loop

        while (t < tf) {
            newstate = step(t, oldstate, derivs);
            for (int i = 0; i < neqns; i++) {
                oldstate[i] = newstate[i];
            }

            t = t + dt;

            if ((t + dt) > tf) {
                dt = tf - t;
            }
        }
        return oldstate;
    }
    
    /** Integrate from t0 to tf.
     * @param t0    initial time or independent variable
     * @param x0    double[] containing the initial state.
     * @param tf    final time or independent variable.
     * @param eom   Object implementing the EquationsOfMotion interface
     * @param print_switch  Boolean (true = print, false = don't print).
     * @return      double[] containing the new state
     */
    public double[] integrate(double t0, double[] x0, double tf, SecondDerivatives eom, boolean print_switch) {
        int neqns = x0.length;

        double dt = this.step_size;
        double t = t0;
        double[] newstate = new double[neqns];
        double[] oldstate = new double[neqns];

        // put initial conditions into the previous state array

        for (int i = 0; i < neqns; i++) {
            oldstate[i] = x0[i];
        }

        if ((t + dt) > tf) {
            dt = tf - t;
        }

        if (print_switch) {
            _lp.print(t, oldstate);
        }

        // main integration loop

        while (t < tf) {
            newstate = step(t, oldstate, eom);
            for (int i = 0; i < neqns; i++) {
                oldstate[i] = newstate[i];
            }

            t = t + dt;

            if (print_switch) {
                _lp.print(t, oldstate);
            }

            if ((t + dt) > tf) {
                dt = tf - t;
            }
        }
        return oldstate;
    }
	/** Method for testing the class.
	 * @param args Input arguments (not used).
	 */

	public static void main(String args[])
	{
		Nystrom4 nystrom = new Nystrom4(10.0);
		double t0 = 0.0;
		// 90 minute orbit
		double tf = 90.0 * 60.0;
		double mu = 398600.4415; // GM in km^3/s^2 (value from JGM-3
		double a = Math.pow((mu*(tf/(2.0*Math.PI))*(tf/(2.0*Math.PI))), 1.0/3.0);
		TwoBody tb = new TwoBody(mu, a, 0.0, 0.0, 0.0, 0.0, 0.0);
		double[] x0 = tb.rv.getArray();
		double[] xf = nystrom.integrate(t0, x0, tf, tb, false);
		VectorN xinit = new VectorN(x0);
		VectorN xfin = new VectorN(xf);
		VectorN diff = xfin.minus(xinit);
		diff.print("diff");
		
		tb.propagate(t0, tf);
		VectorN rf = tb.getR();
		VectorN vf = tb.getV();
		VectorN rvf = rf.append(vf);
		VectorN diff2 = xfin.minus(rvf);
		diff2.print("diff2");
	}
    
}
