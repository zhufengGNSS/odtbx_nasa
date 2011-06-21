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

/**
 * <P>
 * RungeKutta8 is a fixed stepsize Runge-Kutta 8th order integrator.
 * Integrator parameter values came from:
 * Fehlberg, Erwin, NASA TR R-287.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RungeKutta4 {

    /** Step size.
     */
    public double step_size;

    /** Default constructor
     */
    public RungeKutta4() {
        this.step_size = 1.0;
    }

    /** Construct RungeKutta8 integrator with specified step size.
     * @param s Step size.
     */

    public RungeKutta4(double s) {
        this.step_size = s;
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
    public double[] step(double x, double[] y, Derivatives derivs) {
        int n = y.length;
        double f[][] = new double[4][n];

        double yout[] = new double[n];
        double ytmp[] = new double[n];

        double h = this.step_size;
        
        double hh = h * 0.5;
        double h6 = h/6.0;
        double xh = x + hh;

        // build f matrix
        f[0] = derivs.derivs(x, y);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[0][i];
        }
        f[1] = derivs.derivs(xh, ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[1][i];
        }
        f[2] = derivs.derivs(xh, ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * f[2][i];
        }
        f[3] = derivs.derivs(x + h, ytmp);


        // construct solutions
        for (int i = 0; i < n; i++) {
            yout[i] = y[i] + h6*(f[0][i] + 2.0 *(f[1][i] + f[2][i]) +f[3][i]);
        }

        return yout;
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

    public double[] integrate(double t0, double[] x0, double tf, Derivatives derivs, Printable pr, boolean print_switch) {
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

    public double[] integrate(double t0, double[] x0, double tf, Derivatives derivs) {
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

    public double[] integrate(double t0, double[] x0, double tf, EquationsOfMotion eom, boolean print_switch) {
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
            eom.print(t, oldstate);
        }

        // main integration loop

        while (t < tf) {
            newstate = step(t, oldstate, eom);
            for (int i = 0; i < neqns; i++) {
                oldstate[i] = newstate[i];
            }

            t = t + dt;

            if (print_switch) {
                eom.print(t, oldstate);
            }

            if ((t + dt) > tf) {
                dt = tf - t;
            }
        }
        return oldstate;
    }

    
}
