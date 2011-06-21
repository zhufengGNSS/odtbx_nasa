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

public class RungeKutta8 {

    static final double [] a = { 0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 0.5,
    5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0 };

    static final double [] c = {41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0.0, 0.0};

    static final double [] chat = {0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0};

    static final double [][] b = new double [13][12];
    static {
        for (int i = 0; i < 13; i++) {
            for (int j = 0; j < 12; j++) {
                b[i][j] = 0.0;
            }
        }

        b[1][0] = 2.0/27.0;
        b[2][0] = 1.0/36.0;
        b[2][1] = 1.0/12.0;
        b[3][0] = 1.0/24.0;
        b[3][2] = 1.0/8.0;
        b[4][0] = 5.0/12.0;
        b[4][2] = -25.0/16.0;
        b[4][3] = 25.0/16.0;
        b[5][0] = 1.0/20.0;
        b[5][3] = 0.25;
        b[5][4] = 0.2;
        b[6][0] = -25.0/108.0;
        b[6][3] = 125.0/108.0;
        b[6][4] = -65.0/27.0;
        b[6][5] = 125.0/54.0;
        b[7][0] = 31.0/300.0;
        b[7][4] = 61.0/225.0;
        b[7][5] = -2.0/9.0;
        b[7][6] = 13.0/900.0;
        b[8][0] = 2.0;
        b[8][3] = -53.0/6.0;
        b[8][4] = 704.0/45.0;
        b[8][5] = -107.0/9.0;
        b[8][6] = 67.0/90.0;
        b[8][7] = 3.0;
        b[9][0] = -91.0/108.0;
        b[9][3] = 23.0/108.0;
        b[9][4] = -976.0/135.0;
        b[9][5] = 311.0/54.0;
        b[9][6] = -19.0/60.0;
        b[9][7] = 17.0/6.0;
        b[9][8] = -1.0/12.0;
        b[10][0] = 2383.0/4100.0;
        b[10][3] = -341.0/164.0;
        b[10][4] = 4496.0/1025.0;
        b[10][5] = -301.0/82.0;
        b[10][6] = 2133.0/4100.0;
        b[10][7] = 45.0/82.0;
        b[10][8] = 45.0/164.0;
        b[10][9] = 18.0/41.0;
        b[11][0] = 3.0/205.0;
        b[11][5] = -6.0/41.0;
        b[11][6] = -3.0/205.0;
        b[11][7] = -3.0/41.0;
        b[11][8] = 3.0/41.0;
        b[11][9] = 6.0/41.0;
        b[12][0] = -1777.0/4100.0;
        b[12][3] = -341.0/164.0;
        b[12][4] = 4496.0/1025.0;
        b[12][5] = -289.0/82.0;
        b[12][6] = 2193.0/4100.0;
        b[12][7] = 51.0/82.0;
        b[12][8] = 33.0/164.0;
        b[12][9] = 12.0/41.0;
        b[12][11] = 1.0;
    }

    /** Step size.
     */
    public double step_size;

    /** Default constructor
     */
    public RungeKutta8() {
        this.step_size = 1.0;
    }

    /** Construct RungeKutta8 integrator with specified step size.
     * @param s Step size.
     */

    public RungeKutta8(double s) {
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
     * @param derivs   Object containing the Equations of Motion
     * @return      double[] containing the new state
     */
    public double[] step(double x, double[] y, Derivatives derivs) {
        int n = y.length;
        double f[][] = new double[13][n];

        double y8th[] = new double[n];
        double ytmp[] = new double[n];

        double h = this.step_size;

        double xeval[] = new double [13];

        for (int i = 0; i < 13; i++)  // find times for function evals
        {
            xeval[i] = x + a[i] * h;
        }

        // build f matrix
        f[0] = derivs.derivs(x, y);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * b[1][0] * f[0][i];
        }
        f[1] = derivs.derivs(xeval[1], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[2][0]*f[0][i] + b[2][1]*f[1][i]);
        }
        f[2] = derivs.derivs(xeval[2], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[3][0]*f[0][i] + b[3][2]*f[2][i]);
        }
        f[3] = derivs.derivs(xeval[3], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[4][0]*f[0][i] + b[4][2]*f[2][i] + b[4][3]*f[3][i]);
        }
        f[4] = derivs.derivs(xeval[4], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[5][0]*f[0][i] + b[5][3]*f[3][i] + b[5][4]*f[4][i]);
        }
        f[5] = derivs.derivs(xeval[5], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[6][0]*f[0][i] + b[6][3]*f[3][i] + b[6][4]*f[4][i] + b[6][5]*f[5][i]);
        }
        f[6] = derivs.derivs(xeval[6], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[7][0]*f[0][i] + b[7][4]*f[4][i] + b[7][5]*f[5][i] + b[7][6]*f[6][i]);
        }
        f[7] = derivs.derivs(xeval[7], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[8][0]*f[0][i] + b[8][3]*f[3][i] + b[8][4]*f[4][i] + b[8][5]*f[5][i] + b[8][6]*f[6][i] + b[8][7]*f[7][i]);
        }
        f[8] = derivs.derivs(xeval[8], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[9][0]*f[0][i] + b[9][3]*f[3][i] + b[9][4]*f[4][i] + b[9][5]*f[5][i] + b[9][6]*f[6][i] + b[9][7]*f[7][i]
            + b[9][8]*f[8][i]);
        }
        f[9] = derivs.derivs(xeval[9], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[10][0]*f[0][i] + b[10][3]*f[3][i] + b[10][4]*f[4][i] + b[10][5]*f[5][i] + b[10][6]*f[6][i] + b[10][7]*f[7][i]
            + b[10][8]*f[8][i] + b[10][9]*f[9][i]);
        }
        f[10] = derivs.derivs(xeval[10], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[11][0]*f[0][i] + b[11][5]*f[5][i] + b[11][6]*f[6][i] + b[11][7]*f[7][i]
            + b[11][8]*f[8][i] + b[11][9]*f[9][i]);
        }
        f[11] = derivs.derivs(xeval[11], ytmp);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * (b[12][0]*f[0][i] + b[12][3]*f[3][i] + b[12][4]*f[4][i] + b[12][5]*f[5][i] + b[12][6]*f[6][i] + b[12][7]*f[7][i]
            + b[12][8]*f[8][i] + b[12][9]*f[9][i] + f[11][i]);
        }
        f[12] = derivs.derivs(xeval[12], ytmp);

        // construct solutions
        for (int i = 0; i < n; i++) {
            y8th[i] = y[i] + h*(chat[5]*f[5][i] + chat[6]*f[6][i] + chat[7]*f[7][i] + chat[8]*f[8][i] + chat[9]*f[9][i] + chat[11]*f[11][i] + chat[12]*f[12][i]);
        }

        return y8th;
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
