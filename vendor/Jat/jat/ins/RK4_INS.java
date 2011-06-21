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

package jat.ins;

//import jat.matvec.data.*;
//import java.io.*;
import jat.alg.integrators.*;
//import jat.gps_ins.*;
import jat.gps_ins.Derivs;
import jat.gps_ins.EOM;

/**
 * <P>
 * RK4_INS is a fixed stepsize Runge-Kutta 84th order integrator.
 * Integrator parameter values came from:
 * Fehlberg, Erwin, NASA TR R-287. It is specially designed
 * for using INS measurements to provide the trajectory.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RK4_INS {

    /** Step size.
     */
    public double step_size;
    
	private int index = 0;
	
	private INS_MeasurementList insData;

    /** Default constructor
     */
    public RK4_INS(INS_MeasurementList ins) {
        this.step_size = 1.0;
        this.insData = ins;        
    }

    /** Construct RungeKutta8 integrator with specified step size.
     * @param s Step size.
     */

    public RK4_INS(double s, INS_MeasurementList ins) {
        this.step_size = s;
        this.insData = ins;
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
    public double[] step(double x, double[] y, Derivs derivs, INS_Measurement measl_1, INS_Measurement measl) {
        int n = y.length;
        double f[][] = new double[4][n];

        double yout[] = new double[n];
        double ytmp[] = new double[n];

        double h = this.step_size;
        
        double hh = h * 0.5;
        double h6 = h/6.0;
        double xh = x + hh;

        double xeval[] = new double [13];


        // build f matrix
        f[0] = derivs.derivs(x, y, measl_1, measl, 0);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[0][i];
        }
        f[1] = derivs.derivs(xh, ytmp, measl_1, measl, 1);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[1][i];
        }
        f[2] = derivs.derivs(xh, ytmp, measl_1, measl, 1);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * f[2][i];
        }
        f[3] = derivs.derivs(x + h, ytmp, measl_1, measl, 2);


        // construct solutions
        for (int i = 0; i < n; i++) {
            yout[i] = y[i] + h6*(f[0][i] + 2.0 *(f[1][i] + f[2][i]) +f[3][i]);
        }

        return yout;
    }


    /** Take a single integration step. Note: it is assumed that any printing
     * will be done by the calling function.
     * @param x     time or independent variable
     * @param y     double[] containing needed inputs (usually the state)
     * @param neqns number of equations to integrate
     * @param derivs   Object containing the Equations of Motion
     * @return      double[] containing the new state
     */
    public double[] step(double x, double[] y, Derivs derivs) {
    	
		// get the two INS measurements
		INS_Measurement insMeasl_1 = this.insData.get(index);
		index = index + 1;
		INS_Measurement insMeasl = this.insData.get(index);
		index = index + 1;
    	
        int n = y.length;
        double f[][] = new double[4][n];

        double yout[] = new double[n];
        double ytmp[] = new double[n];

        double h = this.step_size;
        
        double hh = h * 0.5;
        double h6 = h/6.0;
        double xh = x + hh;

        double xeval[] = new double [13];


        // build f matrix
        f[0] = derivs.derivs(x, y, insMeasl_1, insMeasl, 0);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[0][i];
        }
        f[1] = derivs.derivs(xh, ytmp, insMeasl_1, insMeasl, 1);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + hh * f[1][i];
        }
        f[2] = derivs.derivs(xh, ytmp, insMeasl_1, insMeasl, 1);

        for (int i = 0; i < n; i++) {
            ytmp[i] = y[i] + h * f[2][i];
        }
        f[3] = derivs.derivs(x + h, ytmp, insMeasl_1, insMeasl, 2);


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

    public double[] integrate(double t0, double[] x0, double tf, Derivs derivs, Printable pr, boolean print_switch) {
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
			// get the two INS measurements
			INS_Measurement insMeasl_1 = this.insData.get(index);
			index = index + 1;
			INS_Measurement insMeasl = this.insData.get(index);
			index = index + 1;
        	
            newstate = step(t, oldstate, derivs, insMeasl_1, insMeasl);
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

    public double[] integrate(double t0, double[] x0, double tf, Derivs derivs) {
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
			// get the two INS measurements
			INS_Measurement insMeasl_1 = this.insData.get(index);
			index = index + 1;
			INS_Measurement insMeasl = this.insData.get(index);
			index = index + 1;
        	
            newstate = step(t, oldstate, derivs, insMeasl_1, insMeasl);
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
     * @param eom   Object implementing the EOM interface
     * @param print_switch  Boolean (true = print, false = don't print).
     * @return      double[] containing the new state
     */

    public double[] integrate(double t0, double[] x0, double tf, EOM eom, boolean print_switch) {
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
			// get the two INS measurements
			INS_Measurement insMeasl_1 = this.insData.get(index);
			index = index + 1;
			INS_Measurement insMeasl = this.insData.get(index);
			index = index + 1;
        	
            newstate = step(t, oldstate, eom, insMeasl_1, insMeasl);
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
