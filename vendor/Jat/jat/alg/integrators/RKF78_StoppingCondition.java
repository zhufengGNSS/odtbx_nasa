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

/** Implements a Runge-Kutta-Fehlberg adaptive step size integrator
 * from Numerical Recipes. Modified to RK78 from the original RK45 in NR.
 * RK78 values from Erwin Fehlberg, NASA TR R-287
 * gets derivs via Derivatives interface
 * 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RKF78_StoppingCondition extends RungeKuttaFehlberg78 implements Printable{

    private double minStepSize_;
    private double stepSize_;
    private double accuracy_;
    private double currentStepSize_;
    private boolean verbose;
    private boolean adaptive_;
    
    private double stoptime = 0.0;

    /** Default constructor.
     */
    public RKF78_StoppingCondition() {
    	super();
    }

    /** Construct a RungeKuttaFehlberg78 integrator with user specified accuracy.
     * @param accuracy Desired accuracy.
     */
    public RKF78_StoppingCondition(double accuracy) {
    	super(accuracy);
    }

    /** Get the stop time.
     * @return stop time.
     */
    public double stopTime(){
    	return this.stoptime;
    }

    private int[] nok = new int[1];
    private int[] nbad = new int[1];


    /** Integrate the equations of motion.
     * @param y Initial state values.
     * @param start Initial time.
     * @param end Final time.
     * @param dv Equations of Motion.
     * @param pr Printable.
     * @param print_switch Print flag. True = call the print method.
     * @return the final state.
     */
    public double[] integrate(double[] y, double start, double end,
    Derivatives dv, Printable pr, boolean print_switch,
    StoppingCondition sc, boolean sc_switch) {
    	
         double h = stepSize_;
        double hmin = minStepSize_;
        nok[0] = nbad[0] = 0;
        if (adaptive_) {
            odeint(y, start, end, accuracy_, h, hmin, nok, nbad, dv, pr, print_switch, sc, sc_switch);
            if (verbose) {
                System.out.println("nok = " + nok[0] + "\tnbad = " + nbad[0]);
            }
        } else {
            rkdumb(y, start, end, h, dv, pr, print_switch, sc, sc_switch);
        }
        return y;
    }

    /** Integrate the equations of motion. No printing/plotting interface provided.
     * @param y Initial state values.
     * @param start Initial time.
     * @param end Final time.
     * @param dv Equations of Motion.
     * @param pr Printable.
     * @param print_switch Print flag. True = call the print method.
     * @return the final state.
     */
    public double[] integrate(double[] y, double start, double end, Derivatives dv,
        StoppingCondition sc, boolean sc_switch) {
        boolean print_switch = false;
        double h = stepSize_;
        double hmin = minStepSize_;
        nok[0] = nbad[0] = 0;
        if (adaptive_) {
            odeint(y, start, end, accuracy_, h, hmin, nok, nbad, dv, this, print_switch, sc, sc_switch);
            if (verbose) {
                System.out.println("nok = " + nok[0] + "\tnbad = " + nbad[0]);
            }
        } else {
            rkdumb(y, start, end, h, dv, this, print_switch, sc, sc_switch);
        }
        return y;
    }

    private void rkdumb(double[] ystart, double start, double end,
    double h, Derivatives dv, Printable pr, boolean print_switch,
        StoppingCondition sc, boolean sc_switch) {
        	
        boolean stop = false;
        this.stoptime = start;
        int nvar = ystart.length;

        int nSteps = (int) Math.abs((end - start) / h);
        if (nSteps < 1)
            nSteps = 1;
        h = (end - start) / nSteps;
        double[] dydx = new double[nvar];
        double[] yend = new double[nvar];
        double[] yerr = new double[nvar];
        for (int step = 0; step < nSteps; step++) {
            double x = start + step * h;
            dydx = dv.derivs(x, ystart);
            rkck(ystart, dydx, x, h, yend, yerr, dv);
            for (int n = 0; n < nvar; n++){
                ystart[n] = yend[n];
            }
            
            // print
            if (print_switch) {
                pr.print(x, ystart);
            }
            
            // check stopping condition
            if (sc_switch){
            	stop = sc.stoppingCondition(x, ystart);
            }
            if (stop) {
            	this.stoptime = x;
            	return;
            }
            this.stoptime = x;
            	
        }

    }

    private static final int MAXSTP = 1000000;
    private static final double TINY = 1.0e-30;

    private int kmax;
    private int kount;
    private double[] xp;
    private double[][] yp;
    private double dxsav;

    private void odeint(double[] ystart, double x1, double x2,
    double eps, double h1, double hmin, int[] nok,
    int[] nbad, Derivatives dv, Printable pr, boolean print_switch,
    StoppingCondition sc, boolean sc_switch) {
    	boolean stop = false;
    	this.stoptime = x1;
        int nvar = ystart.length;
        double[] x = new double[1];
        double[] hnext = new double[1];
        double[] hdid = new double[1];
        double[] yscal = new double[nvar];
        double[] y = new double[nvar];
        double[] dydx = new double[nvar];
        x[0] = x1;
        double h = Math.abs(h1);
        if (x2 < x1)
            h = -h;
        nok[0] = nbad[0] = kount = 0;
        for (int i = 0; i < nvar; i++)
            y[i] = ystart[i];
        double xsav = 0;
        if (kmax > 0)
            xsav = x[0] - dxsav * 2.0;
        for (int nstp = 1; nstp <= MAXSTP; nstp++) {
            dydx = dv.derivs(x[0], y);
            for (int i = 0; i < nvar; i++)
                yscal[i] = Math.abs(y[i]) + Math.abs(dydx[i] * h) + TINY;
            if (kmax > 0 && kount < kmax - 1 &&
            Math.abs(x[0] - xsav) > Math.abs(dxsav)) {
                xp[++kount] = x[0];
                for (int i = 0; i < nvar; i++)
                    yp[i][kount] = y[i];
                xsav = x[0];
            }
            if ( (x[0] + h - x2) * (x[0] + h - x1) > 0.0)
                h = x2 - x[0];
            rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, dv);
            if (hdid[0] == h)
                ++nok[0];
            else
                ++nbad[0];

            // invoke the print method, if required
            if (print_switch) {
                pr.print(x[0], y);
            }
            
            // update stoptime
            this.stoptime = x[0];
            
            // check for stopping condition            
            if (sc_switch){
            	stop = sc.stoppingCondition(x[0], y);
            }
            if (stop) {
            	return;
            }
            

            if ( (x[0] - x2) * (x2 - x1) >= 0.0 ) {   // are we done?
                // save off the data
                for (int i = 0; i < nvar; i++){
                    ystart[i] = y[i];
                }
                if (kmax != 0) {
                    xp[++kount] = x[0];
                    for (int i = 0; i < nvar; i++)
                        yp[i][kount] = y[i];
                }
                return;
            }
            if (Math.abs(hnext[0]) <= hmin) {
                error("Step size too small in odeint");
                System.out.println("h = "+hnext[0]);
            }
            h = hnext[0];
            currentStepSize_ = h;            // added for comphys
            //            System.out.println("Current Step Size = "+h);
        }
        error("Too many steps in routine odeint");
        System.out.println("step size = "+currentStepSize_);
    }

    private static final double SAFETY = 0.9;
    private static final double PGROW = -1.0/8.0;
    private static final double PSHRNK = -1.0/7.0;
    private static final double ERRCON = 2.56578451395034701E-8;

    private static final double [] a = { 0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 0.5,
    5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0 };

    private static final double [][] b = new double [13][12];
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

    private static final double [] c = {41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0.0, 0.0};

    private static final double [] chat = {0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0};


    public void print(double t, double [] y){
        // do nothing
    }

    private void error(String msg) {
        System.err.println("RungeKuttaFehlberg78: " + msg);
    }

}

