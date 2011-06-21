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

package jat.demo.IntegratorExample;
import jat.plot.*;
import jat.alg.integrators.*;

/**
 * <P>
 * The SimpleIntegrator Class provides an example of how the RungeKutta8 class can be used to
 * integrate a set of ODEs and plot the solution.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RKF78Integrator implements Derivatives, Printable {

    TwoPlots traj_plot = new TwoPlots();

    /** Creates a new instance of SimpleIntegrator */
    public RKF78Integrator() {
        // set up the trajectory plot
        traj_plot.setTitle("Simple Oscillator");
        traj_plot.topPlot.setXLabel("t");
        traj_plot.topPlot.setYLabel("x");
        traj_plot.bottomPlot.setXLabel("t");
        traj_plot.bottomPlot.setYLabel("y");

    }
    /** Compute the derivatives.
     * @param t    double containing time or the independent variable.
     * @param x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
    public double[] derivs(double t, double[] x)
    {
        double [] out = new double[2];
        out[0] = x[1];
        out[1] = -x[0];
        return out;
    }

    /** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    public void print(double t, double [] y){

        // handle the first variable for plotting - this is a little mystery but it works
        boolean first = true;
        if (t == 0.0) first = false;

        // print to the screen for warm fuzzy feeling
        System.out.println(t+" "+y[0]+" "+y[1]+" "+first);

        // add data point to the plot
        traj_plot.topPlot.addPoint(0, t, y[0], first);
        traj_plot.bottomPlot.addPoint(0, t, y[1], first);

    }

    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args)
    {

        // create an RungeKuttaFehlberg78 integrator
        RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78();
        rk78.setNonAdaptive();
        rk78.setStepSize(0.01);

        // create an instance
        RKF78Integrator si = new RKF78Integrator();
        

        // initialize the variables
        double [] x0 = new double[2];
        x0[0] = 1.0;
        x0[1] = 0.0;

        // set the final time
        double tf = 10.0;

        // set the initial time to zero
        double t0 = 0.0;

        // integrate the equations
        rk78.integrate(t0, x0, tf, si, si, true);

        // make the plot visible
        si.traj_plot.setVisible(true);
    }
}
