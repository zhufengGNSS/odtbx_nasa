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

package jat.demo.AttitudeExample;

import jat.plot.*;
import jat.alg.integrators.*;

/**
 * <P>This class demonstrates a way to do a simple s/c attitude simulation.
 * The simulated spacecraft is assumed to be rigid with constant external torque.
 * </P>
 * If the main()is assumed to be a client to JAT, the steps taken to do the 
 * simulation is as follows:<br>
 * 1.	Create an object of an integrator (RungeKutta8)<br>
 * 2.	Create an object of the class containing EOM <br>
 * 3.	Initialize variables <br>
 * 4.	Set the start & end time of the simulation <br>
 * 5.	Integrate the EOM <br>
 * 6.	Gets data at each time-step <br>
 * 7.	Make plots visible <br>
 * 
 * @author <a href="mailto:ntakada@users.sourceforge.net"> Noriko Takada
 * @version 1.6 
 */
public class RigidConstantTorque implements EquationsOfMotion{
ThreePlots rotation_plot = new ThreePlots();
//TwoPlots quarternion_plot1 = new TwoPlots();
//TwoPlots quarternion_plot2 = new TwoPlots();
ThreePlots angle_plot = new ThreePlots();
SinglePlot quarternion_check = new SinglePlot();
TwoPlots quarternion_plot1 = new TwoPlots();
TwoPlots quarternion_plot2 = new TwoPlots();

    /** Creates a new instance of SimpleIntegrator */
    public RigidConstantTorque()
    {
        // set up the trajectory plot
        rotation_plot.setTitle("Angular Velocity");
        rotation_plot.topPlot.setXLabel("t(sec)");
        rotation_plot.topPlot.setYLabel("w1(rad/sec)");
		rotation_plot.middlePlot.setXLabel("t(sec)");
		rotation_plot.middlePlot.setYLabel("w2(rad/sec)");
        rotation_plot.bottomPlot.setXLabel("t(sec)");
        rotation_plot.bottomPlot.setYLabel("w3(rad/sec)");
        
        quarternion_plot1.setTitle("Quaternions: e1 & e2");
        quarternion_plot1.topPlot.setXLabel("t");
        quarternion_plot1.topPlot.setYLabel("e1");
        quarternion_plot1.bottomPlot.setXLabel("t");
        quarternion_plot1.bottomPlot.setYLabel("e2");
        
        quarternion_plot2.setTitle("Quaternions: e3 & e4");
        quarternion_plot2.topPlot.setXLabel("t");
        quarternion_plot2.topPlot.setYLabel("e1");
        quarternion_plot2.bottomPlot.setXLabel("t");
        quarternion_plot2.bottomPlot.setYLabel("e2");
        
        angle_plot.setTitle("Angles between B and A frame");
        angle_plot.topPlot.setXLabel("t(sec)");
        angle_plot.topPlot.setYLabel("angle between a1 & b1 (rad)");
		angle_plot.middlePlot.setXLabel("t(sec)");
		angle_plot.middlePlot.setYLabel("angle between a2 & b2 (rad)");
        angle_plot.bottomPlot.setXLabel("t(sec)");
        angle_plot.bottomPlot.setYLabel("angle between a3 & b3 (rad)");
        
        quarternion_check.setTitle("Quarternion Check");
        quarternion_check.plot.setXLabel("t(sec)");
        quarternion_check.plot.setYLabel("e1^2 + e2^2 + e3^2 + e4^2");
       
    }
    /** Compute the derivatives.
     * Equations of Motion
     * @param t    double containing time or the independent variable.
     * @param x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
    public double[] derivs(double t, double[] x)
    {
        double M1 = 1;           // External torque
        double M2 = 0;
        double M3 = 0; 
        double I1 = 10.42;       // Spacecraft Principle inertia
        double I2 = 35.42;
        double I3 = 41.67;
        
        double [] out = new double[7];
        out[0] = ((I2-I3)/I1)*x[1]*x[2] + M1/I1 ;
        out[1] = ((I3-I1)/I2)*x[0]*x[2] + M2/I2 ;
        out[2] = ((I1-I2)/I3)*x[0]*x[1] + M3/I3 ;
        out[3] = -0.5* (-x[2]*x[4] + x[1]*x[5] - x[0]*x[6]);
        out[4] = -0.5* (x[2]*x[3]  - x[0]*x[5] - x[1]*x[6]);
        out[5] = -0.5* (-x[2]*x[6] + x[0]*x[4] - x[1]*x[3]);
        out[6] = -0.5* (x[2]*x[5]  + x[1]*x[4] + x[0]*x[3]);
        
        
        return out;
    }

    /** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    public void print(double t, double [] y){

        
        boolean first = true;
        if (t == 0.0) first = false;
        
        double w1 = y[0];
        double w2 = y[1];
        double w3 = y[2];
        double e1 = y[3];
        double e2 = y[4];
        double e3 = y[5];
        double e4 = y[6];
        double quat_check = e1*e1+e2*e2+e3*e3+e4*e4;
        
        // Calculate Transformation matrix elements
        double c11 = 1- 2*(e2*e2 + e3*e3);
        double c12 = 2* (e1*e2 + e3*e4);
        double c13 = 2* (e1*e3 - e2*e4);
        double c21 = 2* (e2*e1 - e3*e4);
        double c22 = 1- 2*(e3*e3 + e1*e1);
        double c23 = 2* (e2*e3 + e1*e4);
        double c31 = 2* (e3*e1 + e2*e4);
        double c32 = 2* (e3*e2 - e1*e4);
        double c33 = 1- 2*(e1*e1 + e2*e2);
        double angle11 = Math.toDegrees(Math.acos(c11));
        double angle22 = Math.toDegrees(Math.acos(c22));
        double angle33 = Math.toDegrees(Math.acos(c33));
        
        // add data point to the plot
        rotation_plot.topPlot.addPoint(0, t, w1, first);
		rotation_plot.middlePlot.addPoint(0, t, w2, first);
        rotation_plot.bottomPlot.addPoint(0, t, w3, first);
        quarternion_check.plot.addPoint(0, t, quat_check, first);
        angle_plot.topPlot.addPoint(0, t, angle11, first);
        angle_plot.middlePlot.addPoint(0, t, angle22, first);
        angle_plot.bottomPlot.addPoint(0, t, angle33, first);
        quarternion_plot1.topPlot.addPoint(0, t, e1, first);
        quarternion_plot1.bottomPlot.addPoint(0,t,e2, first);
        quarternion_plot2.topPlot.addPoint(0,t,e3, first);
        quarternion_plot2.bottomPlot.addPoint(0,t,e4,first);
        

        // also print to the screen 
        System.out.println(t+" "+y[0]+" "+y[1]+" "+y[2]);
    }

    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args)
    {

        // create an RK8 integrator with step-size of 0.1
        RungeKutta8 rk8 = new RungeKutta8(0.1);

        // create an instance
        RigidConstantTorque si = new RigidConstantTorque();

        // initialize the variables
        double [] x0 = new double[7];
        x0[0] = 0.0;
        x0[1] = 0.1;
		x0[2] = 0.0;
		x0[3] = 0.0;
		x0[4] = 0.0;
		x0[5] = 0.0;
		x0[6] = 1.0;
		
        // set the final time
        double tf = 30.0;

        // set the initial time to zero
        double t0 = 0.0;

        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);

        // make the plot visible
        si.rotation_plot.setVisible(true);
        si.angle_plot.setVisible(true);
        si.quarternion_check.setVisible(true);
    	si.quarternion_plot1.setVisible(true);
    	si.quarternion_plot2.setVisible(true);
    
    }
		

}
