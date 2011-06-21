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
 * <p>This class demonstrates a way to do a simple s/c attitude simulation.
 * The simulated s/c is assumed to be rigid. 
 * An eccentric orbit with gravity gradient effect is the torque environment for 
 * the spacecraft. </p>
 * 
 * If the main()is assumed to be a client to JAT, the steps taken to do the 
 * simulation is as follows: <br>
 * 1.	Create an object of an integrator (RungeKutta8) <br>
 * 2.	Create an object of the class containing EOM <br>
 * 3.	Initialize variables <br>
 * 4.	Set the start & end time of the simulation <br>
 * 5.	Integrate the EOM <br>
 * 6.	Gets data at each time-step <br>
 * 7.	Make plots visible <br>
 * 
 * @author <a href="mailto:ntakada@users.sourceforge.net"> Noriko Takada
 * @version 1.4 
 * 
 * */

public class RigidGGEccentricOrbit implements EquationsOfMotion{

ThreePlots rotation_plot = new ThreePlots();
ThreePlots angle_plot = new ThreePlots();
SinglePlot quarternion_check = new SinglePlot();
TwoPlots quarternion_plot1 = new TwoPlots();
TwoPlots quarternion_plot2 = new TwoPlots();

public static final double PI = 3.14159;
//public static final double e = 0.5;

    /** Creates a new instance of SimpleIntegrator */
    public RigidGGEccentricOrbit()
    {
        // set up the trajectory plot
        rotation_plot.setTitle("Angular Velocity");
        rotation_plot.topPlot.setXLabel("# of orbits");
        rotation_plot.topPlot.setYLabel("w1(rad/sec)");
		rotation_plot.middlePlot.setXLabel("# of orbits");
		rotation_plot.middlePlot.setYLabel("w2(rad/sec)");
        rotation_plot.bottomPlot.setXLabel("# of orbits");
        rotation_plot.bottomPlot.setYLabel("w3(rad/sec)");
        
        quarternion_plot1.setTitle("Quaternions: e1 & e2");
        quarternion_plot1.topPlot.setXLabel("# of orbits");
        quarternion_plot1.topPlot.setYLabel("e1");
        quarternion_plot1.bottomPlot.setXLabel("# of orbits");
        quarternion_plot1.bottomPlot.setYLabel("e2");
        
        quarternion_plot2.setTitle("Quaternions: e3 & e4");
        quarternion_plot2.topPlot.setXLabel("# of orbits");
        quarternion_plot2.topPlot.setYLabel("e3");
        quarternion_plot2.bottomPlot.setXLabel("# of orbits");
        quarternion_plot2.bottomPlot.setYLabel("e4");
        
        angle_plot.setTitle("Angles between B and A frame");
        angle_plot.topPlot.setXLabel("# of orbits");
        angle_plot.topPlot.setYLabel("angle between a1 & b1 (rad)");
		angle_plot.middlePlot.setXLabel("# of orbits");
		angle_plot.middlePlot.setYLabel("angle between a2 & b2 (rad)");
        angle_plot.bottomPlot.setXLabel("# of orbits");
        angle_plot.bottomPlot.setYLabel("angle between a3 & b3 (rad)");
        
        quarternion_check.setTitle("Quarternion Check");
        quarternion_check.plot.setXLabel("# of orbits");
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
        
        double I1 = 10.42;
        double I2 = 35.42;
        double I3 = 41.67;
        double e  = 0.3;
        
        
        double c11 = 1- 2*( x[4]*x[4] + x[5]*x[5]);
        double c21 = 2* (x[3]*x[4]-x[5]*x[6]);
        double c31 = 2* (x[3]*x[5]+x[4]*x[6]); 
        
        /* The non-dimensionalized equations of motion for this series are:

         // *** Angular Velocity Equations ***
         w1dot = (2*pi*A)*((iyy-izz)/ixx)*(w2*w3-3*c21*c31*B);
         w2dot = (2*pi*A)*((izz-ixx)/iyy)*(w1*w3-3*c31*c11*B);
         w3dot = (2*pi*A)*((ixx-iyy)/izz)*(w1*w2-3*c11*c21*B);

         // *** Quaternions ***
         e1dot = (-pi*A*(-(w3+1/A)*e2 + w2*e3 - w1*e4));
         e2dot = (-pi*A*((w3+1/A)*e1 - w1*e3 - w2*e4));
         e3dot = (-pi*A*(-w2*e1 + w1*e2 - (w3-1/A)*e4));
         e4dot = (-pi*A*(w1*e1 + w2*e2 + (w3-1/A)*e3));

         c11 = 1 - 2*(e2^2 + e3^2)
         c21 = 2*(e1*e2 - e3*e4)
         c31 = 2*(e1*e3 + e2*e4)
			
		 A = ((1-e^2)^(3/2)/(1+e*cos(2*pi*t))^2
		 B = (1+e*cos(2*pi*t))^3/(1-e^2)^3	
         */

         double A=Math.pow((1-e*e),1.5)/Math.pow((1+ e*Math.cos(2*PI*t)),2);
         double B=Math.pow((1+e*Math.cos(2*PI*t)),3)/Math.pow((1-e*e), 3);
         
         double [] out = new double[7];
         // *** Angular Velocity Equations (dw/dnu (orbits)) *** //
         out[0] = A*2*PI*((I2-I3)/I1)*(x[1]*x[2]-3*c21*c31*B);
         out[1] = A*2*PI*((I3-I1)/I2)*(x[0]*x[2]-3*c31*c11*B);
         out[2] = A*2*PI*((I1-I2)/I3)*(x[0]*x[1]-3*c11*c21*B);
         
         // *** Quaternions (de/dnu (orbits)) ***
         out[3] = -PI*A*(-(x[2]+1/A)*x[4] + x[1]*x[5] - x[0]*x[6]);
         out[4] = -PI*A*((x[2]+1/A)*x[3] - x[0]*x[5] - x[1]*x[6]);
         out[5] = -PI*A*(-x[1]*x[3] + x[0]*x[4] - (x[2]-1/A)*x[6]);
         out[6] = -PI*A*(x[0]*x[3] + x[1]*x[4] + (x[2]-1/A)*x[5]);
         //System.out.println("e is "+ e+ " "+ x[0]+" "+A+" " +B);
              
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
        //System.out.println(t+" "+y[0]+" "+y[1]+" "+y[2]);
        
    }

    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args)
    {

        // create an RK8 integrator with step-size of 0.1
        RungeKutta8 rk8 = new RungeKutta8(0.1);

        // create an instance
        RigidGGEccentricOrbit si = new RigidGGEccentricOrbit();

        // initialize the variables
        double [] x0 = new double[7];
        x0[0] = 0.0;
        x0[1] = 0.0;
		x0[2] = 1.0;
		x0[3] = 0.0;
		x0[4] = 0.0;
		x0[5] = 0.0;
		x0[6] = 1.0;
		
        // set the final time
        double tf = 10.0;

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
