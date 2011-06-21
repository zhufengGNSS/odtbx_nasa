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
 
package jat.attitude;
 
import jat.plot.*;
import jat.alg.integrators.*;

/**
 * @author Noriko Takada
 *
 * This class implements Printable.
 * It gets output at each integration step and make plots.
 * 
 */


public class FourRWManeuverOutput implements Printable
{
		// Create variables for the necessary plots
		ThreePlots rotation_plot;								
		ThreePlots angle_plot;
		SinglePlot quarternion_check;
		TwoPlots quarternion_plot1;
		TwoPlots quarternion_plot2;
		TwoPlots RWplot1;
		TwoPlots RWplot2;
		float quat_values[][];
		//double timeDuration;
		double time_step;
		int plotYes = 0;
		
		
		//constructor: No plot generation
		public FourRWManeuverOutput(double time_step, float quat_values[][])
		{
			this.quat_values = quat_values;
			this.time_step = time_step;
		}
		//constructor
		public FourRWManeuverOutput (ThreePlots rotation, ThreePlots angle
      	                       ,SinglePlot quatCheck, TwoPlots quat1
      	                       ,TwoPlots quat2, TwoPlots RW1, TwoPlots RW2
      	                       ,double time_step, float quat_values[][])
      	{
      		plotYes = 1;
      		rotation_plot = rotation;
      		angle_plot = angle;
      		quarternion_check = quatCheck;
      		quarternion_plot1 = quat1;
      		quarternion_plot2 = quat2;
      		RWplot1 = RW1;
      		RWplot2 = RW2;
      		
      		// Quaternion Values Index
      		// Index	[0  1  2  3  4]
      		// Variable	[t e1 e2 e3 e4]
      		// quat_values = new  float[5][numberOfPoints];
      		this.quat_values = quat_values;
      		//this.timeDuration = timeDuration;
      		this.time_step = time_step;      	
      		
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
        	
        	RWplot1.setTitle("RW speeds 1&2");
        	RWplot1.topPlot.setXLabel("t");
        	RWplot1.topPlot.setYLabel("Omega1");
        	RWplot1.bottomPlot.setXLabel("t");
        	RWplot1.bottomPlot.setYLabel("Omega2");
        
        	RWplot2.setTitle("RW speeds 3&4");
        	RWplot2.topPlot.setXLabel("t");
        	RWplot2.topPlot.setYLabel("Omega3");
        	RWplot2.bottomPlot.setXLabel("t");
        	RWplot2.bottomPlot.setYLabel("Omega4");
      	}// End of constructor	
		
		/** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
     	*  This method is executed by the propagator at each integration step.
     	* @param t Time.
     	* @param y Data array.
     	*/
    	public void print(double t, double [] y)
    	{
        	boolean first = true;
        	if (t == 0.0) first = false;
        	
        	int currentPts = (int)(t/time_step); // This is the quarternion array index.
        
       	 	double w1 = y[0];
        	double w2 = y[1];
        	double w3 = y[2];
        	double e1 = y[3];
        	double e2 = y[4];
        	double e3 = y[5];
        	double e4 = y[6];
        	double Omega1 = y[7];
        	double Omega2 = y[8];
        	double Omega3 = y[9];
        	double Omega4 = y[10];
        	
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
        
        	if (plotYes == 1)
        	{
        	
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
        		RWplot1.topPlot.addPoint(0, t, Omega1, first);
        		RWplot1.bottomPlot.addPoint(0,t,Omega2, first);
        		RWplot2.topPlot.addPoint(0,t,Omega3, first);
        		RWplot2.bottomPlot.addPoint(0,t,Omega4,first);
        	
        		
        	}
        	
			quat_values[0][currentPts] = (float)t; // time value
			quat_values[1][currentPts] = (float)e1; // quaternion 1
			quat_values[2][currentPts] = (float)e2; // quarternion 2
			quat_values[3][currentPts] = (float)e3; // quarternion 3
			quat_values[4][currentPts] = (float)e4; // quarternion 4
        	// also print to the screen 
        	System.out.println(t+" "+y[0]+" "+y[1]+" "+y[2]);
    	}// End of print
    	
}// End of File
    
	
	
	