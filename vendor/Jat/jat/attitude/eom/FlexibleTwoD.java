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
package jat.attitude.eom;

import jat.alg.integrators.*;
import jat.plot.*;
import jat.matvec.data.*;

/**
 * <p>
 * Equations of motion of a flexible spacecraft. 
 * The system being modeled consists of a rigid spacecraft
 * with two massless beams with a concentrated mass at the tip
 * of each beam,cantilevered at the central body on opposite sides.
 * 
 *
 * @author Noriko Takada
 * @version 1.3 (8/15/2004)
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion
 *  
 */
public class FlexibleTwoD implements EquationsOfMotion
{
	//		time_step:			Time step of numerical integration
 	//		quat_values[][]		Two dimensional array that contains quarternions from simulation
 	//		I1,I2,I3			Principal moments of inertia
 	//		m					Mass of the concentrated mass at the beam tips
 	//		a					Distance from the center of mass to the beam attachment point
 	//		L					Length of the beams
 	//		EI					Flexural Rigidity
 	//		rotation_plot		Plots of angular velocities
 	//		angle_plot			Plots of euler angles
 	//		quarternion_check	Plot of e1^2 + e2^2 + e3^2 +e4^2
 	//		quat_plot			Plots of quaternions
 	//		energy_and_angularM	Plots of system energy and angular momentum
 	// 
	double  time_step;
	private float quat_values[][];
	private float quatBeam1[][];
	private float quatBeam2[][];
	
	private double J = 13.35;
    private double m = 1;
    private double a = 0.40;
    private double L = 1.80;
    private double EI = 15.52;
    
    private double wn = 1;
    private double damping = 1;
    private double K = wn*wn*J;
    private double Kd =J*2*damping*wn;
    private double K_bang = K;
    private double Kd_bang = Kd;
    private double torque_level=1;
    private double no_torque_zone=0.001;
    private double alpha_com = 2*Math.PI/180;
    
    // Create variables for the necessary plots
	private ThreePlots positionPlot = new ThreePlots();								
	private ThreePlots velocityPlot = new ThreePlots();
	private FourPlots hPlot = new FourPlots();
	private ThreePlots energyPlot = new ThreePlots();
	private FourPlots quaternionPlot = new FourPlots();
	private TwoPlots	beamAnglePlot = new TwoPlots();
	
	/**
	 * Constructor
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 */
	public FlexibleTwoD(double time_step, float quat_values[][],float quatBeam1[][],
						float quatBeam2[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.quatBeam1 = quatBeam1;
		this.quatBeam2 = quatBeam2;
	}
	
	
	/** Construct a RConstantTorque object
	 * @param		time_step:			Time step of numerical integration
 	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
 	 * @param		I1,I2,I3			Principal moments of inertia
 	 * @param		m					Tip mass
 	 * @param		a					Distance from the center of mass to the beam attachment point
 	 * @param		L					Length of the massless beam
 	 * @param		EI					Flexural rigidity
	 */
	public FlexibleTwoD(double time_step, double m, double a, double L, double EI,
											  double J, double K, double Kd, double alpha_com,
											  float quat_values[][],
											  float quatBeam1[][], float quatBeam2[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.quatBeam1 = quatBeam1;
		this.quatBeam2 = quatBeam2;
		this.m = m;
		this.a = a;
		this.L = L;
		this.EI = EI;	
		this.J=J;
		
		this.K = K;
		this.Kd = Kd;
		this.alpha_com = alpha_com*Math.PI/180;
		
	}
	
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		/*
		private ThreePlots positionPlot = new ThreePlots();								
		private ThreePlots velocityPlot = new ThreePlots();
		private FourPlots hPlot = new FourPlots();
		private TwoPlots EandQPlot = new TwoPlots();
		private FourPlots quaternionPlot = new FourPlots();
		*/
		
		// Setup plots
		positionPlot.setTitle("Position plot");
		positionPlot.topPlot.setXLabel("t(sec)");
		positionPlot.topPlot.setYLabel("alpha (rad)");
		positionPlot.middlePlot.setXLabel("t (sec)");
		positionPlot.middlePlot.setYLabel("u1 (m)");
		positionPlot.bottomPlot.setXLabel("t (sec)");
		positionPlot.bottomPlot.setYLabel("u2 (m)");
		
		velocityPlot.setTitle("Velocity plot");
		velocityPlot.topPlot.setXLabel("t (sec)");
		velocityPlot.topPlot.setYLabel("alphaDot (rad/sec)");
		velocityPlot.middlePlot.setXLabel("t (sec)");
		velocityPlot.middlePlot.setYLabel("u1Dot (m/sec)");
		velocityPlot.bottomPlot.setXLabel("t (sec)");
		velocityPlot.bottomPlot.setYLabel("u2Dot (m/sec)");
		
		hPlot.setTitle("Angular momentum");
		hPlot.firstPlot.setXLabel("t (sec)");
		hPlot.firstPlot.setYLabel("Hm1");
		hPlot.secondPlot.setXLabel("t (sec)");
		hPlot.secondPlot.setYLabel("Hm2");
		hPlot.thirdPlot.setXLabel("t (sec)");
		hPlot.thirdPlot.setYLabel("Hsc");
		hPlot.fourthPlot.setXLabel("t (sec)");
		hPlot.fourthPlot.setYLabel("Total Angular Momentum");
		
		energyPlot.setTitle("Energy Plot");
		energyPlot.topPlot.setXLabel("t (sec)");
		energyPlot.topPlot.setYLabel("Potential");
		energyPlot.middlePlot.setXLabel("t (sec)");
		energyPlot.middlePlot.setYLabel("Kinetic");
		energyPlot.bottomPlot.setXLabel("t (sec)");
		energyPlot.bottomPlot.setYLabel("Total Energy");
		
		quaternionPlot.setTitle("Quaternions");
		quaternionPlot.firstPlot.setXLabel("t (sec)");
		quaternionPlot.firstPlot.setYLabel("q1");
		quaternionPlot.secondPlot.setXLabel("t (sec)");
		quaternionPlot.secondPlot.setYLabel("q2");
		quaternionPlot.thirdPlot.setXLabel("t (sec)");
		quaternionPlot.thirdPlot.setYLabel("q3");
		quaternionPlot.fourthPlot.setXLabel("t (sec)");
		quaternionPlot.fourthPlot.setYLabel("q4");
		
		beamAnglePlot.setTitle("Beam Angle plot");
		beamAnglePlot.topPlot.setXLabel("t (sec)");
		beamAnglePlot.topPlot.setYLabel("theta1 (degree)");
		beamAnglePlot.bottomPlot.setXLabel("t (sec)");
		beamAnglePlot.bottomPlot.setYLabel("theta2 (degree)");
		        
	}
	
	
	/** Compute the derivatives.
     * Equations of Motion
     * @param t    double containing time or the independent variable.
     * @param x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
		
	public double[] derivs(double t, double[] x)
    {
       	//! The method derivs has the limit of 65535 byts        
       	
       	double [] out = new double[6];
       	// Define state variables
       	double	alpha 		= x[0];
       	double u1 			= x[1];
       	double u2	 		= x[2];
       	double alphaDot	= x[3];
       	double u1Dot 		= x[4];
       	double u2Dot 		= x[5];
       	
       	// Define
       	double u1square = u1*u1;
       	double u2square = u2*u2;
       	double Lfifth = L*L*L*L*L;
       	double Lfour = L*L*L*L;
       	double Lcube = L*L*L;
       	double msquare = m*m;
       	double asquare = a*a;
       	double alphaDotSquare = alphaDot*alphaDot;
       	double Lsquare = L*L;
       	double Me = 0;
       	double Mc = K*(alpha_com-alpha)-Kd*alphaDot;;
       	
       	out[0]=alphaDot;
       	out[1]=u1Dot;
       	out[2]=u2Dot;
       	// alphaDoubleDot Equation
       	out[3] = (-a*m*u1*Lcube*alphaDotSquare+a*m*u2*Lcube*alphaDotSquare+3*a*EI*u1-3*a*EI*u2-2*alphaDot*m*u1*u1Dot*Lcube-2*alphaDot*m*u2*u2Dot*Lcube-m*u1*Lfour*alphaDotSquare+m*u2*Lfour*alphaDotSquare+Me*Lcube+Mc*Lcube+(3*EI*u1-3*EI*u2)*L)/(m*u2square*Lcube+m*u1square*Lcube+J*Lcube);
       	// u1DoubleDot Equation
       	out[4] = (-3*m*EI*u1*asquare+u1*J*m*Lcube*alphaDotSquare-m*a*Mc*Lcube-m*Mc*Lfour-m*Me*Lfour+u1*u2square*Lcube*msquare*alphaDotSquare-2*a*Lfour*msquare*u2*alphaDotSquare-Lcube*msquare*u2*alphaDotSquare*asquare+u1*u1square*Lcube*msquare*alphaDotSquare+Lcube*msquare*u1*alphaDotSquare*asquare-3*m*EI*u1*u1square-3*m*EI*u1*u2square+3*m*EI*u2*asquare-3*EI*u1*J+2*a*Lfour*msquare*u1*alphaDotSquare-Lfifth*u2*msquare*alphaDotSquare-m*a*Me*Lcube+2*a*alphaDot*Lcube*u1*u1Dot*msquare+Lfifth*u1*msquare*alphaDotSquare+2*a*alphaDot*Lcube*u2*u2Dot*msquare+2*alphaDot*Lfour*u1*u1Dot*msquare+2*alphaDot*Lfour*u2*u2Dot*msquare+(6*a*m*EI*u2-6*a*m*EI*u1)*L+(-3*m*EI*u1+3*m*EI*u2)*Lsquare)/(u2square*Lcube*msquare+u1square*Lcube*msquare+J*m*Lcube);
       	// u2DoubleDot Equation
       	out[5] = (-3*m*EI*u2*u2square-3*EI*u2*J+m*a*Mc*Lcube+u2*u2square*Lcube*msquare*alphaDotSquare+u2*u1square*Lcube*msquare*alphaDotSquare+m*Mc*Lfour+m*Me*Lfour-2*a*alphaDot*Lcube*u1*u1Dot*msquare+Lcube*msquare*u2*alphaDotSquare*asquare-Lcube*msquare*u1*alphaDotSquare*asquare-3*m*EI*u2*asquare+3*m*EI*u1*asquare-3*m*EI*u2*u1square-2*a*Lfour*msquare*u1*alphaDotSquare+2*a*Lfour*msquare*u2*alphaDotSquare-2*a*alphaDot*Lcube*u2*u2Dot*msquare+m*a*Me*Lcube-Lfifth*u1*msquare*alphaDotSquare+Lfifth*u2*msquare*alphaDotSquare-2*alphaDot*Lfour*u2*u2Dot*msquare-2*alphaDot*Lfour*u1*u1Dot*msquare+J*m*u2*Lcube*alphaDotSquare+(6*a*m*EI*u1-6*a*m*EI*u2)*L+(3*m*EI*u1-3*m*EI*u2)*Lsquare)/(u2square*Lcube*msquare+u1square*Lcube*msquare+J*m*Lcube);
       	
       	return out;
       	
     }// End of derivs
    	
    /** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    	public void print(double t, double [] y)
    	{
        	boolean first = true;
        	if (t == 0.0) first = false;
        	
        	int currentPts = (int)(t/time_step); // This is the array index
        	
       	 	double alpha = y[0];
        	double u1 = y[1];
        	double u2 = y[2];
        	double alphaDot = y[3];
        	double u1Dot = y[4];
        	double u2Dot = y[5];
        	
        	double theta2 = Math.asin(u2/L);
     
        		theta2=-theta2;
        	
        	double theta1 = Math.asin(u1/L);
        	// Extract quaternion values
        	// See Bong wie (p.318)
        	// Euler axis
        	double e1=0;
        	double e2=0;
        	double e3=1;
        	// Define Euler Parameters
        	// Make the alpha twice as much to exaggerate it
        	double q1=e1*Math.sin(alpha*2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q2=e2*Math.sin(alpha*2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q3=e3*Math.sin(alpha*2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q4=Math.cos(alpha/2);
        	
        	// Define Euler Parameters for Beam1 rotaion
        	double q1Beam1 = e1*Math.sin(theta1/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q2Beam1 = e2*Math.sin(theta1/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q3Beam1 = e3*Math.sin(theta1/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q4Beam1 = Math.cos(theta1/2);
        	
        	//Define Euler Parameters for Beam2 rotation
        	double q1Beam2 = e1*Math.sin(theta2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q2Beam2 = e2*Math.sin(theta2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q3Beam2 = e3*Math.sin(theta2/2)/Math.sqrt(e1*e1+e2*e2+e3*e3);
        	double q4Beam2 = Math.cos(theta2/2);
        	
        	// Calculate Transformation matrix elements
        	// Transform from A to B see Bong Wie (p.318)
        	double c11 = 1- 2*(q2*q2 + q3*q3);
        	double c12 = 2* (q1*q2 + q3*q4);
        	double c13 = 2* (q1*q3 - q2*q4);
        	double c21 = 2* (q2*q1 - q3*q4);
        	double c22 = 1- 2*(q3*q3 + q1*q1);
        	double c23 = 2* (q2*q3 + q1*q4);
        	double c31 = 2* (q3*q1 + q2*q4);
        	double c32 = 2* (q3*q2 - q1*q4);
        	double c33 = 1- 2*(q1*q1 + q2*q2);
        	double angle11 = Math.toDegrees(Math.acos(c11));
        	double angle22 = Math.toDegrees(Math.acos(c22));
        	double angle33 = Math.toDegrees(Math.acos(c33));
        	double quat_check = q1*q1 + q2*q2 + q3*q3 + q4*q4;
        	
        	        	// Build 3 x 4 Q matrix
    		double[][] array = new double[3][4];
    		// Row1
    		array[0][0] = c11;
    		array[0][1] = c12;
    		array[0][2] = c13;
    		// Row2
    		array[1][0] = c21;
    		array[1][1] = c22;
    		array[1][2] = c23;
    		// Row3
    		array[2][0] = c31;
    		array[2][1] = c32;
    		array[2][2] = c33;
    		    	
    		//Construct Transformation Matrix
    		Matrix  T = new Matrix(array,3,3);
    		Matrix  T_transpose = new Matrix(T.transpose().A, 3,3);
    		//Matrix  q_times_qtrans = new Matrix(q_matrix.times(q_transpose).A, 3,3);
    		//Matrix  inv_q_times_qtrans = new Matrix(q_times_qtrans.invert().A, 3,3);
    		//double[][] pseudoinv_q = q_transpose.times(inv_q_times_qtrans).A;
    		
    		//Energy Calculation
    		double Potential = 3*EI*(u1*u1)/(2*L*L*L) + 3*EI*(u2*u2)/(2*L*L*L);
    		double Kinetic = (0.5)*J*(alphaDot*alphaDot) +  (0.5)*m*((-u1*alphaDot)*(-u1*alphaDot) + ((a+L)*alphaDot + u1Dot)*((a+L)*alphaDot + u1Dot))  +  (0.5)*m*((-u2*alphaDot)*(-u2*alphaDot) + (-(a+L)*alphaDot + u2Dot)*(-(a+L)*alphaDot + u2Dot));
    		double energy= Kinetic+Potential;
			//double Potential = 3*EI*(u1*u1)/(2*L*L*L) + 3*EI*(u2*u2)/(2*L*L*L);
			//double Kinetic = (1/2)*J*(alphaDot*alphaDot)+(1/2)*m*((-u1*alphaDot)*(-u1*alphaDot) + ((a+L)*alphaDot+u1Dot)*((a+L)*alphaDot+u1Dot)) + (1/2)*m*((-u2*alphaDot)*(-u2*alphaDot)+(-(a+L)*alphaDot + u2Dot)* (-(a+L)*alphaDot + u2Dot));
			//double Energy = Potential + Kinetic;  

  			//Angular Momentum Calculation
  			//Hsc (i,1)= J*alpha_dot;
   			double Hsc = J*alphaDot;
   			
   			//r1a1 = (a+L)*cos(alpha) + u1*sin(alpha);
   			//r1a2 = -(a+L)*sin(alpha) + u1*cos(alpha);
   			//r1_dot_a1 = -u1*alpha_dot*cos(alpha) + (a+L)*alpha_dot*sin(alpha) + u1_dot*sin(alpha);
   			//r1_dot_a2 = u1*alpha_dot*sin(alpha) + (a+L)*alpha_dot*cos(alpha) + u1_dot*cos(alpha);
   			double r1a1 = (a+L)*Math.cos(alpha) + u1*Math.sin(alpha);
   			double r1a2 = -(a+L)*Math.sin(alpha) + u1*Math.cos(alpha);
   			double r1_dot_a1 = -u1*alphaDot*Math.cos(alpha) + (a+L)*alphaDot*Math.sin(alpha) + u1Dot*Math.sin(alpha);
   			double r1_dot_a2 = u1*alphaDot*Math.sin(alpha) + (a+L)*alphaDot*Math.cos(alpha) + u1Dot*Math.cos(alpha);
   			
   			//r2a1 = -(a+L)*cos(alpha) + u2*sin(alpha); 
   			//r2a2 = (a+L)*sin(alpha) + u2*cos(alpha);
   			//%r2_dot_a1 = alpha_dot*cos(alpha) -(a+L)*alpha_dot*sin(alpha) + u2_dot*sin(alpha);
   			//r2_dot_a1 = -alpha_dot*u2*cos(alpha) -(a+L)*alpha_dot*sin(alpha) + u2_dot*sin(alpha);
   			//%r2_dot_a2 = -alpha_dot*u2*sin(alpha) - (a+L)*alpha_dot*cos(alpha) + u2_dot*cos(alpha);
   			//r2_dot_a2 = alpha_dot*u2*sin(alpha) - (a+L)*alpha_dot*cos(alpha) + u2_dot*cos(alpha);
   
   			double r2a1 = -(a+L)*Math.cos(alpha) + u2*Math.sin(alpha); 
   			double r2a2 = (a+L)*Math.sin(alpha) + u2*Math.cos(alpha);
   			//%r2_dot_a1 = alpha_dot*cos(alpha) -(a+L)*alpha_dot*sin(alpha) + u2_dot*sin(alpha);
   			double r2_dot_a1 = -alphaDot*u2*Math.cos(alpha) -(a+L)*alphaDot*Math.sin(alpha) + u2Dot*Math.sin(alpha);
   			//%r2_dot_a2 = -alpha_dot*u2*sin(alpha) - (a+L)*alpha_dot*cos(alpha) + u2_dot*cos(alpha);
   			double r2_dot_a2 = alphaDot*u2*Math.sin(alpha) - (a+L)*alphaDot*Math.cos(alpha) + u2Dot*Math.cos(alpha);
   
   			//Hm1 (i,1)= r1a1*m*r1_dot_a2 - r1a2*m*r1_dot_a1;
   			//Hm2 (i,1)= r2a1*m*r2_dot_a2 - r2a2*m*r2_dot_a1;
     		double Hm1= r1a1*m*r1_dot_a2 - r1a2*m*r1_dot_a1;
   			double Hm2= r2a1*m*r2_dot_a2 - r2a2*m*r2_dot_a1;
   			
   			//angular_momentum(i,1) = Hsc(i,1)+Hm1(i,1)+Hm2(i,1);
			double angular_momentum = Hsc+Hm1+Hm2;

       	
        	/*
			private ThreePlots positionPlot = new ThreePlots();								
			private ThreePlots velocityPlot = new ThreePlots();
			private FourPlots hPlot = new FourPlots();
			private TwoPlots EandQPlot = new TwoPlots();
			private FourPlots quaternionPlot = new FourPlots();
			*/
        	        	       	
        	positionPlot.topPlot.addPoint(0, t, alpha, first);
        	positionPlot.middlePlot.addPoint(0, t, u1, first);
        	positionPlot.bottomPlot.addPoint(0, t, u2, first);
        	
        	velocityPlot.topPlot.addPoint(0, t, alphaDot, first);
        	velocityPlot.middlePlot.addPoint(0, t, u1Dot, first);
        	velocityPlot.bottomPlot.addPoint(0, t, u2Dot, first);
        	
        	hPlot.firstPlot.addPoint(0, t, Hsc, first);
        	hPlot.secondPlot.addPoint(0, t, Hm1, first);
        	hPlot.thirdPlot.addPoint(0, t, Hm2, first);
        	hPlot.fourthPlot.addPoint(0, t, angular_momentum, first);
        	
        	energyPlot.topPlot.addPoint(0, t, Potential, first);
        	energyPlot.middlePlot.addPoint(0, t,	Kinetic, first );
        	energyPlot.bottomPlot.addPoint(0, t, energy, first);
        	
        	quaternionPlot.firstPlot.addPoint(0, t, q1, first);
        	quaternionPlot.secondPlot.addPoint(0,t,q2, first);
        	quaternionPlot.thirdPlot.addPoint(0, t, q3, first);
        	quaternionPlot.fourthPlot.addPoint(0,t,q4, first);
        	
        	beamAnglePlot.topPlot.addPoint(0, t, theta1, first);
        	beamAnglePlot.bottomPlot.addPoint(0, t, theta2, first);
        	
        	quat_values[0][currentPts] = (float)t; // time value
        	quat_values[1][currentPts] = (float)q1; // quaternion 1
        	quat_values[2][currentPts] = (float)q2; // quarternion 2
        	quat_values[3][currentPts] = (float)q3; // quarternion 3
        	quat_values[4][currentPts] = (float)q4; // quarternion 4
        	
        	quatBeam1[0][currentPts] = (float)t; // time value
        	quatBeam1[1][currentPts] = (float)q1Beam1; // quaternion 1
        	quatBeam1[2][currentPts] = (float)q2Beam1; // quarternion 2
        	quatBeam1[3][currentPts] = (float)q3Beam1; // quarternion 3
        	quatBeam1[4][currentPts] = (float)q4Beam1; // quarternion 4
        	
        	quatBeam2[0][currentPts] = (float)t; // time value
        	quatBeam2[1][currentPts] = (float)q1Beam2; // quaternion 1
        	quatBeam2[2][currentPts] = (float)q2Beam2; // quarternion 2
        	quatBeam2[3][currentPts] = (float)q3Beam2; // quarternion 3
        	quatBeam2[4][currentPts] = (float)q4Beam2; // quarternion 4
        	//also print to the screen 
        	System.out.println(t+" "+y[0]+" "+y[1]+" "+y[2]+" "+y[3]+" "+y[4]+" "+y[5]);
    	}// End of print
    	
   
   /**
    * Return the quarternion values after simulation
    * @author	Noriko Takada
    */
   public float[][] getQuaternion()
   {
   		return quat_values;
   }
    
    /**
    * Return the quarternion values after simulation
    * @author	Noriko Takada
    */
   public float[][] getQuatBeam1()
   {
   		return quatBeam1;
   }
   
   /**
    * Return the quarternion values after simulation
    * @author	Noriko Takada
    */
   public float[][] getQuatBeam2()
   {
   		return quatBeam2;
   }
    	
    	
   /**
    * Make the plots visible after simulation
    * @author	Noriko Takada
    */
   public void makePlotsVisible()
   {
   		/*
		private ThreePlots positionPlot = new ThreePlots();								
		private ThreePlots velocityPlot = new ThreePlots();
		private FourPlots hPlot = new FourPlots();
		private TwoPlots EandQPlot = new TwoPlots();
		private FourPlots quaternionPlot = new FourPlots();
		*/
		
		positionPlot.setVisible(true);
		velocityPlot.setVisible(true);
		hPlot.setVisible(true);
		energyPlot.setVisible(true);
		quaternionPlot.setVisible(true);
		//beamAnglePlot.setVisible(true);
   }
   
   /** Runs the example.
    * @param args Arguments.
    */
    public static void main(String[] args)
    {

        double time_step=0.01;
        double timeDuration=10;
        double tf = 10;
        double t0 = 0.0;
               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
        RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78(1e-6);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
		float quatBeam1[][]= new  float[5][numberOfPts+1];
		float quatBeam2[][]= new  float[5][numberOfPts+1];
		 
        // create an instance
        FlexibleTwoD si = new FlexibleTwoD(time_step,  quat_values, quatBeam1, quatBeam2);
		
        // initialize the variables
        double [] x0 = new double[6];
        x0[0] = 0.0;
        x0[1] = 0.1;
        x0[2] = 0.5;
        x0[3] = 0.0;
        x0[4] = 0.0;
        x0[5] = 0.0;
        

        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        //double [] xf= rk78.integrate(t0, x0, tf, si, si, true);
        
        // make the plot visible
        si.makePlotsVisible();
        
    }
}
