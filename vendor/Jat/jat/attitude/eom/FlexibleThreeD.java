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
 * with two massless beams with a concentrated mass at the tip of
 * each beam, cantilevered at the central body on opposite sides. 
 *
 * @author Noriko Takada
 * @version (1.3): 8/15/2004
 * 	Modification since the last version
 * 		Switched to the interface EquationsOfMotion  
 */
public class FlexibleThreeD implements EquationsOfMotion
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
	double  time_step;
	private float quat_values[][];
	private float quatBeam1[][];
	private float quatBeam2[][];
	
	
	private double Ixx = 45.42;       // Spacecraft Principle inertia
    private double Iyy = 35.42;
    private double Izz = 13.35;
    private double m = 1;
    private double a = 0.40;
    private double L = 1.80;
    private double EI = 15.52;
    
    // Create variables for the necessary plots
	private FourPlots uPlot = new FourPlots();								
	private ThreePlots eulerPlot = new ThreePlots();
	private SinglePlot quaternion_check = new SinglePlot();
	private FourPlots quatPlot = new FourPlots();
	private ThreePlots energyPlot = new ThreePlots();
	private ThreePlots eulerRatePlot = new ThreePlots();
	private ThreePlots angularPlot = new ThreePlots();
	private FourPlots HPlot = new FourPlots();
	
	// For error detection
	private TwoPlots potentialPlot = new TwoPlots();
	private ThreePlots kineticPlot = new ThreePlots();
	private ThreePlots H1_BreakDown = new ThreePlots();
	private ThreePlots H2_BreakDown = new ThreePlots();
	private ThreePlots H3_BreakDown = new ThreePlots();
	
	/**
	 * Constructor
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 */
	public FlexibleThreeD(double time_step, float quat_values[][],float quatBeam1[][],
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
	public FlexibleThreeD(double time_step, double m, double a, double L, double EI,
											  double I1, double I2, double I3,float quat_values[][]
											  ,float quatBeam1[][], float quatBeam2[][])
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
		this.Ixx = I1;
		this.Iyy= I2;
		this.Izz = I3;
	}
	
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		
		/*
		private FourPlots uPlot = new FourPlots();								
		private ThreePlots eulerPlot = new ThreePlots();
		private SinglePlot quarternion_check = new SinglePlot();
		private FourPlots quatPlot = new FourPlots();
		private TwoPlots EandHPlot = new TwoPlots();
		private ThreePlots eulerRatePlot = new ThreePlots();
		private ThreePlots angularPlot = new ThreePlots();
		private ThreePlots HPlot = new ThreePlots();
		*/
		// Setup plots
		uPlot.setTitle("Tip mass position & Velocity");
        uPlot.firstPlot.setXLabel("t(sec)");
        uPlot.firstPlot.setYLabel("u1(m)");
        uPlot.secondPlot.setXLabel("t(sec)");
		uPlot.secondPlot.setYLabel("u2(m)");
        uPlot.thirdPlot.setXLabel("t(sec)");
        uPlot.thirdPlot.setYLabel("u1Dot(m/sec)");
        uPlot.fourthPlot.setXLabel("t(sec)");
        uPlot.fourthPlot.setYLabel("u2Dot(m/sec)");
        
        eulerPlot.setTitle("Euler angles");
        eulerPlot.topPlot.setXLabel("t (sec)");
        eulerPlot.topPlot.setYLabel("psi");
        eulerPlot.middlePlot.setXLabel("t (sec)");
        eulerPlot.middlePlot.setYLabel("theta");
        eulerPlot.bottomPlot.setXLabel("t (sec)");
        eulerPlot.bottomPlot.setYLabel("phi");
        
        quaternion_check.setTitle("Quarternion Check");
        quaternion_check.plot.setXLabel("t(sec)");
        quaternion_check.plot.setYLabel("q1^2 + q2^2 + q3^2 + q4^2");
        
        quatPlot.setTitle("Quaternion Values");
        quatPlot.firstPlot.setXLabel("t(sec)");
        quatPlot.firstPlot.setYLabel("q1");
        quatPlot.secondPlot.setXLabel("t (sec)");
        quatPlot.secondPlot.setYLabel("q2");
        quatPlot.thirdPlot.setXLabel("t (sec)");
        quatPlot.thirdPlot.setYLabel("q3");
        quatPlot.fourthPlot.setXLabel("t (sec)");
        quatPlot.fourthPlot.setYLabel("q4");
        
        energyPlot.setTitle("Energy ");
        energyPlot.topPlot.setXLabel("t(sec)");
        energyPlot.topPlot.setYLabel("Potential");
        energyPlot.middlePlot.setXLabel("t(sec)");
        energyPlot.middlePlot.setYLabel("Kinetic");
        energyPlot.bottomPlot.setXLabel("t (sec)");
        energyPlot.bottomPlot.setYLabel("Potential+Kinetic");
        
        eulerRatePlot.setTitle("Euler angle rates");
        eulerRatePlot.topPlot.setXLabel("t (sec)");
        eulerRatePlot.topPlot.setYLabel("psiDot");
        eulerRatePlot.middlePlot.setXLabel("t (sec)");
        eulerRatePlot.middlePlot.setYLabel("thetaDot");
        eulerRatePlot.bottomPlot.setXLabel("t (sec)");
        eulerRatePlot.bottomPlot.setYLabel("phiDot");
        
        
        angularPlot.setTitle("Angular Velocities");
        angularPlot.topPlot.setXLabel("t(sec)");
        angularPlot.topPlot.setYLabel("w1 (rad/sec)");
		angularPlot.middlePlot.setXLabel("t(sec)");
		angularPlot.middlePlot.setYLabel("w2 (rad/sec)");
        angularPlot.bottomPlot.setXLabel("t(sec)");
        angularPlot.bottomPlot.setYLabel("w3 (rad/sec)");
        
        HPlot.setTitle("Elements of Angular momentum vector");
        HPlot.firstPlot.setXLabel("t (sec)");
        HPlot.firstPlot.setYLabel("H1");
        HPlot.secondPlot.setXLabel("t (sec)");
        HPlot.secondPlot.setYLabel("H2");
        HPlot.thirdPlot.setXLabel("t (sec)");
        HPlot.thirdPlot.setYLabel("H3");
        HPlot.fourthPlot.setXLabel("t (sec)");
        HPlot.fourthPlot.setYLabel("Total Angular Momentum");
               
        potentialPlot.setTitle("Potential energy break-down");
        potentialPlot.topPlot.setXLabel("t (sec)");
        potentialPlot.topPlot.setYLabel("Potential 1");
        potentialPlot.bottomPlot.setXLabel("t (sec)");
        potentialPlot.bottomPlot.setYLabel("Potentil 2");
        
        kineticPlot.setTitle("Kinetic energy break-down");
        kineticPlot.topPlot.setXLabel("t (sec)");
        kineticPlot.topPlot.setYLabel("Kinetic Spacecraft");
        kineticPlot.middlePlot.setXLabel("t (sec)");
        kineticPlot.middlePlot.setYLabel("Kinetic Particle 1");
        kineticPlot.bottomPlot.setXLabel("t (sec)");
        kineticPlot.bottomPlot.setYLabel("Kinetic Particle 2"); 
        
        H1_BreakDown.setTitle("H1 Break Down");
        H1_BreakDown.topPlot.setXLabel("t (sec)");
        H1_BreakDown.topPlot.setYLabel("H1_sc");
        H1_BreakDown.middlePlot.setXLabel("t (sec)");
        H1_BreakDown.middlePlot.setYLabel("H1_m1");
        H1_BreakDown.bottomPlot.setXLabel("t (sec)");
        H1_BreakDown.bottomPlot.setYLabel("H1_m2");
        
        H2_BreakDown.setTitle("H2 Break Down");
        H2_BreakDown.topPlot.setXLabel("t (sec)");
        H2_BreakDown.topPlot.setYLabel("H2_sc");
        H2_BreakDown.middlePlot.setXLabel("t (sec)");
        H2_BreakDown.middlePlot.setYLabel("H2_m1");
        H2_BreakDown.bottomPlot.setXLabel("t (sec)");
        H2_BreakDown.bottomPlot.setYLabel("H2_m2");
        
        H3_BreakDown.setTitle("H3 Break Down");
        H3_BreakDown.topPlot.setXLabel("t (sec)");
        H3_BreakDown.topPlot.setYLabel("H3_sc");
        H3_BreakDown.middlePlot.setXLabel("t (sec)");
        H3_BreakDown.middlePlot.setYLabel("H3_m1");
        H3_BreakDown.bottomPlot.setXLabel("t (sec)");
        H3_BreakDown.bottomPlot.setYLabel("H3_m2");       
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
       	
       	double [] out = new double[10];
       	// Define state variables
       	double	u1 			= x[0];
       	double u2 			= x[1];
       	double psi 		= x[2];
       	double theta 		= x[3];
       	double phi 		= x[4];
       	double u1Dot 		= x[5];
       	double u2Dot 		= x[6];
       	double psiDot 		= x[7];
       	double thetaDot 	= x[8];
       	double phiDot 		= x[9];
       	
       	// These are an exact copy from MATLAB!
       	double Lsquare = L*L;
		double Lcube = L*L*L;
	 	double Lquad = L*L*L*L;
       	double u1square = u1*u1;
       	double u2square = u2*u2;
       	double  asquare = a*a;
       	double  psiDotSquare = psiDot*psiDot;
       	double  thetaDotSquare = thetaDot*thetaDot;
       	double  phiDotSquare = phiDot*phiDot;
       	double  cos_psi_square = Math.cos(psi)*Math.cos(psi);
       	double  sin_psi_square = Math.sin(psi)*Math.sin(psi);
       	double cos_theta_square = Math.cos(theta)*Math.cos(theta);
       	double  sin_theta_square = Math.sin(theta)*Math.sin(theta);
       	double cos_phi_square = Math.cos(phi)*Math.cos(phi);
       	double  sin_phi_square = Math.sin(phi)*Math.sin(phi);
       	double  c_psi = Math.cos(psi);
       	double  s_psi = Math.sin(psi);
       	double  c_theta = Math.cos(theta);
       	double s_theta = Math.sin(theta);
       	double  c_phi = Math.cos(phi);
       	double  s_phi = Math.sin(phi);
       	
      	double  	 A11 = m;
       	double  A12 = 0;
       	double  A13 = (m*a+m*L)*c_theta*c_phi;
       	double  A14 = (-m*a-m*L)*s_phi;
       	double  A15 = 0;
       	double  B1 = (-3*EI*u1+m*u1*psiDotSquare*Lcube*cos_theta_square*cos_phi_square-m*u1*thetaDotSquare*Lcube*cos_phi_square-m*u1*psiDotSquare*Lcube*cos_theta_square+m*u1*Lcube*phiDotSquare+m*u1*Lcube*psiDotSquare+m*u1*thetaDotSquare*Lcube-2*m*phiDot*psiDot*u1*Lcube*s_theta+(m*a*Lcube*psiDotSquare+m*Lquad*psiDotSquare)*c_theta*s_theta*s_phi+(2*m*psiDot*a*thetaDot*Lcube+2*m*psiDot*thetaDot*Lquad)*c_phi*s_theta-2*m*psiDot*u1*thetaDot*Lcube*c_theta*c_phi*s_phi)/Lcube;
       	double A21 = 0;
       	double A22 = m;
       	double A23 = (-m*a-m*L)*c_theta*c_phi;
       	double A24 = (m*a+m*L)*s_phi;
       	double A25 = 0;
       	double  B2 = (-3*EI*u2+m*u2*psiDotSquare*Lcube*cos_theta_square*cos_phi_square-m*u2*thetaDotSquare*Lcube*cos_phi_square-m*u2*psiDotSquare*Lcube*cos_theta_square+m*u2*Lcube*phiDotSquare+m*u2*Lcube*psiDotSquare+m*u2*thetaDotSquare*Lcube+(-2*m*psiDot*a*thetaDot*Lcube-2*m*psiDot*thetaDot*Lquad)*c_phi*s_theta+(-m*a*Lcube*psiDotSquare-m*Lquad*psiDotSquare)*c_theta*s_theta*s_phi-2*m*psiDot*u2*thetaDot*Lcube*c_theta*c_phi*s_phi-2*m*phiDot*psiDot*u2*Lcube*s_theta)/Lcube;
       	double A31 = (m*a+m*L)*c_theta*c_phi;
       	double  A32 = (-m*a-m*L)*c_theta*c_phi;
       	double  A33 = -Ixx*cos_theta_square+Iyy*cos_theta_square+m*u1square+m*u2square+4*m*L*a*cos_theta_square+2*m*asquare*cos_theta_square+2*m*Lsquare*cos_theta_square+m*u1square*cos_theta_square*cos_phi_square-m*u1square*cos_theta_square+m*u2square*cos_theta_square*cos_phi_square-m*u2square*cos_theta_square-Iyy*cos_theta_square*cos_phi_square+Izz*cos_theta_square*cos_phi_square+Ixx+(-2*m*u2*a+2*m*u1*L-2*m*L*u2+2*m*u1*a)*c_theta*s_theta*s_phi;
       	double A34 = (-m*u2*a+m*u1*L-m*L*u2+m*u1*a)*c_phi*s_theta+(Iyy-Izz-m*u1square-m*u2square)*c_theta*c_phi*s_phi;
       	double A35 = (-Ixx-m*u1square-m*u2square)*s_theta+(m*u2*a+m*L*u2-m*u1*L-m*u1*a)*c_theta*s_phi;
       	double B3 = -2*m*u1Dot*psiDot*u1-2*m*u2Dot*psiDot*u2-2*m*u1Dot*psiDot*u1*cos_theta_square*cos_phi_square-2*m*u2Dot*psiDot*u2*cos_theta_square*cos_phi_square+2*m*u1Dot*psiDot*u1*cos_theta_square+2*m*u2Dot*psiDot*u2*cos_theta_square+(2*m*phiDot*psiDot*u2square*cos_theta_square+2*phiDot*psiDot*Izz*cos_theta_square-2*phiDot*psiDot*Iyy*cos_theta_square+2*m*phiDot*psiDot*u1square*cos_theta_square)*c_phi*s_phi+(Iyy*thetaDotSquare-m*thetaDotSquare*u1square-m*thetaDotSquare*u2square-Izz*thetaDotSquare)*c_phi*s_theta*s_phi+(Ixx*thetaDot*phiDot-Izz*phiDot*thetaDot+Iyy*thetaDot*phiDot+2*m*thetaDot*phiDot*u1square*cos_phi_square+2*m*thetaDot*phiDot*u2square*cos_phi_square+2*Izz*phiDot*thetaDot*cos_phi_square-2*Iyy*thetaDot*phiDot*cos_phi_square)*c_theta+(2*m*u1Dot*phiDot*u1+2*m*u2Dot*phiDot*u2)*s_theta+(-4*m*a*psiDot*thetaDot*u1*cos_theta_square-2*m*L*psiDot*thetaDot*u2-4*m*L*psiDot*thetaDot*u1*cos_theta_square-2*m*a*psiDot*thetaDot*u2+4*m*L*psiDot*thetaDot*u2*cos_theta_square+2*m*a*psiDot*thetaDot*u1+4*m*a*psiDot*thetaDot*u2*cos_theta_square+2*m*L*psiDot*thetaDot*u1)*s_phi+(-2*Iyy*psiDot*thetaDot*cos_phi_square+4*m*thetaDot*psiDot*asquare-2*Ixx*psiDot*thetaDot-2*m*thetaDot*psiDot*u2square+2*m*thetaDot*psiDot*u2square*cos_phi_square+2*Izz*psiDot*thetaDot*cos_phi_square+2*Iyy*psiDot*thetaDot-2*m*thetaDot*psiDot*u1square+2*m*thetaDot*psiDot*u1square*cos_phi_square+4*m*thetaDot*psiDot*Lsquare+8*m*L*thetaDot*psiDot*a)*c_theta*s_theta+(m*u1*L*phiDotSquare-m*u1*L*thetaDotSquare+m*u2*L*thetaDotSquare-m*u1*a*thetaDotSquare-m*u2*a*phiDotSquare-m*u2*L*phiDotSquare+m*u1*a*phiDotSquare+m*u2*a*thetaDotSquare)*c_theta*c_phi+(-2*m*L*psiDot*phiDot*u1+2*m*L*psiDot*phiDot*u2+2*m*a*psiDot*phiDot*u2-2*m*a*psiDot*phiDot*u1)*c_theta*c_phi*s_theta+(2*m*u1*u1Dot*thetaDot+2*m*u2*u2Dot*thetaDot)*c_theta*c_phi*s_phi+(2*m*u2Dot*psiDot*L+2*m*u2Dot*psiDot*a-2*m*u1Dot*psiDot*a-2*m*u1Dot*psiDot*L)*c_theta*s_theta*s_phi+(-2*m*phiDot*u2Dot*L+2*m*phiDot*u1Dot*L-2*m*phiDot*u2Dot*a+2*m*phiDot*u1Dot*a)*c_theta*s_phi;
       	double A41 = (-m*a-m*L)*s_phi;
       	double A42 = (m*a+m*L)*s_phi;
       	double A43 = (-m*u2*a+m*u1*L-m*L*u2+m*u1*a)*c_phi*s_theta+(Iyy-Izz-m*u1square-m*u2square)*c_theta*c_phi*s_phi;
       	double A44 = 4*m*L*a+2*m*asquare+2*m*Lsquare-m*u1square*cos_phi_square+m*u1square-m*u2square*cos_phi_square+m*u2square+Iyy*cos_phi_square-Izz*cos_phi_square+Izz;
       	double A45 = (m*u2*a+m*L*u2-m*u1*L-m*u1*a)*c_phi;
       	double B4 = -2*m*u2*u2Dot*thetaDot-2*m*u1*u1Dot*thetaDot+2*m*u2*u2Dot*thetaDot*cos_phi_square+2*m*u1*u1Dot*thetaDot*cos_phi_square+(2*m*a*psiDot*phiDot*u1-2*m*a*psiDot*phiDot*u2+2*m*L*psiDot*phiDot*u1-2*m*L*psiDot*phiDot*u2)*s_phi*s_theta+(-2*thetaDot*phiDot*m*u2square-2*thetaDot*phiDot*m*u1square+2*Iyy*thetaDot*phiDot-2*Izz*phiDot*thetaDot)*s_phi*c_phi+(psiDotSquare*Iyy*cos_phi_square-psiDotSquare*m*u1square*cos_phi_square-2*psiDotSquare*m*Lsquare-psiDotSquare*Iyy-psiDotSquare*m*u2square*cos_phi_square-psiDotSquare*Izz*cos_phi_square+psiDotSquare*m*u1square+psiDotSquare*Ixx+psiDotSquare*m*u2square-4*psiDotSquare*m*L*a-2*psiDotSquare*m*asquare)*s_theta*c_theta+(2*m*u2Dot*psiDot*L+2*m*u2Dot*psiDot*a-2*m*u1Dot*psiDot*a-2*m*u1Dot*psiDot*L)*c_phi*s_theta+(-Ixx*psiDot*phiDot-Izz*psiDot*phiDot+Iyy*phiDot*psiDot+2*m*phiDot*psiDot*u1square*cos_phi_square+2*Izz*psiDot*phiDot*cos_phi_square-2*Iyy*phiDot*psiDot*cos_phi_square+2*m*phiDot*psiDot*u2square*cos_phi_square-2*m*phiDot*psiDot*u1square-2*m*phiDot*psiDot*u2square)*c_theta+(-2*m*phiDot*u2Dot*L+2*m*phiDot*u1Dot*L-2*m*phiDot*u2Dot*a+2*m*phiDot*u1Dot*a)*c_phi+(2*m*L*u1*psiDotSquare*cos_theta_square-m*u1*L*psiDotSquare-m*u1*L*phiDotSquare+m*u2*a*psiDotSquare-2*m*u2*a*psiDotSquare*cos_theta_square-m*u1*a*psiDotSquare+2*m*u1*a*psiDotSquare*cos_theta_square-2*m*L*u2*psiDotSquare*cos_theta_square+m*u2*L*psiDotSquare+m*u2*a*phiDotSquare-m*u1*a*phiDotSquare+m*u2*L*phiDotSquare)*s_phi+(2*m*u1Dot*psiDot*u1+2*m*u2Dot*psiDot*u2)*c_theta*c_phi*s_phi;
       	double A51 = 0;
       	double A52 = 0;
       	double A53 = (-Ixx-m*u1square-m*u2square)*s_theta+(m*u2*a+m*L*u2-m*u1*L-m*u1*a)*c_theta*s_phi;
       	double A54 = (m*u2*a+m*L*u2-m*u1*L-m*u1*a)*c_phi;
       	double A55 = Ixx+m*u1square+m*u2square;
       	double B5 = -2*m*u1Dot*phiDot*u1-2*m*u2Dot*phiDot*u2+(-2*m*a*psiDot*thetaDot*u1+2*m*a*psiDot*thetaDot*u2+2*m*L*psiDot*thetaDot*u2-2*m*L*psiDot*thetaDot*u1)*s_phi*s_theta+(-m*psiDotSquare*u1square*cos_theta_square+m*thetaDotSquare*u2square+m*thetaDotSquare*u1square+psiDotSquare*Iyy*cos_theta_square-psiDotSquare*Izz*cos_theta_square-Iyy*thetaDotSquare+Izz*thetaDotSquare-m*psiDotSquare*u2square*cos_theta_square)*s_phi*c_phi+(m*u1*L*psiDotSquare-m*u2*L*psiDotSquare+m*u1*a*psiDotSquare-m*u2*a*psiDotSquare)*c_theta*c_phi*s_theta+(Ixx*psiDot*thetaDot+Izz*psiDot*thetaDot-Iyy*psiDot*thetaDot-2*m*thetaDot*psiDot*u1square*cos_phi_square-2*m*thetaDot*psiDot*u2square*cos_phi_square+2*Iyy*psiDot*thetaDot*cos_phi_square-2*Izz*psiDot*thetaDot*cos_phi_square+2*m*thetaDot*psiDot*u1square+2*m*thetaDot*psiDot*u2square)*c_theta+(2*m*u1Dot*psiDot*u1+2*m*u2Dot*psiDot*u2)*s_theta;
       	





  	out[0]=x[5];
  	out[1]=x[6];
  	out[2]=x[7];
  	out[3]=x[8];
  	out[4]=x[9];
  
  	//%--------u1DoubleDot Equation---------------%
  	double u1Den= A23*A41*A32*A14*A55-A23*A51*A34*A42*A15-A11*A42*A34*A53*A25+A31*A12*A23*A44*A55-A31*A12*A23*A54*A45+A31*A12*A54*A43*A25-A31*A12*A43*A24*A55+A31*A12*A53*A24*A45-A31*A12*A44*A53*A25-A31*A23*A44*A52*A15-A31*A23*A42*A14*A55+A31*A23*A54*A42*A15+A31*A23*A52*A14*A45-A31*A54*A43*A22*A15+A31*A54*A22*A13*A45-A31*A52*A13*A24*A45-A31*A52*A14*A43*A25-A31*A44*A22*A13*A55+A31*A44*A53*A22*A15-A31*A54*A42*A13*A25+A31*A43*A22*A14*A55-A31*A53*A22*A14*A45+A31*A42*A14*A53*A25-A31*A53*A24*A42*A15-A23*A54*A11*A42*A35+A23*A42*A14*A51*A35+A23*A44*A11*A52*A35+A23*A54*A11*A32*A45-A23*A54*A41*A32*A15-A23*A11*A52*A34*A45+A23*A41*A34*A52*A15+A23*A44*A51*A32*A15-A23*A51*A32*A14*A45-A23*A52*A14*A41*A35-A23*A44*A11*A32*A55-A12*A23*A41*A34*A55+A12*A23*A51*A34*A45-A12*A23*A44*A51*A35+A12*A23*A54*A41*A35+A12*A43*A24*A51*A35-A12*A21*A53*A34*A45-A12*A21*A33*A44*A55+A12*A33*A41*A24*A55+A12*A21*A43*A34*A55-A12*A53*A24*A41*A35+A12*A44*A21*A53*A35+A12*A44*A33*A51*A25+A12*A21*A33*A54*A45-A12*A54*A33*A41*A25+A12*A41*A34*A53*A25-A12*A33*A51*A24*A45-A12*A54*A21*A43*A35-A12*A51*A34*A43*A25-A21*A42*A13*A34*A55-A11*A22*A43*A34*A55+A21*A53*A32*A14*A45+A52*A14*A33*A41*A25-A33*A22*A41*A14*A55+A52*A13*A24*A41*A35+A52*A14*A21*A43*A35+A53*A24*A41*A32*A15+A51*A34*A43*A22*A15+A53*A24*A11*A42*A35+A51*A34*A42*A13*A25-A21*A33*A52*A14*A45-A32*A13*A41*A24*A55-A33*A11*A42*A24*A55-A51*A34*A22*A13*A45-A44*A32*A13*A51*A25+A32*A13*A51*A24*A45+A33*A11*A52*A24*A45-A21*A33*A54*A42*A15-A21*A43*A34*A52*A15+A44*A21*A32*A13*A55-A44*A21*A52*A13*A35+A44*A11*A22*A33*A55-A44*A11*A22*A53*A35+A44*A22*A51*A13*A35+A44*A53*A11*A32*A25-A41*A34*A52*A13*A25+A21*A33*A42*A14*A55+A11*A52*A34*A43*A25-A44*A21*A53*A32*A15-A44*A33*A11*A52*A25-A44*A33*A22*A51*A15+A21*A52*A13*A34*A45+A11*A22*A53*A34*A45-A33*A41*A24*A52*A15-A54*A11*A22*A33*A45+A54*A32*A13*A41*A25+A54*A33*A22*A41*A15+A54*A11*A22*A43*A35-A42*A13*A24*A51*A35-A41*A32*A14*A53*A25+A21*A33*A44*A52*A15-A42*A14*A21*A53*A35-A42*A14*A33*A51*A25+A41*A34*A22*A13*A55-A41*A34*A53*A22*A15-A54*A21*A32*A13*A45-A43*A22*A14*A51*A35+A33*A22*A51*A14*A45+A53*A22*A14*A41*A35+A51*A32*A14*A43*A25-A21*A43*A32*A14*A55-A54*A22*A41*A13*A35-A54*A43*A11*A32*A25+A54*A33*A11*A42*A25+A33*A51*A24*A42*A15-A53*A24*A11*A32*A45+A21*A53*A34*A42*A15+A54*A21*A42*A13*A35+A54*A21*A43*A32*A15+A43*A24*A11*A32*A55-A43*A24*A51*A32*A15-A43*A24*A11*A52*A35+A23*A11*A42*A34*A55+A31*A44*A52*A13*A25+A31*A43*A24*A52*A15+A31*A42*A13*A24*A55;
  	double u1Num= -B5*A12*A33*A24*A45+B1*A24*A33*A45*A52+B5*A12*A44*A33*A25-B1*A43*A55*A34*A22-B1*A43*A54*A25*A32+B1*A43*A35*A54*A22+B1*A43*A34*A25*A52+B5*A12*A43*A24*A35-B1*A43*A24*A35*A52+B1*A43*A24*A55*A32-B1*A35*A23*A54*A42+B1*A53*A45*A34*A22-B1*A35*A44*A53*A22-B1*A34*A25*A53*A42-B1*A44*A33*A25*A52-B5*A12*A34*A43*A25-B1*A23*A34*A45*A52-B1*A33*A45*A54*A22+B1*A54*A25*A33*A42+B1*A23*A44*A35*A52-B1*A24*A55*A33*A42-B5*A12*A23*A44*A35+B4*A54*A22*A15*A33+B2*A12*A43*A55*A34-B1*A24*A53*A45*A32+B1*A24*A35*A53*A42+B5*A23*A44*A32*A15+B4*A12*A24*A33*A55-B2*A12*A43*A35*A54+B5*A12*A23*A34*A45+B4*A55*A13*A22*A34+B1*A55*A23*A34*A42-B1*A55*A23*A44*A32-B4*A34*A25*A52*A13+B1*A23*A54*A45*A32+B1*A44*A53*A25*A32+B1*A55*A44*A22*A33+B2*A52*A13*A45*A34-B2*A32*A15*A44*A53-B2*A55*A42*A13*A34+B2*A43*A35*A52*A14-B4*A35*A22*A13*A54-B2*A43*A52*A15*A34+B4*A55*A23*A32*A14+B2*A12*A35*A44*A53-B2*A43*A55*A32*A14+B2*A43*A32*A15*A54+B2*A12*A33*A45*A54-B2*A12*A53*A45*A34-B2*A12*A55*A33*A44-B4*A34*A22*A15*A53-B2*A35*A42*A14*A53+B2*A55*A32*A13*A44-B2*A33*A42*A15*A54+B2*A53*A42*A15*A34+B2*A52*A15*A33*A44+B5*A23*A42*A14*A35-B4*A55*A33*A14*A22-B2*A35*A44*A52*A13+B2*A53*A45*A32*A14-B2*A33*A45*A52*A14+B4*A52*A15*A23*A34+B4*A54*A25*A32*A13+B2*A35*A54*A42*A13+B2*A55*A33*A42*A14+B4*A24*A53*A32*A15-B4*A32*A15*A23*A54+B4*A24*A52*A13*A35-B4*A24*A32*A13*A55+B4*A52*A14*A25*A33+B4*A12*A35*A23*A54-B4*A24*A33*A52*A15-B4*A32*A14*A25*A53+B4*A12*A34*A25*A53+B4*A35*A53*A22*A14-B4*A12*A55*A23*A34-B4*A35*A23*A52*A14-B4*A12*A54*A25*A33-B5*A44*A33*A15*A22-B4*A12*A24*A53*A35-B5*A34*A22*A13*A45-B3*A53*A24*A42*A15-B5*A43*A24*A32*A15+B3*A42*A14*A53*A25-B5*A43*A22*A14*A35-B3*A53*A22*A14*A45-B5*A42*A14*A33*A25+B3*A43*A22*A14*A55+B5*A33*A24*A42*A15-B3*A54*A42*A13*A25+B5*A34*A42*A13*A25-B5*A42*A13*A24*A35+B3*A42*A13*A24*A55+B5*A34*A43*A22*A15+B3*A43*A24*A52*A15+B3*A44*A52*A13*A25+B5*A44*A13*A22*A35+B3*A44*A53*A22*A15-B3*A44*A22*A13*A55-B5*A44*A32*A13*A25-B3*A52*A14*A43*A25-B3*A52*A13*A24*A45+B3*A54*A22*A13*A45-B3*A54*A43*A22*A15+B3*A23*A54*A42*A15+B3*A23*A52*A14*A45-B3*A23*A42*A14*A55+B3*A12*A53*A24*A45-B3*A12*A44*A53*A25-B3*A23*A44*A52*A15-B3*A12*A43*A24*A55+B3*A12*A54*A43*A25-B3*A12*A23*A54*A45-B2*A32*A13*A45*A54+B5*A32*A13*A24*A45+B3*A12*A23*A44*A55+B5*A32*A14*A43*A25+B5*A33*A14*A22*A45-B5*A23*A34*A42*A15-B5*A23*A32*A14*A45;
  	out[5]=  u1Num/u1Den;
  
  	//%--------u2DoubleDot Equation---------------%
  	double u2Den= A23*A41*A32*A14*A55-A23*A51*A34*A42*A15-A11*A42*A34*A53*A25+A31*A12*A23*A44*A55-A31*A12*A23*A54*A45+A31*A12*A54*A43*A25-A31*A12*A43*A24*A55+A31*A12*A53*A24*A45-A31*A12*A44*A53*A25-A31*A23*A44*A52*A15-A31*A23*A42*A14*A55+A31*A23*A54*A42*A15+A31*A23*A52*A14*A45-A31*A54*A43*A22*A15+A31*A54*A22*A13*A45-A31*A52*A13*A24*A45-A31*A52*A14*A43*A25-A31*A44*A22*A13*A55+A31*A44*A53*A22*A15-A31*A54*A42*A13*A25+A31*A43*A22*A14*A55-A31*A53*A22*A14*A45+A31*A42*A14*A53*A25-A31*A53*A24*A42*A15-A23*A54*A11*A42*A35+A23*A42*A14*A51*A35+A23*A44*A11*A52*A35+A23*A54*A11*A32*A45-A23*A54*A41*A32*A15-A23*A11*A52*A34*A45+A23*A41*A34*A52*A15+A23*A44*A51*A32*A15-A23*A51*A32*A14*A45-A23*A52*A14*A41*A35-A23*A44*A11*A32*A55-A12*A23*A41*A34*A55+A12*A23*A51*A34*A45-A12*A23*A44*A51*A35+A12*A23*A54*A41*A35+A12*A43*A24*A51*A35-A12*A21*A53*A34*A45-A12*A21*A33*A44*A55+A12*A33*A41*A24*A55+A12*A21*A43*A34*A55-A12*A53*A24*A41*A35+A12*A44*A21*A53*A35+A12*A44*A33*A51*A25+A12*A21*A33*A54*A45-A12*A54*A33*A41*A25+A12*A41*A34*A53*A25-A12*A33*A51*A24*A45-A12*A54*A21*A43*A35-A12*A51*A34*A43*A25-A21*A42*A13*A34*A55-A11*A22*A43*A34*A55+A21*A53*A32*A14*A45+A52*A14*A33*A41*A25-A33*A22*A41*A14*A55+A52*A13*A24*A41*A35+A52*A14*A21*A43*A35+A53*A24*A41*A32*A15+A51*A34*A43*A22*A15+A53*A24*A11*A42*A35+A51*A34*A42*A13*A25-A21*A33*A52*A14*A45-A32*A13*A41*A24*A55-A33*A11*A42*A24*A55-A51*A34*A22*A13*A45-A44*A32*A13*A51*A25+A32*A13*A51*A24*A45+A33*A11*A52*A24*A45-A21*A33*A54*A42*A15-A21*A43*A34*A52*A15+A44*A21*A32*A13*A55-A44*A21*A52*A13*A35+A44*A11*A22*A33*A55-A44*A11*A22*A53*A35+A44*A22*A51*A13*A35+A44*A53*A11*A32*A25-A41*A34*A52*A13*A25+A21*A33*A42*A14*A55+A11*A52*A34*A43*A25-A44*A21*A53*A32*A15-A44*A33*A11*A52*A25-A44*A33*A22*A51*A15+A21*A52*A13*A34*A45+A11*A22*A53*A34*A45-A33*A41*A24*A52*A15-A54*A11*A22*A33*A45+A54*A32*A13*A41*A25+A54*A33*A22*A41*A15+A54*A11*A22*A43*A35-A42*A13*A24*A51*A35-A41*A32*A14*A53*A25+A21*A33*A44*A52*A15-A42*A14*A21*A53*A35-A42*A14*A33*A51*A25+A41*A34*A22*A13*A55-A41*A34*A53*A22*A15-A54*A21*A32*A13*A45-A43*A22*A14*A51*A35+A33*A22*A51*A14*A45+A53*A22*A14*A41*A35+A51*A32*A14*A43*A25-A21*A43*A32*A14*A55-A54*A22*A41*A13*A35-A54*A43*A11*A32*A25+A54*A33*A11*A42*A25+A33*A51*A24*A42*A15-A53*A24*A11*A32*A45+A21*A53*A34*A42*A15+A54*A21*A42*A13*A35+A54*A21*A43*A32*A15+A43*A24*A11*A32*A55-A43*A24*A51*A32*A15-A43*A24*A11*A52*A35+A23*A11*A42*A34*A55+A31*A44*A52*A13*A25+A31*A43*A24*A52*A15+A31*A42*A13*A24*A55;
  	double u2Num= -B3*A24*A55*A13*A41+B3*A24*A13*A45*A51-B3*A24*A53*A45*A11-B3*A44*A13*A25*A51+B3*A24*A15*A53*A41-B3*A43*A55*A14*A21+B3*A43*A14*A25*A51+B3*A43*A15*A54*A21-B2*A33*A11*A45*A54+B1*A31*A53*A45*A24-B5*A23*A11*A34*A45+B5*A23*A44*A11*A35-B3*A43*A54*A25*A11+B5*A23*A41*A34*A15-B2*A31*A55*A44*A13+B4*A31*A14*A25*A53-B5*A23*A14*A41*A35-B3*A43*A24*A15*A51-B5*A21*A33*A14*A45+B3*A43*A24*A55*A11-B2*A11*A35*A44*A53+B2*A15*A54*A33*A41+B1*A31*A55*A23*A44+B2*A53*A45*A11*A34-B5*A41*A34*A13*A25+B5*A21*A13*A34*A45+B5*A21*A33*A44*A15+B1*A31*A43*A54*A25-B1*A31*A43*A55*A24+B2*A31*A15*A44*A53-B2*A13*A41*A35*A54-B1*A53*A45*A21*A34-B2*A15*A41*A34*A53+B2*A51*A35*A44*A13+B1*A41*A34*A53*A25+B1*A21*A33*A54*A45-B5*A44*A33*A11*A25+B2*A55*A44*A11*A33-B5*A44*A21*A13*A35+B1*A23*A54*A41*A35+B3*A23*A54*A45*A11-B4*A31*A54*A25*A13+B2*A55*A41*A34*A13+B1*A44*A21*A53*A35+B2*A33*A51*A45*A14-B1*A55*A41*A23*A34-B2*A55*A14*A33*A41-B1*A54*A25*A33*A41-B2*A13*A45*A51*A34+B1*A44*A33*A51*A25-B2*A15*A44*A33*A51+B1*A23*A51*A34*A45+B1*A24*A55*A33*A41+B2*A43*A15*A51*A34-B1*A51*A35*A23*A44+B2*A53*A41*A35*A14-B2*A43*A51*A35*A14+B5*A33*A11*A24*A45-B1*A55*A21*A33*A44+B5*A13*A24*A41*A35-B5*A33*A41*A24*A15+B5*A14*A33*A41*A25+B1*A43*A55*A21*A34-B1*A24*A53*A41*A35-B5*A31*A23*A44*A15-B1*A24*A33*A51*A45+B5*A31*A44*A13*A25-B1*A43*A54*A21*A35-B1*A43*A51*A34*A25-B5*A43*A24*A11*A35-B1*A31*A23*A54*A45+B5*A14*A21*A43*A35-B1*A31*A44*A53*A25+B5*A11*A34*A43*A25-B5*A21*A43*A34*A15+B1*A43*A51*A35*A24+B2*A31*A13*A45*A54-B5*A31*A14*A43*A25-B5*A31*A13*A24*A45+B2*A43*A11*A35*A54-B2*A31*A53*A45*A14-B4*A15*A51*A23*A34-B2*A31*A43*A15*A54+B5*A31*A23*A14*A45-B2*A43*A55*A11*A34+B4*A54*A21*A35*A13+B2*A31*A43*A55*A14+B3*A54*A25*A13*A41+B5*A31*A43*A24*A15-B3*A15*A23*A54*A41-B3*A23*A14*A45*A51+B3*A53*A45*A14*A21+B4*A51*A35*A23*A14+B4*A55*A11*A23*A34-B4*A31*A24*A15*A53-B4*A11*A35*A23*A54+B4*A15*A21*A53*A34+B4*A31*A24*A55*A13+B4*A54*A25*A33*A11-B4*A14*A25*A33*A51-B4*A11*A34*A25*A53-B3*A14*A25*A53*A41-B4*A14*A21*A35*A53-B4*A55*A21*A13*A34+B4*A24*A11*A35*A53-B3*A13*A45*A54*A21+B3*A44*A53*A25*A11-B3*A15*A44*A21*A53-B4*A24*A51*A35*A13-B4*A24*A55*A33*A11-B4*A15*A21*A33*A54+B3*A55*A44*A21*A13+B4*A51*A34*A25*A13-B3*A55*A23*A11*A44+B4*A55*A21*A33*A14+B4*A31*A15*A23*A54+B3*A55*A23*A41*A14-B4*A31*A55*A23*A14+B3*A23*A44*A15*A51+B4*A24*A15*A33*A51;
    
  	out[6]=u2Num/u2Den;
  	//%---------- psiDoubleDot Equation ----------%
	double psiDen= A23*A41*A32*A14*A55-A23*A51*A34*A42*A15-A11*A42*A34*A53*A25+A31*A12*A23*A44*A55-A31*A12*A23*A54*A45+A31*A12*A54*A43*A25-A31*A12*A43*A24*A55+A31*A12*A53*A24*A45-A31*A12*A44*A53*A25-A31*A23*A44*A52*A15-A31*A23*A42*A14*A55+A31*A23*A54*A42*A15+A31*A23*A52*A14*A45-A31*A54*A43*A22*A15+A31*A54*A22*A13*A45-A31*A52*A13*A24*A45-A31*A52*A14*A43*A25-A31*A44*A22*A13*A55+A31*A44*A53*A22*A15-A31*A54*A42*A13*A25+A31*A43*A22*A14*A55-A31*A53*A22*A14*A45+A31*A42*A14*A53*A25-A31*A53*A24*A42*A15-A23*A54*A11*A42*A35+A23*A42*A14*A51*A35+A23*A44*A11*A52*A35+A23*A54*A11*A32*A45-A23*A54*A41*A32*A15-A23*A11*A52*A34*A45+A23*A41*A34*A52*A15+A23*A44*A51*A32*A15-A23*A51*A32*A14*A45-A23*A52*A14*A41*A35-A23*A44*A11*A32*A55-A12*A23*A41*A34*A55+A12*A23*A51*A34*A45-A12*A23*A44*A51*A35+A12*A23*A54*A41*A35+A12*A43*A24*A51*A35-A12*A21*A53*A34*A45-A12*A21*A33*A44*A55+A12*A33*A41*A24*A55+A12*A21*A43*A34*A55-A12*A53*A24*A41*A35+A12*A44*A21*A53*A35+A12*A44*A33*A51*A25+A12*A21*A33*A54*A45-A12*A54*A33*A41*A25+A12*A41*A34*A53*A25-A12*A33*A51*A24*A45-A12*A54*A21*A43*A35-A12*A51*A34*A43*A25-A21*A42*A13*A34*A55-A11*A22*A43*A34*A55+A21*A53*A32*A14*A45+A52*A14*A33*A41*A25-A33*A22*A41*A14*A55+A52*A13*A24*A41*A35+A52*A14*A21*A43*A35+A53*A24*A41*A32*A15+A51*A34*A43*A22*A15+A53*A24*A11*A42*A35+A51*A34*A42*A13*A25-A21*A33*A52*A14*A45-A32*A13*A41*A24*A55-A33*A11*A42*A24*A55-A51*A34*A22*A13*A45-A44*A32*A13*A51*A25+A32*A13*A51*A24*A45+A33*A11*A52*A24*A45-A21*A33*A54*A42*A15-A21*A43*A34*A52*A15+A44*A21*A32*A13*A55-A44*A21*A52*A13*A35+A44*A11*A22*A33*A55-A44*A11*A22*A53*A35+A44*A22*A51*A13*A35+A44*A53*A11*A32*A25-A41*A34*A52*A13*A25+A21*A33*A42*A14*A55+A11*A52*A34*A43*A25-A44*A21*A53*A32*A15-A44*A33*A11*A52*A25-A44*A33*A22*A51*A15+A21*A52*A13*A34*A45+A11*A22*A53*A34*A45-A33*A41*A24*A52*A15-A54*A11*A22*A33*A45+A54*A32*A13*A41*A25+A54*A33*A22*A41*A15+A54*A11*A22*A43*A35-A42*A13*A24*A51*A35-A41*A32*A14*A53*A25+A21*A33*A44*A52*A15-A42*A14*A21*A53*A35-A42*A14*A33*A51*A25+A41*A34*A22*A13*A55-A41*A34*A53*A22*A15-A54*A21*A32*A13*A45-A43*A22*A14*A51*A35+A33*A22*A51*A14*A45+A53*A22*A14*A41*A35+A51*A32*A14*A43*A25-A21*A43*A32*A14*A55-A54*A22*A41*A13*A35-A54*A43*A11*A32*A25+A54*A33*A11*A42*A25+A33*A51*A24*A42*A15-A53*A24*A11*A32*A45+A21*A53*A34*A42*A15+A54*A21*A42*A13*A35+A54*A21*A43*A32*A15+A43*A24*A11*A32*A55-A43*A24*A51*A32*A15-A43*A24*A11*A52*A35+A23*A11*A42*A34*A55+A31*A44*A52*A13*A25+A31*A43*A24*A52*A15+A31*A42*A13*A24*A55;
	double psiNum=-B1*A31*A55*A22*A44+B1*A24*A45*A51*A32-B1*A31*A24*A45*A52+B1*A31*A45*A54*A22-B4*A12*A51*A34*A25+B4*A24*A55*A11*A32+B4*A12*A55*A21*A34+B4*A12*A51*A35*A24-B4*A12*A54*A21*A35-B4*A31*A12*A55*A24+B4*A31*A12*A54*A25+B4*A31*A52*A15*A24-B3*A45*A54*A11*A22-B3*A42*A15*A54*A21+B5*A11*A22*A34*A45+B3*A55*A42*A21*A14+B5*A24*A11*A42*A35-B5*A44*A11*A22*A35-B5*A41*A34*A22*A15-B3*A55*A22*A41*A14+B5*A12*A41*A34*A25-B3*A44*A25*A11*A52+B5*A44*A11*A32*A25-B5*A11*A32*A24*A45+B4*A31*A55*A22*A14-B4*A24*A51*A32*A15+B3*A52*A15*A21*A44-B1*A31*A54*A25*A42+B3*A54*A22*A15*A41-B4*A21*A34*A52*A15-B4*A51*A35*A22*A14+B3*A45*A22*A14*A51+B3*A55*A11*A22*A44-B3*A44*A22*A15*A51-B3*A42*A14*A25*A51-B5*A42*A14*A21*A35-B5*A44*A32*A21*A15+B3*A52*A14*A25*A41-B5*A12*A21*A34*A45+B5*A24*A41*A32*A15-B5*A41*A32*A14*A25-B3*A45*A52*A14*A21+B5*A32*A21*A14*A45+B5*A22*A14*A41*A35+B5*A21*A34*A42*A15+B1*A31*A24*A55*A42+B5*A12*A44*A21*A35-B5*A11*A42*A34*A25-B4*A31*A52*A14*A25-B4*A31*A54*A22*A15+B4*A52*A14*A21*A35+B3*A54*A25*A11*A42-B4*A24*A11*A52*A35-B4*A55*A32*A21*A14+B3*A24*A42*A15*A51+B3*A24*A45*A11*A52+B4*A54*A21*A32*A15-B3*A24*A55*A11*A42-B3*A12*A54*A25*A41+B3*A12*A44*A25*A51-B3*A24*A52*A15*A41-B3*A12*A55*A21*A44+B3*A12*A24*A55*A41+B3*A12*A45*A54*A21+B2*A51*A35*A42*A14-B2*A42*A15*A51*A34+B2*A31*A42*A15*A54+B2*A55*A41*A32*A14-B2*A12*A55*A41*A34+B2*A31*A45*A52*A14-B2*A55*A11*A32*A44-B2*A31*A52*A15*A44+B2*A52*A15*A41*A34-B3*A12*A24*A45*A51+B2*A31*A12*A55*A44-B2*A31*A55*A42*A14+B2*A51*A32*A15*A44+B1*A54*A21*A35*A42+B2*A11*A52*A35*A44-B2*A31*A12*A45*A54-B2*A45*A11*A52*A34-B1*A44*A21*A35*A52+B2*A55*A11*A42*A34-B1*A41*A35*A54*A22+B2*A45*A54*A11*A32-B1*A45*A54*A21*A32-B2*A41*A35*A52*A14-B1*A41*A34*A25*A52+B1*A55*A44*A21*A32-B2*A41*A32*A15*A54+B1*A45*A21*A34*A52+B1*A51*A34*A25*A42-B2*A45*A51*A32*A14+B1*A54*A25*A41*A32-B2*A12*A51*A35*A44-B1*A24*A55*A41*A32-B2*A11*A42*A35*A54+B2*A12*A41*A35*A54+B1*A24*A41*A35*A52+B2*A12*A45*A51*A34-B5*A31*A12*A44*A25-B5*A31*A24*A42*A15-B1*A24*A51*A35*A42+B5*A31*A42*A14*A25-B5*A31*A22*A14*A45+B5*A31*A44*A22*A15+B1*A51*A35*A22*A44-B1*A21*A34*A55*A42+B1*A55*A41*A34*A22-B5*A12*A24*A41*A35-B4*A11*A22*A34*A55+B4*A54*A11*A22*A35-B1*A44*A25*A51*A32+B5*A31*A12*A24*A45+B4*A51*A34*A22*A15-B4*A54*A25*A11*A32-B1*A45*A51*A34*A22+B4*A51*A32*A14*A25+B1*A31*A44*A25*A52+B4*A11*A52*A34*A25;

	out[7] =psiNum/psiDen;
	//%--------- thetaDoubleDot Equation ---------%
	double thetaDen= A23*A41*A32*A14*A55-A23*A51*A34*A42*A15-A11*A42*A34*A53*A25+A31*A12*A23*A44*A55-A31*A12*A23*A54*A45+A31*A12*A54*A43*A25-A31*A12*A43*A24*A55+A31*A12*A53*A24*A45-A31*A12*A44*A53*A25-A31*A23*A44*A52*A15-A31*A23*A42*A14*A55+A31*A23*A54*A42*A15+A31*A23*A52*A14*A45-A31*A54*A43*A22*A15+A31*A54*A22*A13*A45-A31*A52*A13*A24*A45-A31*A52*A14*A43*A25-A31*A44*A22*A13*A55+A31*A44*A53*A22*A15-A31*A54*A42*A13*A25+A31*A43*A22*A14*A55-A31*A53*A22*A14*A45+A31*A42*A14*A53*A25-A31*A53*A24*A42*A15-A23*A54*A11*A42*A35+A23*A42*A14*A51*A35+A23*A44*A11*A52*A35+A23*A54*A11*A32*A45-A23*A54*A41*A32*A15-A23*A11*A52*A34*A45+A23*A41*A34*A52*A15+A23*A44*A51*A32*A15-A23*A51*A32*A14*A45-A23*A52*A14*A41*A35-A23*A44*A11*A32*A55-A12*A23*A41*A34*A55+A12*A23*A51*A34*A45-A12*A23*A44*A51*A35+A12*A23*A54*A41*A35+A12*A43*A24*A51*A35-A12*A21*A53*A34*A45-A12*A21*A33*A44*A55+A12*A33*A41*A24*A55+A12*A21*A43*A34*A55-A12*A53*A24*A41*A35+A12*A44*A21*A53*A35+A12*A44*A33*A51*A25+A12*A21*A33*A54*A45-A12*A54*A33*A41*A25+A12*A41*A34*A53*A25-A12*A33*A51*A24*A45-A12*A54*A21*A43*A35-A12*A51*A34*A43*A25-A21*A42*A13*A34*A55-A11*A22*A43*A34*A55+A21*A53*A32*A14*A45+A52*A14*A33*A41*A25-A33*A22*A41*A14*A55+A52*A13*A24*A41*A35+A52*A14*A21*A43*A35+A53*A24*A41*A32*A15+A51*A34*A43*A22*A15+A53*A24*A11*A42*A35+A51*A34*A42*A13*A25-A21*A33*A52*A14*A45-A32*A13*A41*A24*A55-A33*A11*A42*A24*A55-A51*A34*A22*A13*A45-A44*A32*A13*A51*A25+A32*A13*A51*A24*A45+A33*A11*A52*A24*A45-A21*A33*A54*A42*A15-A21*A43*A34*A52*A15+A44*A21*A32*A13*A55-A44*A21*A52*A13*A35+A44*A11*A22*A33*A55-A44*A11*A22*A53*A35+A44*A22*A51*A13*A35+A44*A53*A11*A32*A25-A41*A34*A52*A13*A25+A21*A33*A42*A14*A55+A11*A52*A34*A43*A25-A44*A21*A53*A32*A15-A44*A33*A11*A52*A25-A44*A33*A22*A51*A15+A21*A52*A13*A34*A45+A11*A22*A53*A34*A45-A33*A41*A24*A52*A15-A54*A11*A22*A33*A45+A54*A32*A13*A41*A25+A54*A33*A22*A41*A15+A54*A11*A22*A43*A35-A42*A13*A24*A51*A35-A41*A32*A14*A53*A25+A21*A33*A44*A52*A15-A42*A14*A21*A53*A35-A42*A14*A33*A51*A25+A41*A34*A22*A13*A55-A41*A34*A53*A22*A15-A54*A21*A32*A13*A45-A43*A22*A14*A51*A35+A33*A22*A51*A14*A45+A53*A22*A14*A41*A35+A51*A32*A14*A43*A25-A21*A43*A32*A14*A55-A54*A22*A41*A13*A35-A54*A43*A11*A32*A25+A54*A33*A11*A42*A25+A33*A51*A24*A42*A15-A53*A24*A11*A32*A45+A21*A53*A34*A42*A15+A54*A21*A42*A13*A35+A54*A21*A43*A32*A15+A43*A24*A11*A32*A55-A43*A24*A51*A32*A15-A43*A24*A11*A52*A35+A23*A11*A42*A34*A55+A31*A44*A52*A13*A25+A31*A43*A24*A52*A15+A31*A42*A13*A24*A55;
	double thetaNum=B1*A43*A25*A51*A32-B1*A31*A55*A23*A42+B1*A31*A25*A53*A42+B1*A31*A43*A55*A22-B1*A31*A53*A45*A22+B2*A32*A13*A51*A45+B2*A53*A41*A32*A15-B2*A43*A51*A32*A15+B2*A12*A55*A33*A41-B2*A43*A11*A52*A35+B2*A43*A55*A11*A32-B2*A12*A33*A51*A45-B2*A12*A53*A41*A35+B2*A31*A42*A13*A55+B2*A12*A51*A35*A43-B3*A43*A55*A11*A22-B3*A43*A52*A15*A21-B3*A25*A53*A11*A42-B3*A22*A13*A45*A51-B3*A23*A42*A15*A51+B4*A31*A12*A55*A23-B2*A31*A53*A42*A15-B3*A25*A52*A13*A41-B4*A31*A55*A22*A13-B2*A31*A52*A13*A45+B4*A31*A22*A15*A53-B4*A21*A35*A52*A13-B4*A55*A11*A32*A23+B4*A12*A21*A35*A53+B4*A25*A53*A11*A32+B4*A55*A11*A22*A33+B4*A51*A35*A22*A13-B4*A11*A22*A35*A53+B5*A12*A43*A31*A25+B2*A31*A12*A53*A45+B1*A31*A23*A45*A52+B1*A43*A21*A35*A52+B2*A31*A52*A15*A43-B2*A31*A12*A55*A43-B1*A23*A45*A51*A32-B1*A41*A23*A35*A52-B1*A31*A43*A25*A52+B3*A12*A25*A53*A41+B3*A43*A25*A11*A52+B2*A33*A51*A42*A15-B2*A32*A13*A41*A55+B1*A51*A35*A23*A42+B3*A52*A13*A45*A21+B1*A53*A41*A35*A22+B3*A43*A22*A15*A51+B3*A55*A22*A41*A13+B3*A55*A11*A42*A23+B3*A53*A42*A15*A21-B1*A21*A35*A53*A42-B3*A55*A42*A21*A13+B3*A25*A42*A13*A51+B3*A52*A15*A41*A23-B3*A23*A45*A11*A52-B4*A31*A12*A25*A53+B1*A55*A21*A33*A42+B3*A53*A45*A11*A22-B4*A12*A55*A21*A33+B4*A12*A25*A33*A51+B4*A31*A25*A52*A13-B3*A22*A15*A53*A41+B4*A52*A15*A21*A33-B4*A31*A52*A15*A23-B4*A12*A51*A35*A23-B1*A21*A33*A45*A52-B2*A53*A11*A32*A45-B2*A33*A41*A52*A15-B2*A33*A11*A42*A55-B4*A25*A33*A11*A52-B4*A25*A32*A13*A51+B2*A33*A11*A52*A45-B2*A42*A13*A51*A35+B4*A55*A32*A21*A13+B3*A12*A23*A45*A51-B3*A12*A53*A45*A21-B5*A12*A31*A23*A45-B5*A43*A22*A31*A15-B4*A21*A32*A15*A53-B4*A22*A15*A33*A51+B4*A51*A32*A15*A23-B3*A12*A55*A41*A23-B3*A12*A43*A25*A51+B3*A12*A43*A55*A21+B2*A52*A13*A41*A35+B2*A53*A11*A42*A35+B5*A33*A22*A41*A15+B5*A32*A13*A41*A25-B5*A22*A41*A13*A35-B5*A21*A32*A13*A45-B5*A11*A22*A33*A45+B5*A21*A42*A13*A35+B1*A25*A33*A41*A52-B5*A11*A42*A23*A35+B5*A33*A11*A42*A25+B5*A11*A32*A23*A45+B1*A55*A41*A23*A32-B5*A41*A23*A32*A15+B5*A11*A22*A43*A35-B5*A43*A11*A32*A25+B5*A21*A43*A32*A15-B5*A21*A33*A42*A15+B1*A33*A51*A45*A22-B1*A25*A33*A51*A42+B5*A12*A41*A23*A35-B1*A25*A53*A41*A32-B5*A12*A33*A41*A25+B1*A53*A45*A21*A32+B5*A12*A21*A33*A45-B5*A42*A13*A31*A25+B5*A22*A31*A13*A45-B1*A43*A51*A35*A22+B5*A31*A23*A42*A15-B1*A55*A33*A41*A22-B5*A12*A21*A43*A35+B4*A11*A52*A35*A23-B1*A43*A55*A21*A32;

	out[8] =  thetaNum/thetaDen;
	//%----------phiDoubleDot Equation -----------%
	double phiDen=A23*A41*A32*A14*A55-A23*A51*A34*A42*A15-A11*A42*A34*A53*A25+A31*A12*A23*A44*A55-A31*A12*A23*A54*A45+A31*A12*A54*A43*A25-A31*A12*A43*A24*A55+A31*A12*A53*A24*A45-A31*A12*A44*A53*A25-A31*A23*A44*A52*A15-A31*A23*A42*A14*A55+A31*A23*A54*A42*A15+A31*A23*A52*A14*A45-A31*A54*A43*A22*A15+A31*A54*A22*A13*A45-A31*A52*A13*A24*A45-A31*A52*A14*A43*A25-A31*A44*A22*A13*A55+A31*A44*A53*A22*A15-A31*A54*A42*A13*A25+A31*A43*A22*A14*A55-A31*A53*A22*A14*A45+A31*A42*A14*A53*A25-A31*A53*A24*A42*A15-A23*A54*A11*A42*A35+A23*A42*A14*A51*A35+A23*A44*A11*A52*A35+A23*A54*A11*A32*A45-A23*A54*A41*A32*A15-A23*A11*A52*A34*A45+A23*A41*A34*A52*A15+A23*A44*A51*A32*A15-A23*A51*A32*A14*A45-A23*A52*A14*A41*A35-A23*A44*A11*A32*A55-A12*A23*A41*A34*A55+A12*A23*A51*A34*A45-A12*A23*A44*A51*A35+A12*A23*A54*A41*A35+A12*A43*A24*A51*A35-A12*A21*A53*A34*A45-A12*A21*A33*A44*A55+A12*A33*A41*A24*A55+A12*A21*A43*A34*A55-A12*A53*A24*A41*A35+A12*A44*A21*A53*A35+A12*A44*A33*A51*A25+A12*A21*A33*A54*A45-A12*A54*A33*A41*A25+A12*A41*A34*A53*A25-A12*A33*A51*A24*A45-A12*A54*A21*A43*A35-A12*A51*A34*A43*A25-A21*A42*A13*A34*A55-A11*A22*A43*A34*A55+A21*A53*A32*A14*A45+A52*A14*A33*A41*A25-A33*A22*A41*A14*A55+A52*A13*A24*A41*A35+A52*A14*A21*A43*A35+A53*A24*A41*A32*A15+A51*A34*A43*A22*A15+A53*A24*A11*A42*A35+A51*A34*A42*A13*A25-A21*A33*A52*A14*A45-A32*A13*A41*A24*A55-A33*A11*A42*A24*A55-A51*A34*A22*A13*A45-A44*A32*A13*A51*A25+A32*A13*A51*A24*A45+A33*A11*A52*A24*A45-A21*A33*A54*A42*A15-A21*A43*A34*A52*A15+A44*A21*A32*A13*A55-A44*A21*A52*A13*A35+A44*A11*A22*A33*A55-A44*A11*A22*A53*A35+A44*A22*A51*A13*A35+A44*A53*A11*A32*A25-A41*A34*A52*A13*A25+A21*A33*A42*A14*A55+A11*A52*A34*A43*A25-A44*A21*A53*A32*A15-A44*A33*A11*A52*A25-A44*A33*A22*A51*A15+A21*A52*A13*A34*A45+A11*A22*A53*A34*A45-A33*A41*A24*A52*A15-A54*A11*A22*A33*A45+A54*A32*A13*A41*A25+A54*A33*A22*A41*A15+A54*A11*A22*A43*A35-A42*A13*A24*A51*A35-A41*A32*A14*A53*A25+A21*A33*A44*A52*A15-A42*A14*A21*A53*A35-A42*A14*A33*A51*A25+A41*A34*A22*A13*A55-A41*A34*A53*A22*A15-A54*A21*A32*A13*A45-A43*A22*A14*A51*A35+A33*A22*A51*A14*A45+A53*A22*A14*A41*A35+A51*A32*A14*A43*A25-A21*A43*A32*A14*A55-A54*A22*A41*A13*A35-A54*A43*A11*A32*A25+A54*A33*A11*A42*A25+A33*A51*A24*A42*A15-A53*A24*A11*A32*A45+A21*A53*A34*A42*A15+A54*A21*A42*A13*A35+A54*A21*A43*A32*A15+A43*A24*A11*A32*A55-A43*A24*A51*A32*A15-A43*A24*A11*A52*A35+A23*A11*A42*A34*A55+A31*A44*A52*A13*A25+A31*A43*A24*A52*A15+A31*A42*A13*A24*A55;
	double phiNum= B2*A31*A44*A52*A13-B2*A31*A54*A42*A13-B2*A31*A12*A44*A53+B2*A31*A42*A14*A53+B2*A31*A12*A54*A43+B1*A21*A53*A34*A42-B2*A31*A52*A14*A43-B1*A21*A43*A34*A52-B1*A21*A33*A54*A42-B1*A44*A21*A53*A32-B2*A12*A54*A33*A41-B1*A44*A33*A51*A22+B2*A12*A41*A34*A53+B3*A12*A23*A54*A41-B3*A12*A23*A44*A51-B2*A11*A42*A34*A53+B2*A54*A33*A11*A42+B2*A54*A32*A13*A41-B2*A41*A32*A14*A53+B2*A11*A52*A34*A43-B2*A54*A43*A11*A32+B2*A51*A34*A42*A13+B2*A52*A14*A33*A41+B2*A51*A32*A14*A43-B2*A42*A14*A33*A51+B3*A23*A42*A14*A51+B3*A23*A44*A11*A52-B2*A44*A33*A11*A52-B2*A41*A34*A52*A13+B3*A12*A43*A24*A51-B3*A23*A54*A11*A42-B2*A12*A51*A34*A43-B3*A12*A54*A21*A43-B3*A12*A53*A24*A41-B3*A54*A22*A41*A13+B3*A12*A44*A21*A53-B1*A43*A24*A51*A32+B3*A54*A11*A22*A43+B3*A53*A22*A14*A41-B3*A23*A52*A14*A41+B2*A12*A44*A33*A51-B3*A42*A13*A24*A51+B3*A52*A14*A21*A43-B2*A44*A32*A13*A51+B3*A52*A13*A24*A41+B1*A31*A23*A54*A42-B3*A42*A14*A21*A53-B3*A43*A24*A11*A52+B3*A44*A22*A51*A13-B3*A43*A22*A14*A51+B4*A12*A51*A23*A34-B4*A12*A21*A53*A34-B3*A44*A11*A22*A53-B4*A53*A22*A31*A14-B4*A52*A13*A31*A24+B4*A12*A53*A31*A24+B4*A31*A23*A52*A14+B4*A22*A31*A13*A54-B3*A44*A52*A21*A13-B4*A12*A31*A23*A54+B3*A54*A42*A21*A13+B3*A53*A24*A11*A42+B4*A33*A11*A52*A24-B4*A21*A32*A13*A54+B4*A21*A52*A13*A34+B4*A11*A22*A53*A34+B4*A33*A22*A51*A14-B4*A21*A33*A52*A14-B4*A11*A22*A33*A54-B4*A53*A11*A32*A24+B4*A21*A53*A32*A14-B4*A22*A51*A13*A34-B4*A11*A52*A23*A34-B4*A51*A23*A32*A14+B4*A11*A32*A23*A54+B4*A12*A21*A33*A54-B4*A12*A33*A51*A24-B5*A12*A43*A31*A24-B5*A12*A41*A23*A34+B5*A21*A33*A42*A14-B5*A21*A42*A13*A34+B5*A12*A31*A23*A44+B5*A12*A33*A41*A24-B5*A22*A31*A13*A44-B5*A33*A11*A42*A24+B4*A32*A13*A51*A24-B5*A12*A21*A33*A44+B5*A12*A21*A43*A34+B5*A11*A22*A33*A44-B5*A32*A13*A41*A24+B5*A43*A11*A32*A24+B5*A41*A23*A32*A14-B5*A11*A32*A23*A44+B5*A22*A41*A13*A34+B5*A21*A32*A13*A44-B5*A31*A23*A42*A14+B5*A42*A13*A31*A24+B5*A11*A42*A23*A34-B5*A11*A22*A43*A34-B5*A33*A22*A41*A14+B5*A43*A22*A31*A14-B5*A21*A43*A32*A14+B2*A44*A53*A11*A32+B1*A33*A51*A24*A42+B1*A53*A24*A41*A32+B1*A21*A33*A44*A52+B1*A51*A34*A43*A22-B1*A31*A53*A24*A42-B1*A31*A23*A44*A52+B1*A31*A44*A53*A22+B1*A31*A43*A24*A52-B1*A31*A54*A43*A22+B1*A54*A21*A43*A32-B1*A33*A41*A24*A52-B1*A41*A34*A53*A22-B1*A23*A54*A41*A32-B1*A23*A51*A34*A42+B1*A23*A44*A51*A32+B1*A23*A41*A34*A52+B1*A54*A33*A41*A22;

	out[9]=phiNum/phiDen;

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
        	
       	 	double u1 = y[0];
        	double u2 = y[1];
        	double psi = y[2];
        	double theta = y[3];
        	double phi = y[4];
        	double u1Dot = y[5];
        	double u2Dot = y[6];
        	double psiDot = y[7];
        	double thetaDot = y[8];
        	double phiDot = y[9];
        	
        	double theta2 = Math.asin(u2/L);
             		theta2=-theta2;
           	double theta1 = Math.asin(u1/L);
           	
           	// Extract quaternion values
        	// See Bong wie (p.318)
        	// Euler axis
        	double e1=0;
        	double e2=0;
        	double e3=1;
           	
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
        	
        	
        	// Calculate quaternions using the method on Bong Wie (p.321)
        	double q1 = Math.sin(phi/2)*Math.cos(theta/2)*Math.cos(psi/2) - Math.cos(phi/2)*Math.sin(theta/2)*Math.sin(psi/2);
        	double q2 = Math.cos(phi/2)*Math.sin(theta/2)*Math.cos(psi/2) + Math.sin(phi/2)*Math.cos(theta/2)*Math.sin(psi/2);
        	double q3 = Math.cos(phi/2)*Math.cos(theta/2)*Math.sin(psi/2) - Math.sin(phi/2)*Math.sin(theta/2)*Math.cos(psi/2);
        	double q4 = Math.cos(phi/2)*Math.cos(theta/2)*Math.cos(psi/2) + Math.sin(phi/2)*Math.sin(theta/2)*Math.sin(psi/2);
        	double quat_check = q1*q1+q2*q2+q3*q3+q4*q4;
        
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
        	
        	// Angular velocities in the body frame
        	double wx = phiDot-psiDot*Math.sin(theta);
   			double wy = thetaDot*Math.cos(phi) + psiDot*Math.cos(theta)*Math.sin(phi);
   			double wz = psiDot*Math.cos(theta)*Math.cos(phi) - thetaDot*Math.sin(phi);
			
			//Energy Calculation
        		double Potential = 3*EI*(u1*u1)/(2*(L*L*L)) + 3*EI*(u2*u2)/(2*(L*L*L));
        		double Potential1 = 3*EI*(u1*u1)/(2*(L*L*L));
        		double Potential2 = 3*EI*(u2*u2)/(2*(L*L*L));
				double Kinetic = (0.5)*(Ixx*(wx*wx) + Iyy*(wy*wy) + Izz*(wz*wz)) + (0.5)*m*((-u1*wz)*(-u1*wz) + ((a+L)*wz+u1Dot)*((a+L)*wz+u1Dot) + (u1*wx-(a+L)*wy)*(u1*wx-(a+L)*wy)) + (0.5)*m*((-u2*wz)*(-u2*wz)+ (-(a+L)*wz+u2Dot)*(-(a+L)*wz+u2Dot)+(u2*wx+(a+L)*wy)*(u2*wx+(a+L)*wy));
  				double KineticSpacecraft = (0.5)*(Ixx*(wx*wx) + Iyy*(wy*wy) + Izz*(wz*wz));
  				double KineticParticle1 = (0.5)*m*((-u1*wz)*(-u1*wz) + ((a+L)*wz+u1Dot)*((a+L)*wz+u1Dot) + (u1*wx-(a+L)*wy)*(u1*wx-(a+L)*wy));
  				double KineticParticle2 = (0.5)*m*((-u2*wz)*(-u2*wz)+ (-(a+L)*wz+u2Dot)*(-(a+L)*wz+u2Dot)+(u2*wx+(a+L)*wy)*(u2*wx+(a+L)*wy));
  			// You can't use (1/2); you have to use (0.5)! WHY??
  				
  				
  				double Energy = Potential + Kinetic;
        	
        	// Build 3 x 4 Q matrix
    		double[][] array = new double[3][3];
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
    		
    		//Construct Angular velocity vector
    		double [] angularRateArray = new double[3];
    		angularRateArray[0] = wx;
    		angularRateArray[1] = wy;
    		angularRateArray[2] = wz;
    		Matrix angularRateVector = new Matrix(angularRateArray, 3); //Column vector
    		
    		//Transform angular rates in body frame to inertial frame
    		Matrix inertialAngularRateVector = new Matrix(T_transpose.times(angularRateVector));
    		double wx_i = inertialAngularRateVector.get(0,0);
    		double wy_i = inertialAngularRateVector.get(1,0);
    		double wz_i = inertialAngularRateVector.get(2,0);
    		        	
        	
  			
  			//Angular Momentum Calculation
  				double Hsc_x=Ixx*wx_i; // Angular momentum in a frame
   				double Hsc_y=Iyy*wy_i; //Angular momentum in a frame
   				double Hsc_z=Izz*wz_i; // Angular momentum in a frame
   
   				//-------------Tip mass 1 position vector-------------------------------
   				//r1b1=(a+L);
   				//r1b2=u1;
   				//r1b3=0;
   				//r1b = [r1b1 r1b2 r1b3]';
   				//Construct r1 position vector
    			double [] r1PositionArray = new double[3];
    			r1PositionArray[0] = a+L;
    			r1PositionArray[1] = u1;
    			r1PositionArray[2] = 0;
    			Matrix r1Position_body = new Matrix(r1PositionArray, 3); //Column vector
   				Matrix r1Position_inertial = new Matrix(T_transpose.times(r1Position_body)); // Column vector
  		 
   				//r1i_skew = [0        -r1i(3)     r1i(2);
            	//			r1i(3)     0          -r1i(1);
            	//   			-r1i(2)    r1i(1)      0 ];
   				double[][] r1_skewArray = new double[3][3];
   				// Row1
    			r1_skewArray[0][0] = 0;
    			r1_skewArray[0][1] = -r1Position_inertial.get(2,0);
    			r1_skewArray[0][2] = r1Position_inertial.get(1,0);
    			// Row2
    			r1_skewArray[1][0] = r1Position_inertial.get(2,0);
    			r1_skewArray[1][1] = 0;
    			r1_skewArray[1][2] = -r1Position_inertial.get(0,0);
    			// Row3
    			r1_skewArray[2][0] = -r1Position_inertial.get(1,0);
    			r1_skewArray[2][1] = r1Position_inertial.get(2,0);
    			r1_skewArray[2][2] = 0;
    			Matrix r1_skewMatrix = new Matrix(r1_skewArray, 3, 3);
    		
    			//-------------Tip mass 2 position vector-------------------------------
   				//r2b1=-(a+L);
   				//r2b2=u2;
   				//r2b3=0;
   				//r2b=[r2b1 r2b2 r2b3]';
   				//Construct r2 position vector
   				double [] r2PositionArray = new double[3];
   				r2PositionArray[0] = -(a+L);
   				r2PositionArray[1] = u2;
   				r2PositionArray[2] = 0;
   				Matrix r2Position_body = new Matrix(r2PositionArray, 3);
   				Matrix r2Position_inertial = new Matrix(T_transpose.times(r2Position_body));
   		
   				//r2i_skew = [0      -r2i(3)      r2i(2);
     			//			 r2i(3)    0          -r2i(1);
      			//			-r2i(2)    r2i(1)       0];
   				double[][] r2_skewArray = new double[3][3];
   				// Row1
    			r2_skewArray[0][0] = 0;
    			r2_skewArray[0][1] = -r2Position_inertial.get(2,0);
    			r2_skewArray[0][2] = r2Position_inertial.get(1,0);
    			// Row2
    			r2_skewArray[1][0] = r2Position_inertial.get(2,0);
    			r2_skewArray[1][1] = 0;
    			r2_skewArray[1][2] = -r2Position_inertial.get(0,0);
    			// Row3
    			r2_skewArray[2][0] = -r2Position_inertial.get(1,0);
    			r2_skewArray[2][1] = r2Position_inertial.get(2,0);
    			r2_skewArray[2][2] = 0;
    			Matrix r2_skewMatrix = new Matrix(r2_skewArray, 3, 3);
    		
    			//-------------Tip mass 1 velocity vector-------------------------------
  				//r1Dot_b1= -Wb(3)*u1;
   				//r1Dot_b2= u1_dot+Wb(3)*(a+L);
   				//r1Dot_b3= u1*Wb(1)-Wb(2)*(a+L);
   				//r1Dot_b = [r1Dot_b1  r1Dot_b2  r1Dot_b3]';
   				//r1Dot_i = T'*r1Dot_b;
   				double [] r1VelocityArray = new double[3];
   				r1VelocityArray[0] = -wz*u1;
   				r1VelocityArray[1] = u1Dot+wz*(a+L);
   				r1VelocityArray[2] = u1*wx-wy*(a+L);
   				Matrix r1Velocity_body = new Matrix(r1VelocityArray, 3);
   				Matrix r1Velocity_inertial = new Matrix(T_transpose.times(r1Velocity_body));
   			
   				//-------------Tip mass 2 velocity vector-------------------------------
   				//r2Dot_b1 = -Wb(3)*u2;
   				//r2Dot_b2 = u2_dot - Wb(3)*(a+L);
   				//r2Dot_b3 = Wb(1) + Wb(2)*(a+L);
   				//r2Dot_b = [r2Dot_b1  r2Dot_b2  r2Dot_b3]';
   				//r2Dot_i = T'*r2Dot_b;
				double [] r2VelocityArray = new double[3];
				r2VelocityArray[0] = -wz*u2;    
				r2VelocityArray[1] = u2Dot - wz*(a+L);
				r2VelocityArray[2] = wx + wy*(a+L);
				Matrix r2Velocity_body = new Matrix(r2VelocityArray, 3);
				Matrix r2Velocity_inertial = new Matrix(T_transpose.times(r2Velocity_body));
			
   				//----------- Take Cross products ---------------------------------------
   				//Hm1=r1i_skew*(m*r1Dot_i);
   				Matrix Hm1 = new Matrix(r1_skewMatrix.times(r1Velocity_inertial.times(m)));
   				
   				double Hm1_a1 = -r1Position_inertial.get(2,0)*m*r1Velocity_inertial.get(1,0) + r1Position_inertial.get(1,0)*m*r1Velocity_inertial.get(2,0);
   				double Hm1_a2 = r1Position_inertial.get(2,0)*m*r1Velocity_inertial.get(0,0) -r1Position_inertial.get(0,0)*m*r1Velocity_inertial.get(2,0);
   				double Hm1_a3 = -r1Position_inertial.get(1,0)*m*r1Velocity_inertial.get(0,0) + r1Position_inertial.get(0,0)*m*r1Velocity_inertial.get(1,0);
   				
   				
   				//Hm2=r2i_skew*(m*r2Dot_i);
   				Matrix Hm2 = new Matrix(r2_skewMatrix.times(r2Velocity_inertial.times(m)));
   				
   				double Hm2_a1 = -r2Position_inertial.get(2,0)*m*r2Velocity_inertial.get(1,0) + r2Position_inertial.get(1,0)*m*r2Velocity_inertial.get(2,0);
   				double Hm2_a2 = r2Position_inertial.get(2,0)*m*r2Velocity_inertial.get(0,0) -r2Position_inertial.get(0,0)*m*r2Velocity_inertial.get(2,0);
   				double Hm2_a3 = -r2Position_inertial.get(1,0)*m*r2Velocity_inertial.get(0,0) + r2Position_inertial.get(0,0)*m*r2Velocity_inertial.get(1,0);
   				
   				
   				//double Htotal_1 =Hsc_x + Hm1.get(0,0) + Hm2.get(0,0);
   				//double Htotal_2 =Hsc_y + Hm1.get(1,0) + Hm2.get(1,0);
   				//double Htotal_3 =Hsc_z + Hm1.get(2,0) + Hm2.get(2,0);
   				
   				double Htotal_1 =Hsc_x + Hm1_a1 + Hm2_a1;
   				double Htotal_2 =Hsc_y + Hm1_a2 + Hm2_a2;
   				double Htotal_3 =Hsc_z + Hm1_a3 + Hm2_a3;
   				double Htotal = (Htotal_1)*(Htotal_1)+ (Htotal_2)*(Htotal_2) + (Htotal_3)*(Htotal_3);
       	
        	/*
			private FourPlots uPlot = new FourPlots();								
			private ThreePlots eulerPlot = new ThreePlots();
			private SinglePlot quarternion_check = new SinglePlot();
			private FourPlots quatPlot = new FourPlots();
			private TwoPlots EandHPlot = new TwoPlots();
			private ThreePlots eulerRatePlot = new ThreePlots();
			private ThreePlots angularPlot = new ThreePlots();
			private ThreePlots HPlot = new ThreePlots();
			*/
        	
        	uPlot.firstPlot.addPoint(0, t, u1, first);
        	uPlot.secondPlot.addPoint(0,t, u2, first);
        	uPlot.thirdPlot.addPoint(0, t, u1Dot, first);
        	uPlot.fourthPlot.addPoint(0, t, u2Dot, first);
        	
        	eulerPlot.topPlot.addPoint(0, t, psi, first);
        	eulerPlot.middlePlot.addPoint(0, t, theta, first);
        	eulerPlot.bottomPlot.addPoint(0, t, phi, first);
        	
        	quaternion_check.plot.addPoint(0, t, quat_check, first);
        	
        	quatPlot.firstPlot.addPoint(0, t, q1, first);
        	quatPlot.secondPlot.addPoint(0, t, q2, first);
        	quatPlot.thirdPlot.addPoint(0, t, q3, first);
        	quatPlot.fourthPlot.addPoint(0, t, q4, first);
        	
        	energyPlot.topPlot.addPoint(0, t, Potential, first);
        	energyPlot.middlePlot.addPoint(0, t, Kinetic, first);
        	energyPlot.bottomPlot.addPoint(0, t, Energy, first);
        	
        	eulerRatePlot.topPlot.addPoint(0, t, psiDot, first);
        	eulerRatePlot.middlePlot.addPoint(0, t, thetaDot, first);
        	eulerRatePlot.bottomPlot.addPoint(0, t, phiDot, first);
        	
        	angularPlot.topPlot.addPoint(0, t, wx, first);
        	angularPlot.middlePlot.addPoint(0, t, wy, first);
        	angularPlot.bottomPlot.addPoint(0, t, wz, first);
        	
        	HPlot.firstPlot.addPoint(0, t, Htotal_1, first);
        	HPlot.secondPlot.addPoint(0, t, Htotal_2, first);
        	HPlot.thirdPlot.addPoint(0, t, Htotal_3, first);
        	HPlot.fourthPlot.addPoint(0,t,Htotal, first);
        	
        	potentialPlot.topPlot.addPoint(0,t, Potential1, first);
        	potentialPlot.bottomPlot.addPoint(0, t, Potential2, first);
        	
        	kineticPlot.topPlot.addPoint(0, t, KineticSpacecraft, first);
        	kineticPlot.middlePlot.addPoint(0, t, KineticParticle1, first);
        	kineticPlot.bottomPlot.addPoint(0, t, KineticParticle2, first);
        	
        	H1_BreakDown.topPlot.addPoint(0, t, Hsc_x, first);
        	H1_BreakDown.middlePlot.addPoint(0, t, Hm1_a1, first);
        	H1_BreakDown.bottomPlot.addPoint(0, t, Hm2_a1, first);
        	
        	H2_BreakDown.topPlot.addPoint(0, t, Hsc_y, first);
        	H2_BreakDown.middlePlot.addPoint(0, t, Hm1_a2,first);
        	H2_BreakDown.bottomPlot.addPoint(0, t, Hm2_a2, first);
        	
        	H3_BreakDown.topPlot.addPoint(0, t, Hsc_z, first);
        	H3_BreakDown.middlePlot.addPoint(0, t, Hm1_a3, first);
        	H3_BreakDown.bottomPlot.addPoint(0, t, Hm2_a3,first);
        	
        	
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
        	//System.out.println("m= "+m);//+" "+y[1]+" "+y[2]);
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
   		
		uPlot.setVisible(true);
		eulerPlot.setVisible(true);
		//quaternion_check.setVisible(true);
		//quatPlot.setVisible(true);
		energyPlot.setVisible(true);
		//eulerRatePlot.setVisible(true);
		//angularPlot.setVisible(true);
		HPlot.setVisible(true);
		//potentialPlot.setVisible(true);
		//kineticPlot.setVisible(true);
		//H1_BreakDown.setVisible(true);
		//H2_BreakDown.setVisible(true);
		//H3_BreakDown.setVisible(true);
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
        RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78(1E-12);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
		float quatBeam1[][]= new  float[5][numberOfPts+1];
		float quatBeam2[][]= new  float[5][numberOfPts+1];
		 
        // create an instance
        FlexibleThreeD si = new FlexibleThreeD(time_step,  quat_values,quatBeam1, quatBeam2);
		
        // initialize the variables
        double [] x0 = new double[10];
        x0[0] = 0.5;
        x0[1] = 0.1;
        x0[2] = 0.0;
        x0[3] = 0.0;
        x0[4] = 0.0;
        x0[5] = 0.0;
        x0[6] = 0.0;
        x0[7] = 0.0;
        x0[8] = 0.0;
        x0[9] = 0.0;
              

        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        //rk78.integrate(t0,x0,tf,si,si,true);
        // make the plot visible
        si.makePlotsVisible();
        
    }
}
