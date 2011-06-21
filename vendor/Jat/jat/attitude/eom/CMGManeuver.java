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
import jat.matvec.data.*;
import jat.plot.*;
import jat.attitude.*;

/**
 * <P>
 * This class contains the equations of motion for the simulation of a rigid spacecraft
 * performing attitude maneuver using 3 single-gimbal control moment gyros.
 * The class implements jat.alg.integrators.Derivatives and
 * jat.alg.integrateors.Printble, so an outside application code can perform a numerical simulation
 * of equations of motion defined in this class and get the output from the simulation.
 * 
 * @author Noriko Takada
 * @version Last Modified: 08/15/2004 
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion 
 * 
 */
public class CMGManeuver implements EquationsOfMotion
{
 //		time_step:			Time step of numerical integration
 //		quat_values[][]		Two dimensional array that contains quarternions from simulation
 //		M1,M2,M3			External Torques about 1,2,and 3 axis
 //		I1,I2,I3			Principal moments of inertia
 //		J1,J2,J3			Moments of inertia of CMG1, CMG2, CMG3
 //		A1,A2,A3			Nomial rotor speed of CMG1, CMG2, CMG3
 //		Kx Ky Kz 	Proportional gains
 //		Kxd Kyd Kzd	Derivative gains
 //		psi, theta, phi	Amount of attitude maneuver in terms of euler angles
 //		q
 //		rotation_plot		Plots of angular velocities
 //		angle_plot			Plots of euler angles
 //		quarternion_check	Plot of e1^2 + e2^2 + e3^2 +e4^2
 //		CMGPlot				Plot of CMGs rotor speeds
	
	double time_step = 0.1;
	private float quat_values[][];
	
	// Create variables for the necessary plots
	private ThreePlots gimbal_angle_plot = new ThreePlots();
	private ThreePlots gimbal_rate_plot = new ThreePlots();								
	private ThreePlots angular_velocity_plot = new ThreePlots();
	private SinglePlot quarternion_check = new SinglePlot();
	private ThreePlots euler_angle_plot = new ThreePlots();
	private FourPlots angular_momentum_plot = new FourPlots();
	private SinglePlot energy_plot = new SinglePlot();
	private FourPlots quaternion_plot = new FourPlots();
	
	private double I1 = 10.42;       
    private double I2 = 35.42;
    private double I3 = 41.67;
    
    private double M1 = 0;
    private double M2 = 0;
    private double M3 = 0;
    
    double J1 = 10;
    double J2 = 10;
    double J3 = 10;
     
    double OMEGA1=500;
    double OMEGA2=500;
    double OMEGA3=500;
    							  
	double Kx=-I1;
	double Kxd=4*I1;
	double Ky=-I2;
	double Kyd=4*I2;
	double Kz=-I3;
	double Kzd=4*I3;
	
	// Amount of attitude maneuver in terms of Euler angles
	// These have to be converted into target quarternions.
	private double psi=10*Math.PI/180;
	private double theta=60*Math.PI/180;
	private double phi=30*Math.PI/180;
	
	// Target Quarternions
	double[] q = convertEulerToTargetQ(psi, theta, phi);


	/**
	 * Constructor 1
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 */
	public CMGManeuver(double time_step, float quat_values[][])
	{
		// Setup plots
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
	}
	
	
	/** 
	 *	Constructor 2
	 *	@param				time_step:			Time step of numerical integration
	 * 	@param				quat_values[][]		Two dimensional array that contains quarternions from simulation
	 *  @param I1 I2 I3	Principle moments of inertia about 1, 2, and 3 axis respectively
	 *  @param J			Moment of inertia of CMG
	 * 	@param angle		Orientation angle of CMG
	 * 	@param psi			Target euler angle
	 *  @param phi			Target euler angle
	 *  @param theta		Target euler angle
	 */
	public CMGManeuver(double time_step,  double I1, double I2, double I3
							, double J, double A, double psi, double phi, double theta,
							float quat_values[][] )
	{
		setupPlots();
		
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
		
		this.J1 = J;
		this.J2 = J;
		this.J3 = J;
		this.OMEGA1 = A;
		this.OMEGA2 = A;
		this.OMEGA3 = A;
		
		this.psi = psi;
		this.phi = phi;
		this.theta = theta;
		// a,b, and c are angles in degrees.
		 psi=psi*Math.PI/180;
		 theta=theta*Math.PI/180;
		 phi=phi*Math.PI/180;
		
		q = convertEulerToTargetQ(psi, theta, phi);
	}
	
	/**
	 * Constructor 3
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 * @param 		a		Target psi angle
	 * @param 		b		Target theta angle
	 * @param 		c 		Target phi angle
	 */
	public CMGManeuver(double time_step, double a, double b, double c,
						float quat_values[][])
	{
		 setupPlots();
		 this.time_step = time_step;
		 this.quat_values = quat_values;
		 // a,b, and c are angles in degrees.
		 psi=a*Math.PI/180;
		 theta=b*Math.PI/180;
		 phi=c*Math.PI/180;
		
		q = convertEulerToTargetQ(psi, theta, phi);
	}

	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		// FourPlots wheel_velocity_plot = new FourPlots();								
		// ThreePlots angullar_velocity_plot = new ThreePlots();
		// SinglePlot quarternion_check = new SinglePlot();
		// ThreePlots euler_angle_plot = new ThreePlots();
		// FourPlots angular_momentum_plot = new FourPlots();
		// SinglePlot energy_plot = new SinglePlot();
		
		
		// Setup plots
		gimbal_rate_plot.setTitle("CMG Gimbal Rates");
		gimbal_rate_plot.topPlot.setXLabel("t (sec)");
		gimbal_rate_plot.topPlot.setYLabel("CMG1 Gimbal rate (rad/sec)");
		gimbal_rate_plot.middlePlot.setXLabel("t (sec)");
		gimbal_rate_plot.middlePlot.setYLabel("CMG2 Gimbal rate (rad/sec)");
		gimbal_rate_plot.bottomPlot.setXLabel("t (sec)");
		gimbal_rate_plot.bottomPlot.setYLabel("CMG3 Gimbal rate (rad/sec)");
		
		gimbal_angle_plot.setTitle("CMG Gimbal Angles");
		gimbal_angle_plot.topPlot.setXLabel("t (sec)");
		gimbal_angle_plot.topPlot.setYLabel("CMG1 Gimbal angle (rad)");
		gimbal_angle_plot.middlePlot.setXLabel("t (sec)");
		gimbal_angle_plot.middlePlot.setYLabel("CMG2 Gimbal angle (rad)");
		gimbal_angle_plot.bottomPlot.setXLabel("t (sec)");
		gimbal_angle_plot.bottomPlot.setYLabel("CMG3 Gimbal angle (rad)");
		
		
		angular_velocity_plot.setTitle("Spacecraft Angular Velocities");
		angular_velocity_plot.topPlot.setXLabel("t (sec)");
		angular_velocity_plot.topPlot.setYLabel("w1 (rad/sec)");
		angular_velocity_plot.middlePlot.setXLabel("t (sec)");
		angular_velocity_plot.middlePlot.setYLabel("w2 (rad/sec)");
		angular_velocity_plot.bottomPlot.setXLabel("t (sec)");
		angular_velocity_plot.bottomPlot.setYLabel("w3 (rad/sec)");
		
		quarternion_check.setTitle("Quarternion Check");
        quarternion_check.plot.setXLabel("t(sec)");
        quarternion_check.plot.setYLabel("e1^2 + e2^2 + e3^2 + e4^2");
		
		euler_angle_plot.setTitle("Euler Angles");
        euler_angle_plot.topPlot.setXLabel("t(sec)");
        euler_angle_plot.topPlot.setYLabel("Psi(rad)");
		euler_angle_plot.middlePlot.setXLabel("t(sec)");
		euler_angle_plot.middlePlot.setYLabel("Theta(rad)");
        euler_angle_plot.bottomPlot.setXLabel("t(sec)");
        euler_angle_plot.bottomPlot.setYLabel("Phi(rad)");
		
		angular_momentum_plot.setTitle("Angular Momentum");
		angular_momentum_plot.firstPlot.setXLabel("t (sec)");
		angular_momentum_plot.firstPlot.setYLabel("H (a1)");
		angular_momentum_plot.secondPlot.setXLabel("t (sec)");
		angular_momentum_plot.secondPlot.setYLabel("H (a2)");
		angular_momentum_plot.thirdPlot.setXLabel("t (sec)");
		angular_momentum_plot.thirdPlot.setYLabel("H (a3)");
		angular_momentum_plot.fourthPlot.setXLabel("t (sec)");
		angular_momentum_plot.fourthPlot.setYLabel("H total");
		
		energy_plot.setTitle("RotationalKinetic Energy");
		energy_plot.plot.setXLabel("t (sec)");
		energy_plot.plot.setYLabel("Energy");
		
		quaternion_plot.setTitle("Quaternions");
		quaternion_plot.firstPlot.setXLabel("t (sec)");
		quaternion_plot.firstPlot.setYLabel("q1");
		quaternion_plot.secondPlot.setXLabel("t (sec)");
		quaternion_plot.secondPlot.setYLabel("q2");
		quaternion_plot.thirdPlot.setXLabel("t (sec)");
		quaternion_plot.thirdPlot.setYLabel("q3");
		quaternion_plot.fourthPlot.setXLabel("t (sec)");
		quaternion_plot.fourthPlot.setYLabel("q4");
	}
	
	/** Compute the derivatives.
     * Equations of Motion
     * @params t    double containing time or the independent variable.
     * @params x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
		
	public double[] derivs(double t, double[] x)
    {
       	// Error quarternions (Necessary to determine the required control torque)
       	double qe1 =  q[3]*x[3]  + q[2]*x[4] - q[1]*x[5]  - q[0]*x[6];
		double qe2 = -q[2]*x[3]  + q[3]*x[4] + q[0]*x[5]  - q[1]*x[6];
		double qe3 =  q[1]*x[3]  - q[0]*x[4] + q[3]*x[5]  - q[2]*x[6];
		double qe4 =  q[0]*x[3]  + q[1]*x[4] + q[2]*x[5]  + q[3]*x[6];
       	
       	// Control Torque
        double Mc1 = 2*Kx*qe1*qe4  - Kxd*x[0];
		double Mc2 = 2*Ky*qe2*qe4  - Kyd*x[1];
		double Mc3 = 2*Kz*qe3*qe4  - Kzd*x[2];	
       
       	double w1 = x[0];		
        double w2 = x[1];		
        double w3 = x[2];		
        double q1 = x[3];		
        double q2 = x[4];		
        double q3 = x[5];		
        double q4 = x[6];		
        double delta1 = x[7];		
        double delta2 = x[8];		
        double delta3 = x[9];		
       	
       	
       	double [] out = new double[11];
       	out[0] = ((I2-I3)/I1)*x[1]*x[2] + M1/I1 + Mc1/I1;
       	out[1] = ((I3-I1)/I2)*x[0]*x[2] + M2/I2 + Mc2/I2;
       	out[2] = ((I1-I2)/I3)*x[0]*x[1] + M3/I3 + Mc3/I3;
       	out[3] = -0.5* (-x[2]*x[4] + x[1]*x[5] - x[0]*x[6]);
       	out[4] = -0.5* (x[2]*x[3]  - x[0]*x[5] - x[1]*x[6]);
       	out[5] = -0.5* (-x[2]*x[6] + x[0]*x[4] - x[1]*x[3]);
       	out[6] = -0.5* (x[2]*x[5]  + x[1]*x[4] + x[0]*x[3]);
       	
       
       	// x[7] = Gimbal rate for CMG1
       	// x[8] = Gimbal rate for CMG2
       	// x[9] = Gimbal rate for CMG3
       	
       	Matrix Q = new Matrix(3,3);
       	Matrix P = new Matrix(3,1);
       	Matrix Tc = new Matrix(3,1);
       	Matrix Q_inv = new Matrix(3,3);
       	Matrix Omega = new Matrix(3,1);
       	
       	// Set the matrix elements
       	// Don't forget that the index starts at 0!!
       	Q.set(0,0,0);
       	Q.set(0,1, -J2*(w3+OMEGA2*Math.sin(delta2)));
       	Q.set(0,2, J3*(w2-OMEGA3*Math.cos(delta3)));
       	Q.set(1,0, J1*(w3-OMEGA1*Math.cos(delta1)));
       	Q.set(1,1,0);
       	Q.set(1,2, -J3*(w1+OMEGA3*Math.sin(delta3)));
       	Q.set(2,0,-J1*(w2+OMEGA1*Math.sin(delta1)));
       	Q.set(2,1, J2*(w1-OMEGA2*Math.cos(delta2)));
       	Q.set(2,2,0);
       	
		P.set(0,0, (-J1*OMEGA1*(w2*Math.cos(delta1) + w3*Math.sin(delta1)) + J2*OMEGA2*w2*Math.sin(delta2) + J3*OMEGA3*w3*Math.cos(delta3) - Mc1));
		P.set(1,0,  (J1*OMEGA1*w1*Math.cos(delta1) - J2*OMEGA2*(w3*Math.cos(delta2) + w1*Math.sin(delta2)) + J3*OMEGA3*w3*Math.sin(delta3) -Mc2));
		P.set(2,0, (J1*OMEGA1*w1*Math.sin(delta1) + J2*OMEGA2*w2*Math.cos(delta2) - J3*OMEGA3*(w1*Math.cos(delta3) + w2*Math.sin(delta3))-Mc3));
       	
       	
       	
       	Q_inv=Q.inverse();
       	Omega=Q_inv.times(P);
       	
       	out[7] = Omega.get(0,0);
       	out[8] = Omega.get(1,0);
       	out[9] = Omega.get(2,0);
       	

        // Build  3x3 matrix Q and 3x1 matrix P
            
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
		
        System.out.println(t+" "+y[0]+" "+y[1]+" "+first);

		// Define state variables
		double w1 = y[0];
        double w2 = y[1];
        double w3 = y[2];
        double q1 = y[3];
        double q2 = y[4];
        double q3 = y[5];
        double q4 = y[6];
        double delta1 = y[7];
        double delta2 = y[8];
        double delta3 = y[9];
        
        QuatToDeg tester = new QuatToDeg(q1, q2, q3, q4);
    	
    	double[] angle = new double[3];
    	angle = tester.calculateAngle();
    	double Theta = angle[0];
    	double Psi = angle[1];
    	double Phi = angle[2];
    	
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
    	
    	// Build Transformation Matrix
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
        
   		// ************ Angular Momentum **************************//	
   		//Define Q matrix
   		double Q11 = 0;
  		double Q12 = -J2*(w3+OMEGA2*Math.sin(delta2));
   		double Q13 = J3*(w2-OMEGA3*Math.cos(delta3));
   		double Q21 = J1*(w3-OMEGA1*Math.cos(delta1)) ;
   		double Q22 = 0;
   		double Q23 = -J3*(w1+OMEGA3*Math.sin(delta2));
   		double Q31 = -J1*(w2+OMEGA1*Math.sin(delta1));
   		double Q32 = J2*(w1-OMEGA2*Math.cos(delta2));
   		double Q33 = 0;
   		
   		double qt1 = q[0];
   		double qt2 = q[1];
   		double qt3 = q[2];
   		double qt4 = q[3];
   
   		double qe1 = qt4*q1  +qt3*q2  -qt2*q3  -qt1*q4;
		double qe2 = -qt3*q1   +qt4*q2  +qt1*q3  -qt2*q4;
		double qe3 = qt2*q1  -qt1*q2  +qt4*q3  -qt3*q4;
		double qe4 = qt1*q1   +qt2*q2  +qt3*q3  +qt4*q4;

   		double Mc1=(2*Kx*qe1*qe4  - Kxd*w1);
		double Mc2=(2*Ky*qe2*qe4  - Kyd*w2);
   		double Mc3=(2*Kz*qe3*qe4  - Kzd*w3);
   

		double P1 =  -J1*OMEGA1*(w2*Math.cos(delta1) + w3*Math.sin(delta1)) + J2*OMEGA2*w2*Math.sin(delta2) + J3*OMEGA3*w3*Math.cos(delta3) - Mc1;
		double P2 = J1*OMEGA1*w1*Math.cos(delta1) - J2*OMEGA2*(w3*Math.cos(delta2) + w1*Math.sin(delta2)) + J3*OMEGA3*w3*Math.sin(delta3) -Mc2;
		double P3 = J1*OMEGA1*w1*Math.sin(delta1) + J2*OMEGA2*w2*Math.cos(delta2) - J3*OMEGA3*(w1*Math.cos(delta3) + w2*Math.sin(delta3))-Mc3;
		
		double Delta1Dot=  (P1*Q22*Q33-P1*Q23*Q32-P2*Q12*Q33+P2*Q13*Q32+P3*Q12*Q23-P3*Q13*Q22)/(Q11*Q22*Q33-Q11*Q23*Q32-Q21*Q12*Q33+Q21*Q13*Q32+Q31*Q12*Q23-Q31*Q13*Q22);

		double Delta2Dot=  -(P1*Q21*Q33-P1*Q23*Q31-P2*Q11*Q33+P2*Q13*Q31+P3*Q11*Q23-P3*Q13*Q21)/(Q11*Q22*Q33-Q11*Q23*Q32-Q21*Q12*Q33+Q21*Q13*Q32+Q31*Q12*Q23-Q31*Q13*Q22);

		double Delta3Dot=  (P1*Q21*Q32-P1*Q22*Q31-P2*Q11*Q32+P2*Q12*Q31+P3*Q11*Q22-P3*Q12*Q21)/(Q11*Q22*Q33-Q11*Q23*Q32-Q21*Q12*Q33+Q21*Q13*Q32+Q31*Q12*Q23-Q31*Q13*Q22);

	   //Relative Angular Momentum
		double Hb1_cmg= J2*OMEGA2*Math.cos(delta2)-J3*OMEGA3*Math.sin(delta3)    ;
   		double Hb2_cmg=  -J1*OMEGA1*Math.sin(delta1) + J3*OMEGA3*Math.cos(delta3) ;
   		double Hb3_cmg=  J1*OMEGA1*Math.cos(delta1) - J2*OMEGA2*Math.sin(delta2) ;

   
   		double Hb1_sc=I1*w1;
   		double Hb2_sc=I2*w2;
   		double Hb3_sc=I3*w3;
   
   		double Hb1=Hb1_cmg+Hb1_sc;
   		double Hb2=Hb2_cmg+Hb2_sc;
   		double Hb3=Hb3_cmg+Hb3_sc;
   
   		double Hb_total=Math.sqrt(Hb1*Hb1+Hb2*Hb2+Hb3*Hb3);
   		
   		double [] angular_momentum_array = new double[3];
			angular_momentum_array[0] = Hb1;
			angular_momentum_array[1] = Hb2;
			angular_momentum_array[2] = Hb3;
			Matrix angular_momentum_vector = new Matrix(angular_momentum_array, 3);
		Matrix inertial_angular_momentum = new Matrix(T_transpose.times(angular_momentum_vector));	
		double Hi1 = inertial_angular_momentum.get(0,0);
		double Hi2 = inertial_angular_momentum.get(1,0);
		double Hi3 = inertial_angular_momentum.get(2,0);
   
   		double Hi_total = Math.sqrt(Hi1*Hi1+Hi2*Hi2+Hi3*Hi3);
   		
  	 	//************* Rotational Kinetic Energy *******************
   		double total_velocity1=OMEGA1*Math.cos(delta1)-OMEGA2*Math.sin(delta2)+w1;
   		double total_velocity2=OMEGA2*Math.cos(delta2)+OMEGA3*Math.sin(delta3)+w2;
   		double total_velocity3=-OMEGA1*Math.sin(delta1)+OMEGA3*Math.cos(delta3)+w3;
   
   		double total_H1=Hb1_cmg+I1*w1;
   		double total_H2=Hb2_cmg+I2*w2;
   		double total_H3=Hb3_cmg+I3*w3;
   
   		double energy=0.5*(total_velocity1*total_H1 + total_velocity2*total_H2 + total_velocity3*total_H3);

        
        
       // add data point to the plot
        // FourPlots wheel_velocity_plot = new FourPlots();								
		// ThreePlots angullar_velocity_plot = new ThreePlots();
		// SinglePlot quarternion_check = new SinglePlot();
		// ThreePlots euler_angle_plot = new ThreePlots();
		// FourPlots angular_momentum_plot = new FourPlots();
		// SinglePlot energy_plot = new SinglePlot();
		
		gimbal_rate_plot.topPlot.addPoint(0, t, Delta1Dot, first);
		gimbal_rate_plot.middlePlot.addPoint(0, t, Delta2Dot, first);
		gimbal_rate_plot.bottomPlot.addPoint(0, t, Delta3Dot, first);
		
		gimbal_angle_plot.topPlot.addPoint(0, t, y[7], first);
		gimbal_angle_plot.middlePlot.addPoint(0, t, y[8], first); 
		gimbal_angle_plot.bottomPlot.addPoint(0, t, y[9], first);
		//wheel_velocity_plot.fourthPlot.addPoint(0, t, y[10], first);
		
		angular_velocity_plot.topPlot.addPoint(0, t,y[0], first);
		angular_velocity_plot.middlePlot.addPoint(0, t, y[1], first);
		angular_velocity_plot.bottomPlot.addPoint(0, t, y[2], first);
		
        quarternion_check.plot.addPoint(0, t, quat_check, first);
        
        euler_angle_plot.topPlot.addPoint(0, t, Psi, first);
		euler_angle_plot.middlePlot.addPoint(0, t, Theta, first);
		euler_angle_plot.bottomPlot.addPoint(0, t, Phi, first);  
		
		angular_momentum_plot.firstPlot.addPoint(0, t, Hi1, first);
		angular_momentum_plot.secondPlot.addPoint(0, t, Hi2, first);
		angular_momentum_plot.thirdPlot.addPoint(0, t, Hi3, first);
		angular_momentum_plot.fourthPlot.addPoint(0, t, Hi_total, first);
		
		energy_plot.plot.addPoint(0, t, energy, first);		
		
		quaternion_plot.firstPlot.addPoint(0, t, q1, first);
		quaternion_plot.secondPlot.addPoint(0, t, q2, first);
		quaternion_plot.thirdPlot.addPoint(0, t, q3, first);
		quaternion_plot.fourthPlot.addPoint(0, t, q4, first); 
       
        
        // Store quarternion values for use in animation
        quat_values[0][currentPts] = (float)t; // time value
        quat_values[1][currentPts] = (float)q1; // quaternion 1
        quat_values[2][currentPts] = (float)q2; // quarternion 2
        quat_values[3][currentPts] = (float)q3; // quarternion 3
        quat_values[4][currentPts] = (float)q4; // quarternion 4
        
    }
	
	/**
     * This method converts euler angle representation of attitude to that of
     * quarternions
     * @param a		psi angle
     * @param b		theta angle
     * @param c		phi angle
     * @return		A vector containing quaternions 
     */
    double[] convertEulerToTargetQ(double a, double b, double c)
    {
    	double psi = a;
    	double theta = b;
    	double phi = c;
    	
    	//P.102 321 rotation Transformation matrix
		double a11=Math.cos(theta)*Math.cos(psi);
		double a12=Math.cos(theta)*Math.sin(psi);
		double a13=-Math.sin(theta);
		double a21=-Math.cos(phi)*Math.sin(psi)+Math.sin(phi)*Math.sin(theta)*Math.cos(psi);
		double a22=Math.cos(phi)*Math.cos(psi)+Math.sin(phi)*Math.sin(theta)*Math.sin(psi);
		double a23=Math.sin(phi)*Math.cos(theta);
		double a31=Math.sin(phi)*Math.sin(psi)+Math.cos(phi)*Math.sin(theta)*Math.cos(psi);
		double a32=-Math.sin(phi)*Math.cos(psi)+Math.cos(phi)*Math.sin(theta)*Math.sin(psi);
		double a33=Math.cos(phi)*Math.cos(theta);

		//P.323- Definitions of Quaternion elements
		double trace=a11+a22+a33;
		double cos_alpha=(trace-1)/2;
		double alpha=Math.acos(cos_alpha);

		double e1=(a23-a32)/(2*Math.sin(alpha));
		double e2=(a31-a13)/(2*Math.sin(alpha));
		double e3=(a12-a21)/(2*Math.sin(alpha));

		double [] out = new double[4];
		
		out[0] = e1*Math.sin(alpha/2);
		out[1] = e2*Math.sin(alpha/2);
		out[2] = e3*Math.sin(alpha/2);
		out[3] =  Math.cos(alpha/2);
    
    	return out;
    	
    }// End of convertEulerToTargetQ
    
    /**
    * Return the quarternion values after simulation
    * @author	Noriko Takada
    */
   public float[][] getQuaternion()
   {
   		return quat_values;
   }
    	
   /**
    * Make the plots visible after simulation
    * @author	Noriko Takada
    */
   public void makePlotsVisible()
   {
   		// FourPlots wheel_velocity_plot = new FourPlots();								
		// ThreePlots angullar_velocity_plot = new ThreePlots();
		// SinglePlot quarternion_check = new SinglePlot();
		// ThreePlots euler_angle_plot = new ThreePlots();
		// FourPlots angular_momentum_plot = new FourPlots();
		// SinglePlot energy_plot = new SinglePlot();
   		//gimbal_rate_plot.setVisible(true);
        angular_velocity_plot.setVisible(true);
        quarternion_check.setVisible(true);
        euler_angle_plot.setVisible(true);
        angular_momentum_plot.setVisible(true);
        //energy_plot.setVisible(true);
        quaternion_plot.setVisible(true);
   }
    
    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args)
    {

        double time_step=0.1;
        double timeDuration=10;
        double tf = 20;
        double t0 = 0.0;
               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        
        // create an instance
        CMGManeuver si = new CMGManeuver(time_step, quat_values);
		
        // initialize the variables
        double [] x0 = new double[10];
        x0[0] = 0.0;
        x0[1] = 0.0;
        x0[2] = 0.0;
        x0[3] = 0.0;
        x0[4] = 0.0;
        x0[5] = 0.0;
        x0[6] = 1.0;
        x0[7] = 0.0;
        x0[8] = 0.0;
        x0[9] = 0.0;
       

        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        
        // make the plot visible
        si.makePlotsVisible();
    }
}
