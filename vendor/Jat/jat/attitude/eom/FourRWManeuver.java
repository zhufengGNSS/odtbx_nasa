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
import jat.attitude.QuatToDeg;
import jat.matvec.data.*;
import jat.plot.*;

/**
 * <P>
 * FourRWManeuver contains the equations of motion for the simulation
 * of a rigid spacecraft performing attitude reorientation maneuver 
 * using 4 reaction wheels in a pyramid configuration.
 * The class implements jat.alg.integrators.Derivatives and
 * jat.alg.integrateors.Printble, so an outside application code can 
 * perform a numerical simulation of equations of motion defined 
 * in this class and get the output from the simulation.
 * 
 * @author Noriko Takada
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion 
 */

public class FourRWManeuver implements EquationsOfMotion
{
//time_step:		Time step of numerical integration
 	//quat_values[][]	Two dimensional array that contains quarternions from 
//simulation
 	//I1,I2,I3		Principal moments of inertia about 1,2, and 3 axes
 	//M1,M2,M3		Extermal moment applied about 1,2, and 3 axes
 	//J			Reaction-wheel inertia (all 4 RWs have the same J.)
 	//angle 		All RW is oriented at this angle with respect to the z-axis
 	//b1 b2 b3		Vector components of the moment vector of RW1	
 	//g1 g2 g3		Vector components of the moment vector of RW2
 	//a1 a2 a3		Vector components of the moment vector of RW3
	//d1 d2 d3		Vector componetns of the moment vector of RW4
 	//Kx Ky Kz		Proportional gain
 	//Kxd Kyd Kzd		Derivative gain
 	//psi, theta, phi	Amount of attitude maneuver in terms of euler angles
 	//q			Target quarternions (converted from the euler angles)
	double  time_step;
	private float quat_values[][];
	public static final int THREE_RW=3;
	public static final int FOUR_RW=4;
	int number_of_RW = THREE_RW;
	
    	// Create variables for the necessary plots
	private FourPlots wheel_velocity_plot = new FourPlots();					private ThreePlots angular_velocity_plot = new ThreePlots();
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
    
    	double J = 10; 
    	double angle =25*Math.PI/180; 
    							  
	double b1=-Math.sin(angle);
	double b2=0;
	double b3=Math.cos(angle);

	double g1=0;
	double g2=-Math.sin(angle);
	double g3=Math.cos(angle);

	double a1=Math.sin(angle);
	double a2=0;
	double a3=Math.cos(angle);

	double d1=0;
	double d2=Math.sin(angle);
	double d3=Math.cos(angle);
	
	// Gains (self-determined...)
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
	private double phi=50*Math.PI/180;
	
	// Target Quarternions
	double[] q = convertEulerToTargetQ(psi, theta, phi);//new double[4];
		
	// Convert Euler angles to target quarternions
	//double[] target = convertEulerToTargetQ(psi, theta, phi);
	
	
	/** 
	 *  Constructor 1
	 *  @param time_step:		Time step of numerical integration
    	 *  @param quat_values[][]	Two dimensional array that contains quarternions 
from simulation	
	 *  @param I1 I2 I3		Principle moments of inertia about 1, 2, and 3 axis 
respectively
	 *  @param J	Moment of inertia of reaction wheel
	 *  @param angle	Orientation angle of reaction wheel
	 *  @param psi	Target euler angle
	 *  @param phi	Target euler angle
	 *  @param theta	Target euler angle
	 */
	public FourRWManeuver(double time_step, double I1, double I2, double I3
			, double J, double angle, double a, double b, double c,
			float quat_values[][] )
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
		this.J = J;
		this.angle = angle*Math.PI/180;
		//this.psi = psi;
		//this.phi = phi;
		//this.theta = theta;
		// a,b, and c are angles in degrees.
		 psi=a*Math.PI/180;
		 theta=b*Math.PI/180;
		 phi=c*Math.PI/180;
		
		q = convertEulerToTargetQ(psi, theta, phi);
	}
	
	/**
	 * Constructor 2
	 * @param time_step:		Time step of numerical integration
    	 * @param quat_values[][]	Two dimensional array that contains 
quarternions from simulation
	 * @param a			Target psi angle
	 * @param b			Target theta angle
	 * @param c 			Target phi angle
	 */
	public FourRWManeuver(double time_step,double a, double b, double c, float 
quat_values[][])
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
	 * Constructor 3
	 * @param time_step:		Time step of numerical integration
    	 * @param quat_values[][]	Two dimensional array that contains quarternions 
from simulation
    	 */
	public FourRWManeuver(double time_step, float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
	}
	
	
	/**
	 * Sets the number of RW to either three or four 
	 * @param	a 	(int)	# of RW
	 */
	void setNumberOfRW( int a)
	{
		number_of_RW = a;	
	}
	
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		// FourPlots wheel_velocity_plot = new FourPlots();							// ThreePlots angullar_velocity_plot = new ThreePlots();
		// SinglePlot quarternion_check = new SinglePlot();
		// ThreePlots euler_angle_plot = new ThreePlots();
		// FourPlots angular_momentum_plot = new FourPlots();
		// SinglePlot energy_plot = new SinglePlot();
		
		// Setup plots
		wheel_velocity_plot.setTitle("Reaction Wheel Velocities");
		wheel_velocity_plot.firstPlot.setXLabel("t (sec)");
		wheel_velocity_plot.firstPlot.setYLabel("RW 1 speed (rad/sec)");
		wheel_velocity_plot.secondPlot.setXLabel("t (sec)");
		wheel_velocity_plot.secondPlot.setYLabel("RW 2 speed (rad/sec)");
		wheel_velocity_plot.thirdPlot.setXLabel("t (sec)");
		wheel_velocity_plot.thirdPlot.setYLabel("RW 3 speed (rad/sec)");
		wheel_velocity_plot.fourthPlot.setXLabel("t (sec)");
		wheel_velocity_plot.fourthPlot.setYLabel("RW 4 speed (rad/sec)");
		
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
		
         		double [] out = new double[11];
       		out[0] = ((I2-I3)/I1)*x[1]*x[2] + M1/I1 + Mc1/I1;
       		out[1] = ((I3-I1)/I2)*x[0]*x[2] + M2/I2 + Mc2/I2;
       		out[2] = ((I1-I2)/I3)*x[0]*x[1] + M3/I3 + Mc3/I3;
       		out[3] = -0.5* (-x[2]*x[4] + x[1]*x[5] - x[0]*x[6]);
       		out[4] = -0.5* (x[2]*x[3]  - x[0]*x[5] - x[1]*x[6]);
       		out[5] = -0.5* (-x[2]*x[6] + x[0]*x[4] - x[1]*x[3]);
       		out[6] = -0.5* (x[2]*x[5]  + x[1]*x[4] + x[0]*x[3]);
             
       		//RW angular speed
       		if (number_of_RW ==FOUR_RW)
       		{
       			// RWs speed 
       			// x[7] = RW1 spped
       			// x[8] = RW2 speed
       			// x[9] = RW3 speed
       			// x[10] = RW4 speed
       	
double p1 = J*x[7]*(b3*x[1]-b2*x[2])+J*x[8]*(g3*x[1]-
g2*x[2])+J*x[9]*(a3*x[1]-a2*x[2])+J*x[10]*(d3*x[1]-d2*x[2]);
double p2 = J*x[7]*(b1*x[2]-b3*x[0])+J*x[8]*(g1*x[2]-
g3*x[0])+J*x[9]*(a1*x[2]-a3*x[0])+J*x[10]*(d1*x[2]-d3*x[0]);
double p3 = J*x[7]*(b2*x[0]-b1*x[1])+J*x[8]*(g2*x[0]-
g1*x[1])+J*x[9]*(a2*x[0]-a1*x[1])+J*x[10]*(d2*x[0]-d1*x[1]);

	        		// P matrix
    	    		double [] p_array = new double[3];
        		p_array[0] = p1;
        		p_array[1] = p2;
        		p_array[2] = p3;
        		Matrix p_matrix = new Matrix(p_array, 3); // Column vector
        	
        		// Q matrix
        		Matrix pseudoQ_matrix = takePseudoInverse( J,  b1,  b2, b3
    							 , g1,  g2,  g3,
    							   a1,  a2,  a3,
    							   d1,  d2,  d3);
        
        		//double[][] q_array = new double[3][4];
    			// Row1
    			//q_array[0][0] = J*b1;
    			//q_array[0][1] = J*g1;
    			//q_array[0][2] = J*a1;
    			//q_array[0][3] = J*d1;
    			// Row2
    			//q_array[1][0] = J*b2;
    			//q_array[1][1] = J*g2;
    			//q_array[1][2] = J*a2;
    			//q_array[1][3] = J*d2;
    			// Row3
    			//q_array[2][0] = J*b3;
    			//q_array[2][1] = J*g3;
    			//q_array[2][2] = J*a3;
    			//q_array[2][3] = J*d3;	    	
    			//Construct Transformation Matrix
    			//Matrix  Q = new Matrix(q_array,3,4);
    			//Matrix 	pseudoQ = Q.inverse();
    	
    			double [] Tc_array = new double[3];
    			Tc_array[0] = -Mc1;
    			Tc_array[1] = -Mc2;
    			Tc_array[2] = -Mc3;
    			Matrix Tc_matrix = new Matrix(Tc_array,3);
    		
    			Matrix omega_matrix = 
pseudoQ_matrix.times(Tc_matrix.minus(p_matrix));
    	
    			//Matrix  Q_transpose = new Matrix(Q.transpose().A, 4,3);
    	        	// Pseudoinverse by manual calculation....
                       // matrix element
       			double q11 = -0.5/(J*Math.sin(angle));
       			double q13 = 0.25/(J*Math.cos(angle));
       			double q22 = -0.5/(J*Math.sin(angle));
       			double q23 = 0.25/(J*Math.cos(angle));
       			double q31 = 0.5/(J*Math.sin(angle));
       			double q33 = 0.25/(J*Math.cos(angle));
       			double q42 = 0.5/(J*Math.sin(angle));
       			double q43 = 0.25/(J*Math.cos(angle));
       	
       			// matrix element
       			double a11 = -Mc1-p1;
       			double a12 = -Mc2-p2;
       			double a13 = -Mc3-p3;
       	
       	       		out[7] = q11*a11 + q13*a13;
       			out[8] = q22*a12 + q23*a13;
       			out[9] = q31*a11 + q33*a13;
       			out[10] = q42*a12 + q43*a13;
       		}
       	
       		else if (number_of_RW == THREE_RW)
       		{       			
       			Matrix Q = new Matrix(3,3);
       			Matrix P = new Matrix(3,1);
       			Matrix Tc = new Matrix(3,1);
       			Matrix Q_inv = new Matrix(3,3);
       			Matrix Omega = new Matrix(3,1);
       	
       			// Set the matrix elements
       			// Don't forget that the index starts at 0!!
       			Q.set(0,0,J*b1);
       			Q.set(0,1, J*g1);
       			Q.set(0,2, J*a1);
       			Q.set(1,0, J*b2);
       			Q.set(1,1,J*g2);
       			Q.set(1,2, J*a2);
       			Q.set(2,0, J*b3);
       			Q.set(2,1, J*g3);
       			Q.set(2,2,J*a3);
       	
			P.set(0,0, J*x[7]*(b3*x[1]-b2*x[2])+J*x[8]*(g3*x[1]-
g2*x[2])+J*x[9]*(a3*x[1]-a2*x[2]) );
			P.set(1,0, J*x[7]*(b1*x[2]-b3*x[0])+J*x[8]*(g1*x[2]-
g3*x[0])+J*x[9]*(a1*x[2]-a3*x[0]));
			P.set(2,0, J*x[7]*(b2*x[0]-b1*x[1])+J*x[8]*(g2*x[0]-
g1*x[1])+J*x[9]*(a2*x[0]-a1*x[1]));
       	
       			Tc.set(0,0, -Mc1);
       			Tc.set(1,0, -Mc2);
       			Tc.set(2,0, -Mc3);
       		
       			Q_inv=Q.inverse();
       			Omega=Q_inv.times((Tc.minus(P)));
       	
       			out[7] = Omega.get(0,0);
       			out[8] = Omega.get(1,0);
       			out[9] = Omega.get(2,0);
       			out[10] = 0;
       		}
       		return out;
    	}// End of derivs
    
     	/** Implements the Printable interface to get the data out of the propagator and 
pass it to the plot.
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
        	double Omega1 = y[7];
        	double Omega2 = y[8];
        	double Omega3 = y[9];
        	double Omega4 = y[10];
        
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
    		//Matrix  q_times_qtrans = new Matrix(q_matrix.times(q_transpose).A, 3,3);
    		//Matrix  inv_q_times_qtrans = new Matrix(q_times_qtrans.invert().A, 3,3);
    		//double[][] pseudoinv_q = q_transpose.times(inv_q_times_qtrans).A;
        
         	double Hb1_wheel=J*(Omega1*b1+Omega2*g1+Omega3*a1+Omega4*d1);
   		double Hb2_wheel=J*(Omega1*b2+Omega2*g2+Omega3*a2+Omega4*d2);
   		double Hb3_wheel=J*(Omega1*b3+Omega2*g3+Omega3*a3+Omega4*d3);
   		double Wheelb1=(Omega1*b1+Omega2*g1+Omega3*a1+Omega4*d1);
   		double Wheelb2=(Omega1*b2+Omega2*g2+Omega3*a2+Omega4*d2);
        double Wheelb3=(Omega1*b3+Omega2*g3+Omega3*a3+Omega4*d3);
		double energy=0.5*((Wheelb1+w1)*(J*Wheelb1+I1*w1)+(Wheelb2+w2)*(J*Wheelb2+I2*w2)+(Wheelb3+w3)*(J*Wheelb3+I3*w3));
   
   		double Hb1_sc=I1*w1;
   		double Hb2_sc=I2*w2;
   		double Hb3_sc=I3*w3;
   
   		double Hb1=Hb1_wheel+Hb1_sc;
   		double Hb2=Hb2_wheel+Hb2_sc;
   		double Hb3=Hb3_wheel+Hb3_sc;
        
        	//Construct Angular velocity vector
    		//double [] angularRateArray = new double[3];
    		//angularRateArray[0] = wx;
    		//angularRateArray[1] = wy;
    		//angularRateArray[2] = wz;
    		//Matrix angularRateVector = new Matrix(angularRateArray, 3); 
//Column vector
    		//Transform angular rates in body frame to inertial frame
    		//Matrix inertialAngularRateVector = new 
//Matrix(T_transpose.times(angularRateVector));
    		//double wx_i = inertialAngularRateVector.get(0,0);
    		//double wy_i = inertialAngularRateVector.get(1,0);
    		//double wz_i = inertialAngularRateVector.get(2,0);
        
		double [] angular_momentum_array = new double[3];
		angular_momentum_array[0] = Hb1;
		angular_momentum_array[1] = Hb2;
		angular_momentum_array[2] = Hb3;
		Matrix angular_momentum_vector = new Matrix(angular_momentum_array, 3);
		Matrix inertial_angular_momentum = new 
Matrix(T_transpose.times(angular_momentum_vector));	
		double Hi1 = inertial_angular_momentum.get(0,0);
		double Hi2 = inertial_angular_momentum.get(1,0);
		double Hi3 = inertial_angular_momentum.get(2,0);
		double Hi_total = Math.sqrt(Hi1*Hi1+Hi2*Hi2+Hi3*Hi3);
                
        	// add data point to the plot
        	// FourPlots wheel_velocity_plot = new FourPlots();								
		// ThreePlots angullar_velocity_plot = new ThreePlots();
		// SinglePlot quarternion_check = new SinglePlot();
		// ThreePlots euler_angle_plot = new ThreePlots();
		// FourPlots angular_momentum_plot = new FourPlots();
		// SinglePlot energy_plot = new SinglePlot();
		
		wheel_velocity_plot.firstPlot.addPoint(0, t, y[7], first);
		wheel_velocity_plot.secondPlot.addPoint(0, t, y[8], first); 
		wheel_velocity_plot.thirdPlot.addPoint(0, t, y[9], first);
		wheel_velocity_plot.fourthPlot.addPoint(0, t, y[10], first);
		
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
    	 * @param a	psi angle
    	 * @param b	theta angle
    	 * @param c	phi angle
    	 * @return	A vector containing quaternions 
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
		double a21=-
Math.cos(phi)*Math.sin(psi)+Math.sin(phi)*Math.sin(theta)*Math.cos(psi);
double 
a22=Math.cos(phi)*Math.cos(psi)+Math.sin(phi)*Math.sin(theta)*Math.sin(psi);
		double a23=Math.sin(phi)*Math.cos(theta);
		double 
a31=Math.sin(phi)*Math.sin(psi)+Math.cos(phi)*Math.sin(theta)*Math.cos(psi);
		double a32=-
Math.sin(phi)*Math.cos(psi)+Math.cos(phi)*Math.sin(theta)*Math.sin(psi);
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
    	 * This mehtod takes pseudo inverse given the appropriate inputs
    	 * @param inertia
    	 * @param w1 w2 w3
    	 * @param x1 x2 x3
    	 * @param y1 y2 y3
    	 * @param z1 z2 z3
    	 * @return
    	 */
    	Matrix takePseudoInverse(double inertia, double w1, double w2,double w3
    							 ,double x1, double x2, double x3,
    							  double y1, double y2, double y3,
    							  double z1, double z2, double z3)
    	{
    		// Local variable
    		double J = inertia;
    		double b1 = w1;
    		double b2 = w2;
    		double b3 = w3;
    		double g1 = x1;
    		double g2 = x2;
    		double g3 = x3;
    		double a1 = y1;
    		double a2 = y2;
    		double a3 = y3;
    		double d1 = z1;
    		double d2 = z2;
    		double d3 = z3;
    	
    		// Build 3 x 4 Q matrix
    		double[][] array = new double[3][4];
    		// Row1
    		array[0][0] = J*b1;
    		array[0][1] = J*g1;
    		array[0][2] = J*a1;
    		array[0][3] = J*d1;
    		// Row2
    		array[1][0] = J*b1;
    		array[1][1] = J*g2;
    		array[1][2] = J*a2;
    		array[1][3] = J*d2;
    		// Row3
    		array[2][0] = J*b3;
    		array[2][1] = J*g3;
    		array[2][2] = J*a3;
    		array[2][3] = J*d3;
    	
    		//Construct a Matrix object
    		Matrix  q_matrix = new Matrix(array,3,4);
    		Matrix  q_transpose = new Matrix(q_matrix.transpose().A, 4,3);
    		Matrix  q_times_qtrans = new Matrix(q_matrix.times(q_transpose).A, 3,3);
    		Matrix  inv_q_times_qtrans = new Matrix(q_times_qtrans.invert().A, 3,3);
    		double[][] pseudoinv_q = q_transpose.times(inv_q_times_qtrans).A;
    		Matrix pseudoinv_matrix = q_transpose.times(inv_q_times_qtrans);
    	
    		return pseudoinv_matrix;
    	
    	}
    
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
   		wheel_velocity_plot.setVisible(true);
        	angular_velocity_plot.setVisible(true);
        	//quarternion_check.setVisible(true);
        	euler_angle_plot.setVisible(true);
        	angular_momentum_plot.setVisible(true);
        	energy_plot.setVisible(true);
        	quaternion_plot.setVisible(true);
   	}
   
   	/** Runs the example.
    	 * @param args Arguments.
    	 */
    	public static void main(String[] args)
    	{
        	double time_step=0.01;
       	double timeDuration=10;
        	double tf = 20;
        	double t0 = 0.0;
              	double I1 = 10.42;
		double I2 = 35.42;
		double I3 = 41.67;
        	double J = 10;
    		double angle = 25;
    		double psi =10;
    		double theta = 60;
    		double phi = 50;
    		// initialize the variables
               RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	
// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				       	float quat_values[][]= new  float[5][numberOfPts+1];
// +1 is for AnimationWindow
              	// create an instance
        	//FourRWManeuver si = new FourRWManeuver(time_step, quat_values);
		FourRWManeuver si = new FourRWManeuver(time_step,I1,I2,I3, J, angle, psi, 
theta,phi,quat_values );
        	// initialize the variables
        	double [] x0 = new double[11];
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
        	x0[10] = 0.0;
       		// integrate the equations
        	rk8.integrate(t0, x0, tf, si, true);
        	// make the plot visible
        	si.makePlotsVisible();
    	}
    							  	 
    		
}// End of File
