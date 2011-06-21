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
import jat.matvec.data.Matrix;
import jat.plot.*;


/**
 * <p>
 * Equations of motion for One-axis Bang-Bang control with a fixed-level 
 * control torque
 * 
 * @author Noriko Takada
 * @version 1.3 (8/15/2004)
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion
 *   
 */
public class BangBangTwoD implements EquationsOfMotion
{
	//		time_step:			Time step of numerical integration
 	//		quat_values[][]		Two dimensional array that contains 
 	//							quarternions from simulation
 	//		J					Moment of inertia; Control is applied about 
 	//							the axis with this inertia
 	//		wn					Desired Control Natural Frequency
 	//		damping				Desired Damping Coefficient
 	//		K					Proportional Gain (wn^2*J), calculated 
 	//							assuming the control torque is available with 
 	//							no physical constraint
 	//		Kd					Derivative Gain (J*2*damping*wn), calculated 
 	//							assming the control
 	//                          torque is available with no physical 
 	//							constraint
 	//		K_bang				Proportional Gain for Bang-Bang control
 	//		Kd_bang				Derivative Gain for Bang-Bang control
 	//		torque_level		Fixed level of torque for Bang-Bang control
 	//		no_torque_zone		Amount of dead zone in which no control torque 
 	//							is applied
 	//		theta				Current position
 	//		theta_com			Commanded position
 	//		thetaPlot       	Plots of angular position
 	//		torquePlot		    Plots of calculated torque
 	//		quat_check      	Plot of e1^2 + e2^2 + e3^2 +e4^2
 	//		quaternionPlot		Plots of quaternions
 	//  
	double  time_step;
	private float quat_values[][];
	
	private double J = 600;
    private double wn = 1;
    private double damping = 1;
    private double K = wn*wn*J;
    private double Kd =J*2*damping*wn;
    private double K_bang = K;
    private double Kd_bang = Kd;
    private double torque_level=1;
    private double no_torque_zone=0.001;
    private double theta_com = 2*Math.PI/180;    
    
    // Create variables for the necessary plots
	private TwoPlots thetaPlot = new TwoPlots();
	private TwoPlots thetaDotPlot = new TwoPlots();
	private TwoPlots torquePlot = new TwoPlots();
	private TwoPlots phasePlot = new TwoPlots();
	private FourPlots quaternionPlot = new FourPlots();
	private SinglePlot quat_check = new SinglePlot();							
	
	
	/**
	 * Constructor
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 */
	public BangBangTwoD(double time_step, float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
	}
	
	
	/** Construct a RConstantTorque object
	 * @param		time_step:			Time step of numerical integration
 	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
 	 * @param		J					Moment of inertia about the axis that is to be controled
 	 * @param		wn					Desired Control Natural Frequency
 	 * @param		damping				Desired Damping Coefficient
 	 * @param		K					Proportional Gain (wn^2*J), calculated assuming 
 	 *                                  the control torque is available with no physical constraint
 	 * @param		Kd					Derivative Gain (J*2*damping*wn), calculated assming the control
 	 *                                  torque is available with no physical constraint
 	 * @param		K_bang				Proportional Gain for Bang-Bang control
 	 * @param		Kd_bang				Derivative Gain for Bang-Bang control
 	 * @param		torque				Fixed level of torque for Bang-Bang control
 	 * @param		dz					Amount of dead zone in which no control torque is applied
 	 * @param		theta_com			Commanded position
	 */
	public BangBangTwoD(double time_step, double J,double wn, 
						double damping,  double K, double Kd,  
						double K_bang, double Kd_bang,double torque, 
						double dz, double theta_com, float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.wn = wn;
		this.damping = damping;
		this.K = K;
		this.Kd = Kd;	
		this.K_bang=K_bang;
		this.Kd_bang = Kd_bang;
		this.torque_level = torque;
		this.no_torque_zone = dz;
		this.theta_com = theta_com*Math.PI/180;
		this.J=J;
	}
	
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		/*
		private TwoPlots thetaPlot = new TwoPlots();
		private TwoPlots thetaDotPlot = new TwoPlots();
		private TwoPlots torquePlot = new TwoPlots();
		private TwoPlots phasePlot = new TwoPlots();
		private FourPlots quaternionPlot = new FourPlots();
		private SinglePlot quat_check = new SinglePlot();			
		*/
		 
		// Setup plots
		thetaPlot.setTitle("Position Plot");
		thetaPlot.topPlot.setXLabel("t (sec)");
		thetaPlot.topPlot.setYLabel("[No Physical Constraint] theta (degree)");
		thetaPlot.bottomPlot.setXLabel("t (sec)");
		thetaPlot.bottomPlot.setYLabel("[Bang-Bang] theta (degree/sec)");
		
		thetaDotPlot.setTitle("Velocity Plot");
		thetaDotPlot.topPlot.setXLabel("t (sec)");
		thetaDotPlot.topPlot.setYLabel("[No Physical Constraint] thetaDot (deg/sec)");
		thetaDotPlot.bottomPlot.setXLabel("t (sec)");
		thetaDotPlot.bottomPlot.setYLabel("[Bang-Bang] thetaDot (deg/sec)");
		
		torquePlot.setTitle("Torque  Plot");
		torquePlot.topPlot.setXLabel("t (sec)");
		torquePlot.topPlot.setYLabel("[No Physical Constraint] Torque");
		torquePlot.bottomPlot.setXLabel("t (sec)");
		torquePlot.bottomPlot.setYLabel("[Bang-Bang] Torque");
		
		phasePlot.setTitle("phase Plot");
		phasePlot.topPlot.setXLabel("theta");
		phasePlot.topPlot.setYLabel("thetaDot");
		phasePlot.bottomPlot.setXLabel("theta");
		phasePlot.bottomPlot.setYLabel("thetaDot");
		
		quaternionPlot.setTitle("Quaternions");
		quaternionPlot.firstPlot.setXLabel("t (sec)");
		quaternionPlot.firstPlot.setYLabel("q1");
		quaternionPlot.secondPlot.setXLabel("t (sec)");
		quaternionPlot.secondPlot.setYLabel("q2");
		quaternionPlot.thirdPlot.setXLabel("t (sec)");
		quaternionPlot.thirdPlot.setYLabel("q3");
		quaternionPlot.fourthPlot.setXLabel("t (sec)");
		quaternionPlot.fourthPlot.setYLabel("q4");
		
		quat_check.setTitle("Quaternion Check");
		quat_check.plot.setXLabel("t (sec)");
		quat_check.plot.setYLabel("q1^2 + q2^2 + q3^2 + q4^2 = 1");
		
		
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
       	
       	//-------------- Equation of Motion for Non Bang-Bang control------
       	double [] out = new double[4];
       	// Define state variables
       	double theta			= x[0];
       	double thetaDot  		= x[1];
       	double	thetaBang 		= x[2];
       	double thetaDotBang	= x[3];
       	
       	double Tc=K*(theta_com-theta)-Kd*thetaDot;

		out[0]=thetaDot;
		out[1]=Tc/J;
       	//Tc=torque_level*sign(K_bang*(theta_com-theta)-(Kd_bang)*theta_dot);
		//if abs((theta_com-theta)*180/pi)<no_torque_zone
		//Tc=0;
		//end
		//%Tc=K*(theta_com-theta)-Kd*theta_dot;
		//ydot(1)=theta_dot;
		//ydot(2)=Tc/J;
		
		//--------- Signum Function -----------------
		//For each element of X, SIGN(X) returns 1 if the element
		//is greater than zero, 0 if it equals zero and -1 if it is
		//less than zero.  For complex X, SIGN(X) = X ./ ABS(X). 
		
		double TcBang=torque_level*sign(K_bang*(theta_com-thetaBang)-(Kd_bang)*thetaDotBang);
		if (Math.abs((theta_com-thetaBang)*180/Math.PI)<no_torque_zone)
			TcBang=0;
		
       	out[2]=thetaDotBang;
       	out[3]=TcBang/J;
       	
       	
       	
       	return out;
       	
     }// End of derivs
    	
    /** Implements the Printable interface to get the data out of the propagator 
     *  and pass it to the plot.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    	public void print(double t, double [] y)
    	{
        	boolean first = true;
        	if (t == 0.0) first = false;
        	
        	int currentPts = (int)(t/time_step); // This is the array index
        	
       	 	double theta = y[0];
        	double thetaDot = y[1];
        	double thetaBang = y[2];
        	double thetaDotBang = y[3];
        	double torqueBang =0;
        	if (Math.abs((theta_com-thetaBang)*180/Math.PI)<no_torque_zone)
      			torqueBang=0;
   			else 
      			torqueBang=torque_level*sign(K*(theta_com-thetaBang)-(Kd)*thetaDotBang);   
	        	
        	double torque = K*(theta_com-theta)-Kd*thetaDot;

        	// Extract quaternion values
        	// See Bong wie (p.318)
        	// Euler axis
        	double e1=0;
        	double e2=0;
        	double e3=1;
        	// Define Euler Parameters
        	double q1=e1*Math.sin(thetaBang/2);
        	double q2=e2*Math.sin(thetaBang/2);
        	double q3=e3*Math.sin(thetaBang/2);
        	double q4=Math.cos(thetaBang/2);
        	
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
        	double quat = q1*q1 + q2*q2 + q3*q3 + q4*q4;
        	
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
    	    		
        	/*
			private TwoPlots thetaPlot = new TwoPlots();
			private TwoPlots thetaDotPlot = new TwoPlots();
			private TwoPlots torquePlot = new TwoPlots();
			private TwoPlots phasePlot = new TwoPlots();
			private FourPlots quaternionPlot = new FourPlots();
			private SinglePlot quat_check = new SinglePlot();			
			*/
		        	        	       	
        	thetaPlot.topPlot.addPoint(0, t, (theta*180/Math.PI), first);
        	thetaPlot.bottomPlot.addPoint(0, t, (thetaBang*180/Math.PI), first);
        	
        	thetaDotPlot.topPlot.addPoint(0, t, (thetaDot*180/Math.PI), first);
        	thetaDotPlot.bottomPlot.addPoint(0, t, (thetaDotBang*180/Math.PI), first);
        	
        	torquePlot.topPlot.addPoint(0, t, torque, first);
        	torquePlot.bottomPlot.addPoint(0,t, torqueBang, first);
        	
        	phasePlot.topPlot.addPoint(0, theta, thetaDot, first);
        	phasePlot.bottomPlot.addPoint(0, thetaBang, thetaDotBang, first);
        	
        	quaternionPlot.firstPlot.addPoint(0, t, q1, first);
        	quaternionPlot.secondPlot.addPoint(0,t,q2, first);
        	quaternionPlot.thirdPlot.addPoint(0, t, q3, first);
        	quaternionPlot.fourthPlot.addPoint(0,t,q4, first);
        	
        	quat_check.plot.addPoint(0, t, quat, first);
        	
        	quat_values[0][currentPts] = (float)t; // time value
        	quat_values[1][currentPts] = (float)q1; // quaternion 1
        	quat_values[2][currentPts] = (float)q2; // quarternion 2
        	quat_values[3][currentPts] = (float)q3; // quarternion 3
        	quat_values[4][currentPts] = (float)q4; // quarternion 4
        	//also print to the screen 
        	System.out.println(t+" "+y[0]+" "+y[1]+" "+y[2]);
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
    * Make the plots visible after simulation
    * @author	Noriko Takada
    */
   public void makePlotsVisible()
   {
   		/*
		private TwoPlots thetaPlot = new TwoPlots();
		private TwoPlots thetaDotPlot = new TwoPlots();
		private TwoPlots torquePlot = new TwoPlots();
		private TwoPlots phasePlot = new TwoPlots();
		private FourPlots quaternionPlot = new FourPlots();
		private SinglePlot quat_check = new SinglePlot();			
			*/
		thetaPlot.setVisible(true);
		thetaDotPlot.setVisible(true);
		torquePlot.setVisible(true);
		phasePlot.setVisible(true);
		quaternionPlot.setVisible(true);
		quat_check.setVisible(true);
   }
   
   /** Runs the example.
    * @param args Arguments.
    */
    public static void main(String[] args)
    {

        double time_step=0.01;
       double timeDuration=30;
    
        double t0 = 0.0;
               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
        RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78(1e-6);
		double tf=timeDuration	;// Duration of the simulation time is the same as the final time
		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
		
		 
        // create an instance
        BangBangTwoD si = new BangBangTwoD(time_step,  quat_values);
		
        // initialize the variables
        double [] x0 = new double[4];
        x0[0] = 0.0;
        x0[1] = 0.0;
        x0[2] = 0.0;
        x0[3] = 0.0;
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        //double [] xf= rk78.integrate(t0, x0, tf, si, si, true);
        
        // make the plot visible
        si.makePlotsVisible();        
    }
    
    /**
     * Returnss 1 if x is greater than zero, 0 if it equals zero, and 
     * -1 if it is less than zero
     * @param	x		double value
     */
    int sign(double x)
    {
    	int value=0;
    	if(x>0)
    		value =1;
    	else if (x<=0)
    		value = -1;
    	else
    		value = 0;
    		
        return value;			
    }
    
}

    		