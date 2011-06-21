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
import jat.alg.integrators.RungeKuttaFehlberg78;
import jat.attitude.util.AnimationWindow;
import jat.attitude.util.AnimationWindow2;


/**
 * EonTest class demonstrate how to use eom classes from outside
 * @author 	Noriko Takada
 * @version	1.6 (8/15/2004)
 * 
 * 		scenario = 1	->	Constant torque
 * 		scenario = 2	->	Gravity gradient (circular orbit)
 * 		scenario = 3	->	Gravity gradient (eccentric orbit)
 * 		scenario = 4	->	Rigid spacecraft with a spherical damper
 * 		scenario = 5	->	Four reaction wheel maneuver
 * 		scenario = 6	->	Three single-gimbal control moment gyro maneuver
 * 		scenario = 7	->	Bang-bang control (Under development)
 * 		scenario = 8	->	Simple 2D flexible spacecraft (Under development)
 * 		scenario = 9	->	Simple 3D flexible spacecraft (Under development)
 * 
 * Modification since the last version
 * 		--> Removed: import.jat.attitude	
 */
public class EomTest 
{
	// create an RungeKutta8 integrator with step-size of 0.1
    RungeKutta8 rk8 = new RungeKutta8(0.1);
    RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78(1e-6);
	static EomTest EOM;
	private double time_step = 0.1;
	private double timeDuration = 10;
	
	double i = 10.42;
	double j = 35.42;
	double k = 41.67;
	private double I1 = i;//45.42;
	private double I2 = j;//35.42;
	private double I3 = k;//13.35;
	
	private double J = 600;   // Moment of Inertia for a single axis simulation
	public int animationYes = 1;
	public int plotYes = 1;
	public int number_of_RW = 3;
	public float quat_values[][]= null;
	
	/**
	 * Default Constructor
	 */
	public EomTest(double time_step, double time)
	{
		this.time_step = time_step;
		this.timeDuration = time;
	}
	
	
	/**
	 * Constructor
	 */
	public EomTest(double time_step, double time, double I1, double I2, double I3)
	{
		this.time_step = time_step;
		this.timeDuration = time;
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
	}
	
	/**
	 * Sets up the private variables from outside
	 */
	public void setInertiaThreeD(double I1, double I2, double I3)
	{
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
	}
	
	/** Set moment of inertia for two dimensional simulation
	 */
	public void setInertiaTwoD(double J)
	{
		this.J = J;
	}
	
	public static void main(String[] args)
    {
		 EOM = new EomTest(0.01, 20.0);
		 //EOM = new EomTest(0.1, 10.0, 600);
		 int scenario = 5;
		 
		 if (scenario ==1)
		 {
		 	double	M1 = 1;           // External torque
    		double M2 = 0;
    		double M3 = 0; 
            // initialize the variables
        	double [] x0 = new double[7];
        	x0[0] = 0.0;
        	x0[1] = 0.1;
        	x0[2] = 0.0;
        	x0[3] = 0.0;
        	x0[4] = 0.0;
        	x0[5] = 0.0;
        	x0[6] = 1.0;    
		 	EOM.doConstantTorque(x0, M1, M2, M3);
		 }
		 else if (scenario ==2)
		 {
		 	// initialize the variables
        	double [] x0 = new double[7];
        	x0[0] = 0.1;
        	x0[1] = 0.0;
        	x0[2] = 1.0;
        	x0[3] = 0.0;
        	x0[4] = 0.0;
        	x0[5] = 0.0;
        	x0[6] = 1.0;
		 	EOM.doGGCircular(x0);
		 }
		 else if (scenario == 3)
		 {
		 	double e  = 0.3;
		 	// initialize the variables
        	double [] x0 = new double[7];
        	x0[0] = 0.0;
        	x0[1] = 0.0;
        	x0[2] = 1.0;
        	x0[3] = 0.0;
        	x0[4] = 0.0;
        	x0[5] = 0.0;
        	x0[6] = 1.0;
		 	EOM.doGGEccentric(x0, e);
		 }
		 else if (scenario == 4)
		 {
		 	double c  = 5;        // damping coefficient
    		double j  = 5; 
    		// Initial Condition
    		double [] x0 = new double[10];
        	x0[0] = 0.0;
        	x0[1] = 0.0;
        	x0[2] = 1.0;
        	x0[3] = 0.0;
        	x0[4] = 0.0;
        	x0[5] = 0.0;
        	x0[6] = 1.0;
        	x0[7] = 0.0;
        	x0[8] = 0.0;
        	x0[9] = 0.0;
		 	EOM.doSphericalDamper(x0,c,j);
		 }
		 else if (scenario == 5)
		 {
		 	double J = 10;
    		double angle = 25;
    		double psi =10;
    		double theta = 60;
    		double phi = 50;
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
		 	EOM.doFourRW(x0, J, angle, psi, theta, phi);
		 }
		 else if (scenario == 6)
		 {
		 	double J = 10;
    		double A = 500;
    		double psi =30;
    		double phi = 20;
    		double theta = 50;
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
		 	EOM.doCMGManeuver(x0, J, A, psi, phi, theta);
		 }
		 else if (scenario ==7)
		 {
		 	
    		double J = 600;
    		double wn = 1;
    		double damping = 1;
    		double K = wn*wn*J;
    		double Kd =J*2*damping*wn;
    		double K_bang = K;
    		double Kd_bang = Kd;
    		double torque_level=1;
    		double no_torque_zone=0.001;
    		double theta_com = 2*Math.PI/180;
    		// initialize the variables
        	double [] x0 = new double[4];
        	x0[0] = 0.0;
        	x0[1] = 0.0;
        	x0[2] = 0.0;
        	x0[3] = 0.0;
        	EOM.setInertiaTwoD(J);        	
		 	EOM.doBangBang(x0,  wn,  damping, K, Kd, K_bang, Kd_bang, 
		 					torque_level,  no_torque_zone, theta_com  );
		 }
		 else if (scenario == 8)
		 {
		 	//double J = 13.5;
    		double a = 0.40;
    		double L =1.8;
    		double EI = 15.52;
    		double m = 1;
    		double K = 0;
    		double Kd = 0;
    		double alpha_com = 0;
    		// initialize the variables
        	double [] x0 = new double[6];
        	x0[0] = 0.0;
        	x0[1] = 0.1;
        	x0[2] = 0.1;
        	x0[3] = 0.0;
        	x0[4] = 0.0;
        	x0[5] = 0.0;
        	
		 	EOM.doTwoDFlex(x0, a,L,EI,m, K, Kd, alpha_com);
		 }
		 else if (scenario == 9)
		 {
		 	double a = 0.40;
    		double L =1.8;
    		double EI = 15.52;
    		double m = 1;
    		// initialize the variables
    		double [] x0 = new double[10];
    		x0[0] = 0.1;
    		x0[1] = 0.1;
    		x0[2] = 0.0;
    		x0[3] = 0.0;
    		x0[4] = 0.0;
    		x0[5] = 0.0;
    		x0[6] = 0.0;
    		x0[7] = 0.0;
    		x0[8] = 0.0;
    		x0[9] = 0.0;
    		
		 	EOM.doThreeDFlex(x0, a,L,EI,m); 
		 }	
		 		 	
	}

	/**
	 * Method doThreeDFlex.
	 */
	public void doThreeDFlex(double [] x0,
							  double a, double L, double EI, double m) 
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
       	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
       	float quatBeam1[][]= new  float[5][numberOfPts+1];
		float quatBeam2[][]= new  float[5][numberOfPts+1];
        // create an instance
        FlexibleThreeD si = new FlexibleThreeD(time_step, m, a, L, EI, I1,I2,I3,quat_values, quatBeam1, quatBeam2 );
        
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if (plotYes ==1 )
        si.makePlotsVisible();
        
        // Animation
        quat_values = si.getQuaternion();
        if(animationYes ==1)
        {
			/*
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
														 (float)I1, (float)I2, (float)I3, 
														 numberOfPts, quat_values , "else");
			*/
			AnimationWindow2 theAnimWindow = new AnimationWindow2("Animation", 
															 (float)I3, (float)I3, (float)I3, 
															 numberOfPts, quat_values , "else",
															 quatBeam1, quatBeam2,(float)a, (float)L);
        }	
	}


	/**
	 * Method doTwoDFlex.
	 */
	public void doTwoDFlex(double [] x0,
							double a, double L, double EI, double m, double K, double Kd, double alpha_com) 
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
        RungeKuttaFehlberg78 rk78 = new RungeKuttaFehlberg78(1e-12);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
       	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        float quatBeam1[][]= new  float[5][numberOfPts+1];
		float quatBeam2[][]= new  float[5][numberOfPts+1];
        // create an instance
        FlexibleTwoD si = new FlexibleTwoD(time_step, m, a, L, EI, I3, K, Kd, alpha_com, quat_values, quatBeam1, quatBeam2 );
        
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
       // double [] xf = rk78f.integrate(x0, t0, tf, orbit, lp, true);
		//rk78.integrate(t0, x0, tf, si, si, true);
        // make the plot visible
        if (plotYes ==1 )
        	si.makePlotsVisible();
        
        // Animation
        quat_values = si.getQuaternion();
        quatBeam1 = si.getQuatBeam1();
        quatBeam2 = si.getQuatBeam2();
        if(animationYes ==1)
        {
			
			AnimationWindow2 theAnimWindow = new AnimationWindow2("Animation", 
															 (float)I3, (float)I3, (float)I3, 
															 numberOfPts, quat_values , "else",
															 quatBeam1, quatBeam2,(float)a, (float)L);
															 
			/*
			AnimationWithAppendages theAnimWindow = new AnimationWithAppendages("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "else",
															 quatBeam1, quatBeam2);*/
        }	
	}


	/**
	 * Method doBangBang.
	 */
	public void doBangBang(double [] x0,
							 double wn, double damping, double K, double Kd,
							double K_bang, double Kd_bang, double torque, double dz,
							double theta_com ) 
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
        //BangBangTwoD si = new BangBangTwoD(time_step,  quat_values);
		BangBangTwoD si = new BangBangTwoD(time_step,  J, wn,  damping, 
						  K,  Kd,   K_bang,  Kd_bang, 
						  torque,  dz,  theta_com,quat_values);
        // initialize the variables
        //double [] x0 = new double[4];
        //x0[0] = 0.0;
        //x0[1] = 0.0;
        //x0[2] = 0.0;
        //x0[3] = 0.0;
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        //double [] xf= rk78.integrate(t0, x0, tf, si, si, true);
        
        // make the plot visible
        if (plotYes ==1 )
       		 si.makePlotsVisible();        
        I3 =J;
        I2 =J; // Just for convenience
        I1 =J; // Just for convenience
        if(animationYes ==1)
        {
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "else");
        }	
	}


	/**
	 * Method doFourRW.
	 */
	public void doFourRW(double [] x0, 
						  double J, double angle, double psi, double theta, double phi) 
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
       	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        // create an instance
        FourRWManeuver si = new FourRWManeuver(time_step, I1, I2, I3, J, angle, psi, theta,phi,quat_values );
        
        // Set the number of RW
        si.number_of_RW = number_of_RW;
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if (plotYes ==1 )
        	si.makePlotsVisible();
        
        // Animation
        quat_values = si.getQuaternion();
        if(animationYes ==1)
        {
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "else");
        }		
	}


	/**
	 * Method doSphericalDamper.
	 */
	public void doSphericalDamper(double[] x0,
									double c, double j) 
	{
		double tf = timeDuration;
        double t0 = 0.0;
               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;	
       	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        // create an instance
        RSphericalDamper si = new RSphericalDamper(time_step,I1,I2,I3,c,j, quat_values);
		
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if(plotYes == 1)
           si.makePlotsVisible();
        // Animation
        quat_values = si.getQuaternion();
        if (animationYes == 1)
        {
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
        													 numberOfPts, quat_values , "Else");
        }
	}


	/**
	 * Method doGGEccentric.
	 */
	public void doGGEccentric(double[] x0,
								double e) 
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        // create an instance
        RGGEccentricOrbit si = new RGGEccentricOrbit(time_step,I1,I2,I3,e, quat_values);
		
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if(plotYes ==1)
        	si.makePlotsVisible();
        // Animation
        quat_values = si.getQuaternion();
		
		if (animationYes == 1)
		{
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "Gravity Gradient");
		}
	}


	/**
	 * Method doGGCircular.
	 */
	public void doGGCircular(double[] x0)
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
		float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
        // create an instance
        RGGCircularOrbit si = new RGGCircularOrbit(time_step,I1,I2,I3, quat_values);
		
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if (plotYes == 1)
        	si.makePlotsVisible();
         // Animation
        quat_values = si.getQuaternion();
		if (animationYes ==1)
		{
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "Gravity Gradient");
		}
	}


	/**
	 * Method doConstantTorque.
	 */
	public void doConstantTorque(double[] x0,
							      double M1, double M2, double M3) 
	{
		double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
       	quat_values = new  float[5][numberOfPts+1];// +1 is for AnimationWindow
	    // create an instance
        RConstantTorque si = new RConstantTorque(time_step,M1,M2,M3,I1,I2,I3,  quat_values);
	    
		// integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        // make the plot visible
        if (plotYes ==1 )
        	si.makePlotsVisible();
        // Animation
        quat_values = si.getQuaternion();
		if (animationYes ==1)
		{
			AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "else");
		}			
	}

    
    public void doCMGManeuver(double[] x0,
    							double J, double A, double psi, double phi, double theta)
    {
    	double tf = timeDuration;
        double t0 = 0.0;               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
		float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
    	// create an instance
        CMGManeuver si = new CMGManeuver(time_step, I1, I2, I3, J, A, psi, phi, theta, quat_values);
		
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
		/// make the plot visible
		if (plotYes == 1)
        	si.makePlotsVisible();
        if (animationYes ==1)
        {
        	AnimationWindow theAnimWindow = new AnimationWindow("Animation", 
															 (float)I1, (float)I2, (float)I3, 
															 numberOfPts, quat_values , "else");
        }
    }
    
    
}
