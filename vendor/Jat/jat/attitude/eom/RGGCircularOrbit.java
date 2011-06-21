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
 * This class contains the equations of motion for the simulation
 * of a rigid spacecraft subjected to gravity gradient torque
 * while it orbits around the earth.
 * The class implements jat.alg.integrators.Derivatives and
 * jat.alg.integrateors.Printble, so an outside application code
 * can perform a numerical simulation of equations of motion defined
 * in this class and get the output from the simulation.
 *  
 * @author Noriko Takada
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion * 
 * 
 * @param	time_step:			Time step of numerical integration
 * @param	quat_values[][]		Two dimensional array that contains quarternions from simulation
 * @param	I1,I2,I3			Principal moments of inertia about 1,2, and 3 axes
 * 
 */

public class RGGCircularOrbit implements EquationsOfMotion
{
	double  time_step;
	private float quat_values[][];
	
	// Create variables for the necessary plots
	private ThreePlots angular_velocity_plot = new ThreePlots();								
	private ThreePlots euler_angle_plot = new ThreePlots();
	private SinglePlot quarternion_check = new SinglePlot();
	private FourPlots quaternion_plot = new FourPlots();
	private SinglePlot energy_plot = new SinglePlot();
	private FourPlots angular_momentum_plot = new FourPlots();
	
	
	
	private double I1 = 10.42;       
    private double I2 = 35.42;
    private double I3 = 41.67;
	
	
	/**
	 * Constructor 1
	 * @param	time_step:			Time step of numerical integration
 	 * @param	quat_values[][]		Two dimensional array that contains quarternions from simulation
	 * @param 	I1					Principle moment of inertia about 1 axis
	 * @param 	I2					Principle moment of inertia about 2 axis
	 * @param 	I3					Principle moment of inertia about 3 axis
	 */ 
	public RGGCircularOrbit(double time_step, double I1, double I2, double I3,
							 float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
	}
	
	/**
	 * Constructor 2
	 * @param	time_step:			Time step of numerical integration
 	 * @param	quat_values[][]		Two dimensional array that contains quarternions from simulation
 	 */
	public RGGCircularOrbit(double time_step, float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
	}
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		//ThreePlots angular_velocity_plot = new ThreePlots();								
		//ThreePlots euler_angle_plot = new ThreePlots();
		//SinglePlot quarternion_check = new SinglePlot();
		//FourPlots quaternion_plot = new FourPlots();
		//SinglePlot energy_plot = new SinglePlot();
		//FourPlots angular_momentum_plot = new FourPlots();
		
		// Setup plots
		angular_velocity_plot.setTitle("Angular Velocities");
        angular_velocity_plot.topPlot.setXLabel("t(orbits)");
        angular_velocity_plot.topPlot.setYLabel("w1(rad/orbits)");
        angular_velocity_plot.middlePlot.setXLabel("t(sec)");
		angular_velocity_plot.middlePlot.setYLabel("w2(rad/orbits)");
        angular_velocity_plot.bottomPlot.setXLabel("t(orbits)");
        angular_velocity_plot.bottomPlot.setYLabel("w3(rad/orbits)");
        
        euler_angle_plot.setTitle("Euler angles between body frame & orbit frame");
        euler_angle_plot.topPlot.setXLabel("t(orbits)");
        euler_angle_plot.topPlot.setYLabel("Psi (degrees)");
		euler_angle_plot.middlePlot.setXLabel("t(orbits)");
		euler_angle_plot.middlePlot.setYLabel("Theta (degrees)");
        euler_angle_plot.bottomPlot.setXLabel("t(orbits)");
        euler_angle_plot.bottomPlot.setYLabel("Phi (degrees)");
        
        quarternion_check.setTitle("Quarternion Check");
        quarternion_check.plot.setXLabel("t(sec)");
        quarternion_check.plot.setYLabel("q1^2 + q2^2 + q3^2 + q4^2");
        
        quaternion_plot.setTitle("Quarternions");
        quaternion_plot.firstPlot.setXLabel("t (Orbits)");
        quaternion_plot.firstPlot.setYLabel("q1");
        quaternion_plot.secondPlot.setXLabel("t (Orbits)");
        quaternion_plot.secondPlot.setYLabel("q2");
        quaternion_plot.thirdPlot.setXLabel("t (Orbits) ");
        quaternion_plot.thirdPlot.setYLabel("q3");
        quaternion_plot.fourthPlot.setXLabel("t (Orbits)");
        quaternion_plot.fourthPlot.setYLabel("q4");
        
        energy_plot.setTitle("Energy ");
        energy_plot.plot.setXLabel("t (Orbits)");
        energy_plot.plot.setYLabel("Energy");
        
        angular_momentum_plot.setTitle("Angular Momentum");
        angular_momentum_plot.firstPlot.setXLabel("t (Orbits)");
        angular_momentum_plot.firstPlot.setYLabel("Hi (b1)");
        angular_momentum_plot.secondPlot.setXLabel("t (Orbits)");
        angular_momentum_plot.secondPlot.setYLabel("Hi (b2)");
        angular_momentum_plot.thirdPlot.setXLabel("t (Orbits)");
        angular_momentum_plot.thirdPlot.setYLabel("Hi (b3)");
        angular_momentum_plot.fourthPlot.setXLabel("t (Orbits)");
        angular_momentum_plot.fourthPlot.setYLabel("Hi total");
        
              
	}
	
	
	/** Compute the derivatives.
     * Equations of Motion
     * @params t    double containing time or the independent variable.
     * @params x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
	public double[] derivs(double t, double[] x)
    {
       	double c11 = 1- 2*( x[4]*x[4] + x[5]*x[5]);
        double c21 = 2* (x[3]*x[4]-x[5]*x[6]);
        double c31 = 2* (x[3]*x[5]+x[4]*x[6]); 
        
        double [] out = new double[7];
        out[0] = 2*Math.PI*((I2-I3)/I1)* (x[1]*x[2] - 3*c21*c31);
        out[1] = 2*Math.PI*((I3-I1)/I2)* (x[0]*x[2] - 3*c31*c11) ;
        out[2] = 2*Math.PI*((I1-I2)/I3)* (x[0]*x[1] - 3*c11*c21) ;
        out[3] = -Math.PI* (-(x[2]+1)*x[4] + x[1]*x[5] - x[0]*x[6]);
        out[4] = -Math.PI* ((x[2]+1)*x[3]  - x[0]*x[5] - x[1]*x[6]);
        out[5] = -Math.PI* (-(x[2]-1)*x[6] + x[0]*x[4] - x[1]*x[3]);
        out[6] = -Math.PI* ((x[2]-1)*x[5]  + x[1]*x[4] + x[0]*x[3]);
        
        
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
        
        //ThreePlots angular_velocity_plot = new ThreePlots();								
		//ThreePlots euler_angle_plot = new ThreePlots();
		//SinglePlot quarternion_check = new SinglePlot();
		//FourPlots quaternion_plot = new FourPlots();
		//SinglePlot energy_plot = new SinglePlot();
		//FourPlots angular_momentum_plot = new FourPlots();
        // add data point to the plot
        angular_velocity_plot.topPlot.addPoint(0, t,y[0], first);
        angular_velocity_plot.middlePlot.addPoint(0, t, y[1], first);
		angular_velocity_plot.bottomPlot.addPoint(0, t, y[2], first);
		
		euler_angle_plot.topPlot.addPoint(0, t, Theta, first);
		euler_angle_plot.middlePlot.addPoint(0, t, Psi, first);
		euler_angle_plot.bottomPlot.addPoint(0, t, Phi, first);  
		
		quarternion_check.plot.addPoint(0, t, quat_check, first);
		
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
   		//ThreePlots angular_velocity_plot = new ThreePlots();								
		//ThreePlots euler_angle_plot = new ThreePlots();
		//SinglePlot quarternion_check = new SinglePlot();
		//FourPlots quaternion_plot = new FourPlots();
		//SinglePlot energy_plot = new SinglePlot();
		//FourPlots angular_momentum_plot = new FourPlots();
   		angular_velocity_plot.setVisible(true);
   		euler_angle_plot.setVisible(true);
   		quarternion_check.setVisible(true);
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
        RGGCircularOrbit si = new RGGCircularOrbit(time_step, quat_values);
		
        // initialize the variables
        double [] x0 = new double[7];
        x0[0] = 0.0;
        x0[1] = 0.0;
        x0[2] = 1.0;
        x0[3] = 0.0;
        x0[4] = 0.0;
        x0[5] = 0.0;
        x0[6] = 1.0;
        
        // integrate the equations
        rk8.integrate(t0, x0, tf, si, true);
        
        // make the plot visible
        si.makePlotsVisible();
    }
    							
	   
    	
}// End of File
    		