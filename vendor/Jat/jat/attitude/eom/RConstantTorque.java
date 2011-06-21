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
 * This class contains the equations of motion of a rigid spacecraft
 * subjected to constant torque. The class implements 
 * jat.alg.integrators.Derivatives and jat.alg.integrateors.Printble, 
 * so an outside application code can perform a numerical simulation
 * of equations of motion defined in this class and get the output 
 * from the simulation.
 *
 * @author Noriko Takada
 * @version 1.8 (8/14/2004)
 * Modification since the last version
 * 		Switched to the interface EquationsOfMotion
 * 		 
 */

public class RConstantTorque implements EquationsOfMotion
{
	
	//		time_step:			Time step of numerical integration
 	//		quat_values[][]		Two dimensional array that contains quarternions
 	//							from simulation
 	//		M1,M2,M3			Constant, body fixed torques
 	//		I1,I2,I3			Principal moments of inertia
 	//		rotation_plot		Plots of angular velocities
 	//		angle_plot			Plots of euler angles
 	//		quarternion_check	Plot of e1^2 + e2^2 + e3^2 +e4^2
 	//		quat_plot			Plots of quaternions
 	//		energy_and_angularM	Plots of system energy and angular momentum
 	//  
	double  time_step;
	private float quat_values[][];
	
	private double M1 = 0;           // External torque
    private double M2 = 0;
    private double M3 = 0; 
    private double I1 = 10.42;       // Spacecraft Principle inertia
    private double I2 = 35.42;
    private double I3 = 41.67;
    
    // Create variables for the necessary plots
	private ThreePlots angular_velocity_plot = new ThreePlots();								
	private ThreePlots frame_angle_plot = new ThreePlots();
	private SinglePlot quarternion_check = new SinglePlot();
	private FourPlots quaternion_plot = new FourPlots();
	private SinglePlot energy_plot = new SinglePlot();
	private FourPlots angular_momentum_plot = new FourPlots();
	private ThreePlots nutation_plot = new ThreePlots();
	
	/**
	 * Constructor
	 * @param		time_step:			Time step of numerical integration
	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
	 */
	public RConstantTorque(double time_step, float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
	}
	
	
	/** Construct a RConstantTorque object
	 * @param		time_step:			Time step of numerical integration
 	 * @param		quat_values[][]		Two dimensional array that contains quarternions from simulation
 	 * @param		M1,M2,M3			External Torques about 1,2,and 3 axis
 	 * @param		I1,I2,I3			Principal moments of inertia
	 */
	public RConstantTorque(double time_step,
							double M1, double M2, double M3, double I1, double I2, double I3,
							float quat_values[][])
	{
		setupPlots();
		this.time_step = time_step;
		this.quat_values = quat_values;
		this.M1 = M1;
		this.M2 = M2;
		this.M3 = M3;
		this.I1 = I1;
		this.I2 = I2;
		this.I3 = I3;
	}
	
	
	/**
	 * setupPlots() sets up Plots
	 */
	void setupPlots()
	{
		//ThreePlots angular_velocity_plot = new ThreePlots();								
		//ThreePlots frame_angle_plot = new ThreePlots();
		//SinglePlot quarternion_check = new SinglePlot();
		//FourPlots quaternion_plot = new FourPlots();
		//SinglePlot energy_plot = new SinglePlot();
		//FourPlots angular_momentum_plot = new FourPlots();
		//ThreePlots nutation_plot = new ThreePlots();
		
		// Setup plots
		angular_velocity_plot.setTitle("Angular Velocities");
        angular_velocity_plot.topPlot.setXLabel("t(sec)");
        angular_velocity_plot.topPlot.setYLabel("w1(rad/sec)");
        angular_velocity_plot.middlePlot.setXLabel("t(sec)");
		angular_velocity_plot.middlePlot.setYLabel("w2(rad/sec)");
        angular_velocity_plot.bottomPlot.setXLabel("t(sec)");
        angular_velocity_plot.bottomPlot.setYLabel("w3(rad/sec)");
        
        frame_angle_plot.setTitle("Angles between body and reference frame");
        frame_angle_plot.topPlot.setXLabel("t(sec)");
        frame_angle_plot.topPlot.setYLabel("angle11(degrees)");
		frame_angle_plot.middlePlot.setXLabel("t(sec)");
		frame_angle_plot.middlePlot.setYLabel("angle22(degrees)");
        frame_angle_plot.bottomPlot.setXLabel("t(sec)");
        frame_angle_plot.bottomPlot.setYLabel("angle33(degrees)");
        
        quarternion_check.setTitle("Quarternion Check");
        quarternion_check.plot.setXLabel("t(sec)");
        quarternion_check.plot.setYLabel("e1^2 + e2^2 + e3^2 + e4^2");
        
        quaternion_plot.setTitle("Quarternions");
        quaternion_plot.firstPlot.setXLabel("t (sec)");
        quaternion_plot.firstPlot.setYLabel("q1");
        quaternion_plot.secondPlot.setXLabel("t (sec)");
        quaternion_plot.secondPlot.setYLabel("q2");
        quaternion_plot.thirdPlot.setXLabel("t (sec) ");
        quaternion_plot.thirdPlot.setYLabel("q3");
        quaternion_plot.fourthPlot.setXLabel("t (sec)");
        quaternion_plot.fourthPlot.setYLabel("q4");
        
        energy_plot.setTitle("Energy ");
        energy_plot.plot.setXLabel("t (sec)");
        energy_plot.plot.setYLabel("Energy");
        
        angular_momentum_plot.setTitle("Angular Momentum");
        angular_momentum_plot.firstPlot.setXLabel("t (sec)");
        angular_momentum_plot.firstPlot.setYLabel("Hi (b1)");
        angular_momentum_plot.secondPlot.setXLabel("t (sec)");
        angular_momentum_plot.secondPlot.setYLabel("Hi (b2)");
        angular_momentum_plot.thirdPlot.setXLabel("t (sec)");
        angular_momentum_plot.thirdPlot.setYLabel("Hi (b3)");
        angular_momentum_plot.fourthPlot.setXLabel("t (sec)");
        angular_momentum_plot.fourthPlot.setYLabel("Hi total");
        
        nutation_plot.setTitle("Nutation Angles");
        nutation_plot.topPlot.setXLabel("t (sec)");
        nutation_plot.topPlot.setYLabel("H vector & b3 axis");
        nutation_plot.middlePlot.setXLabel("t (sec)");
        nutation_plot.middlePlot.setYLabel("H vector & b2 axis");
        nutation_plot.bottomPlot.setXLabel("t (sec)");
        nutation_plot.bottomPlot.setYLabel("H vector & b1 axis");
        
          
                
	}
	
	
	/** Compute the derivatives.
     * Equations of Motion
     * @param t    double containing time or the independent variable.
     * @param x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */
		
	public double[] derivs(double t, double[] x)
    {
       	        
       	double [] out = new double[7];
       	out[0] = ((I2-I3)/I1)*x[1]*x[2] + M1/I1 ;
       	out[1] = ((I3-I1)/I2)*x[0]*x[2] + M2/I2 ;
       	out[2] = ((I1-I2)/I3)*x[0]*x[1] + M3/I3 ;
       	out[3] = -0.5* (-x[2]*x[4] + x[1]*x[5] - x[0]*x[6]);
       	out[4] = -0.5* (x[2]*x[3]  - x[0]*x[5] - x[1]*x[6]);
       	out[5] = -0.5* (-x[2]*x[6] + x[0]*x[4] - x[1]*x[3]);
       	out[6] = -0.5* (x[2]*x[5]  + x[1]*x[4] + x[0]*x[3]);
               
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
        	
       	 	double w1 = y[0];
        	double w2 = y[1];
        	double w3 = y[2];
        	double q1 = y[3];
        	double q2 = y[4];
        	double q3 = y[5];
        	double q4 = y[6];
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
        	
        	double angle11 = Math.toDegrees(Math.acos(c11));
        	double angle22 = Math.toDegrees(Math.acos(c22));
        	double angle33 = Math.toDegrees(Math.acos(c33));
        	// See Bong Wie (p. 342)
        	double energy = (0.5)*(I1*w1*w1 + I2*w2*w2 + I3*w3*w3);
        	double Hb1 = I1*w1;
        	double Hb2 = I2*w2;
        	double Hb3 = I3*w3;
        	
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
			
			double nutation_Hb3 = I3*w3/Hi_total;
			double nutation_Hb2 = I2*w2/Hi_total;
			double nutation_Hb1 = I1*w1/Hi_total;
			
			// add data point to the plot
        	//ThreePlots angular_velocity_plot = new ThreePlots();								
			//ThreePlots frame_angle_plot = new ThreePlots();
			//SinglePlot quarternion_check = new SinglePlot();
			//FourPlots quaternion_plot = new FourPlots();
			//SinglePlot energy_plot = new SinglePlot();
			//FourPlots angular_momentum_plot = new FourPlots();
			//ThreePlots nutation_plot = new ThreePlots();
		
        	// Angular Velocity
        	angular_velocity_plot.topPlot.addPoint(0, t, w1, first);
			angular_velocity_plot.middlePlot.addPoint(0, t, w2, first);
        	angular_velocity_plot.bottomPlot.addPoint(0, t, w3, first);
        	
        	quarternion_check.plot.addPoint(0, t, quat_check, first);
        	
        	frame_angle_plot.topPlot.addPoint(0, t, angle11, first);
        	frame_angle_plot.middlePlot.addPoint(0, t, angle22, first);
        	frame_angle_plot.bottomPlot.addPoint(0, t, angle33, first);
        	
        	quaternion_plot.firstPlot.addPoint(0, t, q1, first);
        	quaternion_plot.secondPlot.addPoint(0, t, q2, first);
        	quaternion_plot.thirdPlot.addPoint(0, t, q3, first);
        	quaternion_plot.fourthPlot.addPoint(0, t, q4, first);
        	
        	energy_plot.plot.addPoint(0, t, energy, first);
        	   		
        	angular_momentum_plot.firstPlot.addPoint(0, t, Hi1, first);
			angular_momentum_plot.secondPlot.addPoint(0, t, Hi2, first);
			angular_momentum_plot.thirdPlot.addPoint(0, t, Hi3, first);
			angular_momentum_plot.fourthPlot.addPoint(0, t, Hi_total, first);     	
        	
        	nutation_plot.topPlot.addPoint(0, t, nutation_Hb3, first);
        	nutation_plot.middlePlot.addPoint(0, t, nutation_Hb2, first);
        	nutation_plot.bottomPlot.addPoint(0, t, nutation_Hb1, first);
        	
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
   		angular_velocity_plot.setVisible(true);							
		frame_angle_plot.setVisible(true);
		quarternion_check.setVisible(true);
		quaternion_plot.setVisible(true);
		energy_plot.setVisible(true);
		angular_momentum_plot.setVisible(true);
		nutation_plot.setVisible(true);
   }
   
   /** Runs the example.
    * @param args Arguments.
    */
    public static void main(String[] args)
    {

        double time_step=0.1;
        double timeDuration=10;
        double tf = 10;
        double t0 = 0.0;
               
        RungeKutta8 rk8 = new RungeKutta8(time_step);
		timeDuration=tf;	// Duration of the simulation time is the same as the final time
		
		int numberOfPts = (int)(timeDuration/time_step) +1 ;  				
    
    	float quat_values[][]= new  float[5][numberOfPts+1];// +1 is for AnimationWindow
		
		 
        // create an instance
        RConstantTorque si = new RConstantTorque(time_step,  quat_values);
		
        // initialize the variables
        double [] x0 = new double[7];
        x0[0] = 0.0;
        x0[1] = 0.1;
        x0[2] = 0.0;
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
    		