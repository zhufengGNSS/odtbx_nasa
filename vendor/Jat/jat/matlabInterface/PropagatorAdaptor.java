/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
 *  REVISION HISTORY
 *  Author      		Date         	Comment
 *              	   (MM/DD/YYYY)
 *  Stephen D Metcalfe  09/29/2010   	Mantiss Issue 0000303
 *  									Consolidated interpolators in JATInterpolators
 *  									Removed extraneous import of java.util.Vector 
 *										Commented out unused SimModel JGM2All in runJGMSimulation
 */

package jat.matlabInterface;

import java.util.Arrays;
import jat.alg.integrators.Derivatives;
import jat.alg.integrators.RungeKutta8;
import jat.matvec.data.VectorN;
import jat.sim.SimModel;
import jat.spacecraft.SpacecraftModel;

/**
 * The JatAdaptor.java Class is an adaptor for Matlab to call the JAT
 * Integrators in the same way that Matlab calls its own integrators.
 * 
 * The times specified by the client for output are not necessarily the times
 * that the propagator uses to meet its accuracy requirements as it steps, e.g.
 * its step size.
 * 
 * @author Emergent Space Technologies
 */
public class PropagatorAdaptor extends JATIntegrators{
	
	/* 
	 * Design notes:  Client specification of time and integrator steps
	 * 
	 * The stepSize does not have to be, and is not assumed to be, the same as 
	 * the time spacing in the given time vector
	 *
	 * The stepSize controls the integrator steps. In all logic paths for 
	 * jat.matlabInterface.PropagatorAdaptor's RK8(), JAT uses the 
	 * jat.alg.integrators.RungeKutta8 class. This means that the actual steps 
	 * the integrator takes is directly dependent on the step size.
	 * 
	 * However, the time vector specified from MATLAB is a set of times 
	 * for desired output... not the exact times the integrator steps to as it 
	 * integrates. This decoupling of the integrator and the outputs is also how 
	 * the MATLAB ode** framework is also designed. 
	 * 
	 * After the integration occurs in the PropagatorAdaptor, the 
	 * PropagatorAdaptor calls its checkOutput(), which calls its 
	 * interpolateStates() if the time vector has more elements than [t0 tf].
	 *  
	 * The interpolation that takes place uses the PropagatorAdaptor's 
	 * LagrangianInterpolator(), a "7th Order Lagrangian Interpolator adapted 
	 * from Steve Hughes Matlab Code", to interpolate the states computed at 
	 * each actual integrated time step into states at the requested output 
	 * time points.
	 */

	
	/** This method calls the RK8 integrator based on the initial conditions from Matlab 
	 * for the Matlab defined EOMs
	 * 
	 * @param DerivativesName -- What is the name of the Matlab EOM file
	 * @param initialTime -- Time increments from the user for data output
	 * @param x0 -- Initial State
	 * @param stepSize -- Step size for the integrator (NOT tied to initialTime)
	 * @return -- Return the states based on initialTime points to Matlab
	 */
	 public static double[][] RK8(String DerivativesName, double [] initialTime, double [] x0, double stepSize)
	    {		 	
		    int timeLength = calcStepTimeLength(initialTime, stepSize);
		    
		 	double[] time = new double[timeLength];
		 	double[] stepTime = new double[timeLength];

       		stepTime = createStepTime(initialTime[0], stepSize, timeLength);
       		
		 	if (initialTime.length<3)
	       	{
	       		time = stepTime;
	       	}
	       	else 
	       	{
	       		time = new double[initialTime.length];
	       		time = initialTime;
	       	}
		 	
		    // Initialize start and finish times
	    	int neqns = x0.length;
	    	double [][] oldState = new double [time.length][neqns+1];
	    	
	    	// Use if calling from java
	        //JatDerivatives MD = new JatDerivatives();
	    	// Use if calling from Matlab
	    	MatlabDerivs MD = new MatlabDerivs(DerivativesName);

        	// Send to integrate the equations
        	oldState = jat_RK8(time, x0, stepSize, MD);
	    
        	return oldState;
	    }

	/** This method calls the RK8 integrator based on the initial conditions from Matlab 
	 * for the JGM2 & JGM3 cases
	 * 
	 * @param DerivativesName -- Whether it is JGM2 or JGM3 test case
	 * @param initialTime -- User defined times for data output
	 * @param x0 -- Initial State
	 * @param stepSize -- Step size for the integrator (NOT tied to initialTime)
	 * @param cd -- Coefficient of drag
	 * @param cr -- Coefficient of reflectivity
	 * @param mass -- Mass of the satellite
	 * @param cArea -- Cross sectional area
	 * @param mjd_utc -- Mean Julian Date in UTC time
	 * @param path -- Path to JAT
	 * @param matOrder -- added by DMS on 5/18/07
	 * @param matDegree -- added by DMS on 5/18/07
	 * @return -- Return the states based on initialTime points to Matlab
	 */
	 public static double[][] RK8(String DerivativesName, double [] initialTime, double [] x0, double stepSize, double cd, double cr, double mass, double cArea, double mjd_utc, String path, double matOrder, double matDegree)
	    {
	     	int degree = (int) matDegree;  	// added by DMS 5/18/07
	     	int order = (int) matOrder;		// added by DMS 5/18/07
			
			int timeLength = calcStepTimeLength(initialTime, stepSize);
			
			double[] time = new double[timeLength];
			double[] stepTime = new double[timeLength];

    		stepTime = createStepTime(initialTime[0], stepSize, timeLength);
    		
		 	if (initialTime.length<3)
        	{
        		time = stepTime;
        	}
        	else 
        	{
        		time = new double[initialTime.length];
        		time = initialTime;
        	}
		 	
		    // Initialize start and finish times
	    	int neqns = x0.length;
	    	double[][] output = new double[stepTime.length][neqns+1];
	    	double [][] oldState = new double [time.length][neqns+1];
	    	
	    	boolean use_JGM2 = false;
	        if (DerivativesName.equals("JatUniverseJGM2"))
	        {
	        	use_JGM2 = true;
	        }
	        
        	// add order,degree to runJGMSimulation call (DMS 5/18/07)
        	output = runJGMSimulation(x0,stepTime, mjd_utc, stepSize, cd, cr, mass, cArea, use_JGM2, neqns, path, order, degree);	
        	oldState = checkOutput(output, initialTime, neqns);

		    return oldState;
		        
	    }
	 
	/**
	 * Calculates the length of the step time array used by the integrator
	 * to adequately cover the desired time span.  Used to determine
	 * the length of the array populated by createStepTime().
	 * @param initialTime Array of times specified for output purposes 
	 * @param stepSize The desired integrator step size
	 * @return Length of the time array that should be populated by createStepTime()
	 */
	private static int calcStepTimeLength(double[] initialTime, double stepSize) {
		double tspan = initialTime[initialTime.length-1]-initialTime[0];
		
		// truncates, at least one for holding initialTime[0]:
		int timeLength = (int)(tspan/stepSize) + 1; 
		
		if (timeLength < 2) {
			timeLength = 2; // timeLength must always be > 2, for [t0 t0+stepSize]
		}
		
		if (((timeLength-1)*stepSize + initialTime[0]) < initialTime[initialTime.length-1]) {
		    /* 
		     * Add one more time point to ensure the last propagated time  
			 * is at or after the output's initialTime[initialTime.length-1]
			 */
			timeLength++; 
		}
		return timeLength;
	}
	 
	 /**
	  * Creates the times that the propagator will step at based on the given
	  * initial time, requested array length, and specified step size.  The
	  * time units are not required so long as the time units in initialTime 
	  * and stepSize are consistent.
	  * @param initialTime The initial time point to step from (time units)
	  * @param stepSize The size of each time step for the propagator (time units)
	  * @param timeLength The length of the requested vector (>0)
	  * @return A vector of times with equal stepSize spacing.
	  */
	 private static double[] createStepTime(double initialTime, double stepSize, int timeLength)
	 {
		double[] stepTime = new double[timeLength];

		for (int i = 0; i < timeLength; i++) {
			stepTime[i] = (stepSize * i) + initialTime;
		}
		
		return stepTime;
	 }
	 
	 /**
	  * Interpolates the computed states and times into the requested times for 
	  * the user.  If the requested times, time, is only two elements, then
	  * the computed states are returned instead.
	  * @param output The computed states at specified times from the integrator
	  * @param time The requested times for output
	  * @param neqns The number of equations
	  * @return Interpolated states at requested times.
	  */
	 private static double[][] checkOutput(double[][] output, double[] time, int neqns)
	 {
		 if(time.length<3)
		 {
		 }
		 /* When user specifies output times the code comes here to get the corresponding
	      * states. */
	     else
	     {
	    	 double[] timeArray = new double[output.length];
	    	 double[][] newStates = new double[output.length][neqns];

	    	 for (int i=0; i<output.length; i++)
	    	 {
	    		 timeArray[i] = output[i][0];
	    		 
		    	 for (int j=0; j<neqns; j++)
		    	 {
		    		 newStates[i][j] = output [i][j+1];
		    	 }
		    	 
	    	 }
	    	 output = interpolateStates(time, timeArray, newStates, neqns); 
	     }

		 return output;
	 }

	 // add order, degree DMS 5/18/07
	 private static double[][] runJGMSimulation(double[] x0,double[] time, double mjd_utc, double stepSize, double cd, double cr, double mass, double cArea, boolean use_JGM2, int neqns, String path, int order, int degree)
	 {
		//double start = System.currentTimeMillis();
		double[][] output = new double[neqns][];
//		* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation
//     	boolean[] force_flag = {true,false,false,false,false};  removed DMS 5/18/07
     	boolean[] force_flag = {false,false,false,false,false};  // DMS 5/18/07 for JGM testing
     	
     	VectorN r = new VectorN(x0[0], x0[1], x0[2]);
     	VectorN v = new VectorN(x0[3], x0[4], x0[5]);
     	SimModel sim = new SimModel();
        double t0 = time[0];
        double tf = time[time.length-1];
     	
     	// either NRL or HarrisPriester
     	String drag_model = "NRL";
        SpacecraftModel sm = new SpacecraftModel(r,v,cr,cd,cArea,mass);
     	
//     	SimModel JGM2All = new SimModel();
//     	JGM2All.initialize(sm,t0,tf,mjd_utc, stepsize, 1, out);
//     	JGM2All.initialize(sm,t0,tf,mjd_utc, stepSize);  // added DMS 5/18/07
//     	JGM2All.initializeForces(force_flag, use_JGM2, drag_model, order, degree);  // order, degree added DMS 5/18/07
 
        
        for(int i=0; i<1; i++)
        {
        	sim.initialize(sm,t0,tf,mjd_utc, stepSize);
            String test = drag_model;
            sim.initializeForces(force_flag, use_JGM2, test, order, degree);  // order, degree added DMS 5/18/07
            output = sim.runloopMatlabAdaptor(x0, time);
        }	 
             
         //double elapsed = (System.currentTimeMillis()-start)*0.001/60;
         return output;
	 }
	 
	 /** Actually integrates the equations of motion using the RungeKutta8 step method */
	 private static double[][] jat_RK8(double[] time, double []x0, double stepSize, Derivatives derivs)
	 {
		 int neqns = x0.length;
	     double dt = stepSize;
	     double t = time[0];
	     double tf = time[time.length-1];

	     int timeLength = calcStepTimeLength(time, stepSize);
	
	     double[] newstate = new double[neqns];
	     double[] oldstate = new double[neqns];
	     RungeKutta8 integrate = new RungeKutta8(stepSize);
	     double[] timeArray = new double[timeLength];
	     double[][] newStates = new double[timeLength][neqns];
	     
	     // put initial conditions into the previous state array
	     for (int i = 0; i < neqns; i++) 
	     {
	    	 oldstate[i] = x0[i];
	     }

	     if ((t + dt) > tf) 
	     {
	         dt = (tf - t);
	     }

	     int counter=1;
	     // Set the initial time in the time vector
	     timeArray[0]=t;
	     
	     /* Steps the integrator while t<tf and assigns the states to newState vector and
	      * time to the time vector.
	      * Check both the time and the counter index because floating-point
	      * errors can cause the counter to accidentally overrun */
	     while((t < tf) && (counter < timeArray.length)) {
	    	 newstate = integrate.step(t, oldstate, derivs);
	         for (int i = 0; i < neqns; i++) 
	         {
	        	oldstate[i] = newstate[i];
	         }
	         
	         t = t + dt;
	         
	         timeArray[counter]=t;
	         for (int i = 0; i < neqns; i++)
	         {
	        	 newStates[counter][i]=newstate[i];
	         }
	         counter ++;
	         
	         if ((t + dt) > tf) 
	         {
	        	dt = tf - t;
	         }
	      }
	     
	     // Sets the initial states in the newStates vector.
		 for (int j = 0; j<neqns; j++)
	     {
	    	 newStates[0][j]=x0[j];
	     }

		 // Creates Output matrix that will contain the time and state values. 
	     double[][] newOutput = new double[timeArray.length][neqns+1];

	    	 for (int i = 0; i<timeArray.length; i++)
	    	 {
	    		 for (int k = 0; k<neqns+1; k++)
	    		 {
	    			 if (k==0)
	    			 {
	    				 newOutput[i][k]= timeArray[i];
	    			 }
	    			 else
	    			 {
	    				 newOutput[i][k] = newStates[i][k-1];
	    			 }
	    		 }		 
	    	 }
	    	 newOutput = checkOutput(newOutput, time, neqns);

	      return newOutput;        
}
	 
	 /** Called if user specifies output other than step size based. Figures out if states
	  * and time need to be interpolated and then calls the interpolator to get the states
	  * corresponding to the user defined times. */
	 private static double[][] interpolateStates(double[] OutputTimes, double[] timeArray, double[][] newStates, int neqns)
	 {
    	 int timeSpot = 0;
		 double currentTime; 
    	 double[] interpolatedState = new double[neqns];
    	 double[][] newOutput = new double[OutputTimes.length][neqns+1];

    	 for (int i=0; i<OutputTimes.length; i++)
    		 {
    		 	currentTime = OutputTimes[i];
    		 	//store.clear();
	 			
	 			try
	 			{
	 				timeSpot = Arrays.binarySearch(timeArray, currentTime);
	 				if (timeSpot < 0)
	 				{
	 					//interpolatedState = linearInterpolator(timeArray, currentTime, neqns, newStates);
	 					
	 					interpolatedState = JATIntegrators.LagrangianInterpolator(timeArray, currentTime, neqns, newStates);
	 				}
	 			}
	 			finally
	 			{
	 				
	 			}

	 			for (int k = 0; k<neqns+1; k++)
	 			{
	 				if (k==0)
	 				{
	 					newOutput[i][k]= currentTime;
	 				}
	 				else
	 				{
	 					if(timeSpot<0)
	 					{
	 						newOutput[i][k] = interpolatedState[k-1];
	 					}
	 					else
	 					{
	 						newOutput[i][k] = newStates[timeSpot][k-1];
	 					}
	 				}
	 			}	
    	 }
    	 return newOutput;
	 }
}
