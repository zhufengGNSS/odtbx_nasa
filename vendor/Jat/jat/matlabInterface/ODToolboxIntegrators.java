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
 *  Stephen D Metcalfe  09/27/2010   	Mantiss Issue 0000210
 *  									Added HermiteInterpolator and ddHermite
 *  Stephen D Metcalfe  09/29/2010   	Mantiss Issue 0000303
 *  									Consolidated interpolators in JATInterpolators
 */

package jat.matlabInterface;

import jat.alg.integrators.Derivatives;
import jat.alg.integrators.RungeKutta4;
import jat.alg.integrators.RungeKutta8;
import jat.alg.integrators.RungeKuttaFehlberg78;

import java.util.Arrays;

/**
 * The ODToolboxIntegrators.java Class is an adaptor for Matlab to call the JAT Integrators in the same 
 * way that Matlab calls it's own integrators. These integrators will only work with JAT force models.
 *
 * @author Derek Surka, Emergent Space Technologies
 */

public class ODToolboxIntegrators extends JATIntegrators {

	/** This method calls the RK8 integrator based on the initial conditions from Matlab
	 * @param inputTime - user defined times
	 * @param x0 		- initial state
	 * @param stepSize	- Integrator step size
	 * @param derivs	- ODToolboxJATModel object that contains force information
	 * @return  		- Return the states at the requested input times
	 */
	public static double[][] RK8(double[] inputTime, double[] x0, double stepSize, Derivatives derivs)
	{
		int              neqns = x0.length;
		double[]          xNew = new double[neqns];
		double[]          xOld = new double[neqns];
		RungeKutta8 integrator = new RungeKutta8(stepSize);

		// Generate output time vector based on user-provided inputs
		double tSimulation = inputTime[inputTime.length-1] - inputTime[0];
		int timeLength, i;

		if ( tSimulation%stepSize >= 1.0 )
			timeLength = (int) (tSimulation/stepSize)+2;
		else
			timeLength = (int) (tSimulation/stepSize)+1;

		double[]  outputTime = new double[timeLength];
		double[][] simStates = new double[timeLength][neqns];
		double[]     simTime = createStepTime(inputTime, stepSize, timeLength);

		if (inputTime.length<3)
			outputTime = simTime; // if user did not specify timesteps, use integrator steps
		else {
			outputTime = new double[inputTime.length];
			outputTime = inputTime;
		}

		// store the initial states
		for (i=0; i < neqns; i++) {
			simStates[0][i]= x0[i];
			xOld[i]        = x0[i];
		}

		/* Step the integrator and store new states. */
		for (int tCounter=1; tCounter<timeLength; tCounter++) {
			xNew = integrator.step(simTime[tCounter], xOld, derivs);

			for (i = 0; i < neqns; i++) {
				simStates[tCounter][i] = xNew[i];
				xOld[i]                = xNew[i];
			}
		}

		return interpolateStates(outputTime, simTime, simStates, neqns);

	}

	/** This method calls the RK4 integrator based on the initial conditions from Matlab
	 * @param inputTime - user defined times
	 * @param x0 		- initial state
	 * @param stepSize	- Integrator step size
	 * @param derivs	- ODToolboxJATModel object that contains force information
	 * @return  		- Return the states at the requested input times
	 */
	public static double[][] RK4(double[] inputTime, double[] x0, double stepSize, Derivatives derivs)
	{
		int              neqns = x0.length;
		double[]          xNew = new double[neqns];
		double[]          xOld = new double[neqns];
		RungeKutta4 integrator = new RungeKutta4(stepSize);

		// Generate output time vector based on user-provided inputs
		double tSimulation = inputTime[inputTime.length-1] - inputTime[0];
		int timeLength, i;

		if ( tSimulation%stepSize >= 1.0 )
			timeLength = (int) (tSimulation/stepSize)+2;
		else
			timeLength = (int) (tSimulation/stepSize)+1;

		double[]  outputTime = new double[timeLength];
		double[][] simStates = new double[timeLength][neqns];
		double[]     simTime = createStepTime(inputTime, stepSize, timeLength);

		if (inputTime.length<3)
			outputTime = simTime; // if user did not specify timesteps, use integrator steps
		else {
			outputTime = new double[inputTime.length];
			outputTime = inputTime;
		}

		// store the initial states
		for (i=0; i < neqns; i++) {
			simStates[0][i]= x0[i];
			xOld[i]        = x0[i];
		}

		/* Step the integrator and store new states. */
		for (int tCounter=1; tCounter<timeLength; tCounter++) {
			xNew = integrator.step(simTime[tCounter], xOld, derivs);

			for (i = 0; i < neqns; i++) {
				simStates[tCounter][i] = xNew[i];
				xOld[i]                = xNew[i];
			}
		}

		return interpolateStates(outputTime, simTime, simStates, neqns);

	}

	/** This method calls the RKF78 integrator based on the initial conditions from Matlab
	 * @param inputTime 	- user defined times
	 * @param x0 			- initial state
	 * @param stepSize		- Integrator step size
	 * @param minStepSize	- Integrator minimum step size
	 * @param accuracy		- Integrator accuracy tolerance
	 * @param derivs		- ODToolboxJATModel object that contains force information
	 * @return  			- Return the states at the requested input times
	 */
	public static double[][] RKF78(double[] inputTime, double[] x0, double stepSize, double minStepSize, double accuracy, Derivatives derivs)
	{
		double[] outputTime;
		RungeKuttaFehlberg78 integrator = new RungeKuttaFehlberg78();

		// Set integrator parameters
		integrator.setAccuracy(accuracy);
		integrator.setMinimumStepSize(minStepSize);
		integrator.setStepSize(stepSize);
		integrator.setAdaptive();
		integrator.setSaveSteps(true);

		integrator.integrate(inputTime[0], x0, inputTime[inputTime.length-1], derivs);
		double[]     simTime = integrator.getIntermediateTimes();
		double[][] simStates = integrator.getIntermediateStates();

		// Get output times
		if (inputTime.length<3)
			outputTime = simTime; // if user did not specify timesteps, use integrator steps
		else {
			outputTime = inputTime;
		}

		return interpolateStates(outputTime, simTime, simStates, x0.length);

	}

	private static double[] createStepTime(double[] inputTime, double stepSize, int timeLength)
	{
		double[] stepTime 	= new double[timeLength];
		stepTime[0] 		= inputTime[0];
		double tFinal 		= inputTime[inputTime.length-1];
		double stop 		= tFinal - (stepSize/1000);

		double newStep;
		int counter = 0;
		while (stepTime[counter]<stop) {
			if ( (stepTime[counter]+stepSize) > tFinal ) {
				newStep = tFinal - stepTime[counter];				
			}
			else {
				newStep = stepSize;
			}
			stepTime[counter+1] = stepTime[counter]+newStep;
			counter ++;        			
		}

		return stepTime;
	}

	/** Determines if states and time need to be interpolated and then calls the 
	 *  interpolator to get the states corresponding to the user defined times. 
	 *  Also packages time and sim states into one array formatted per Matlab ode outputs
	 *  	Each row is a different timestep
	 *  	First column is time
	 *  	Remaining columns are states */	 
	private static double[][] interpolateStates(double[] outputTime, double[] simTime, double[][] simStates, int neqns)
	{
		int timeSpot = 0;
		double currentTime; 
		double[] interpolatedState = new double[neqns];
		double[][] newOutput = new double[outputTime.length][neqns+1];

		for (int i=0; i<outputTime.length; i++) {
			currentTime     = outputTime[i];
			newOutput[i][0] = currentTime;

			timeSpot = Arrays.binarySearch(simTime, currentTime);
			if (timeSpot < 0) {
				if (simTime.length < 1) {
					throw new java.lang.IndexOutOfBoundsException("Number of integrator steps = " 
						+ simTime.length + "\nThe State Interpolator requires at least 1 integration time step.");
				} 
				else if (simTime.length < 8) {
					interpolatedState = JATIntegrators.HermiteInterpolator(simTime, currentTime, neqns, simStates);
				}
				else {
					interpolatedState = JATIntegrators.LagrangianInterpolator(simTime, currentTime, neqns, simStates);
				}
				for (int k = 0; k<neqns; k++)
					newOutput[i][k+1] = interpolatedState[k];
			} else {
				for (int k = 0; k<neqns; k++)
					newOutput[i][k+1] = simStates[timeSpot][k];
			}
		}	

		return newOutput;
	}

}
