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
 *  REVISION HISTORY
 *  Author      		Date         	Comment
 *              	   (MM/DD/YYYY)
 *  Stephen D Metcalfe  09/29/2010   	Original - Mantiss Issue 0000303
 *  									Consolidated JAT interpolators here.
 *  									Fixed HermiteInterpolator start selection
 *  									per kberry's sugestion. Added check for
 *                                      NaN/Infinity returned from ddHermite.
 *  Stephen D Metcalfe  10/04/2010		Removed HermiteInterpolator test for
 *  									current time < start time
 */

package jat.matlabInterface;

import java.util.Vector;

public class JATIntegrators {

	/** This is the linear interpolator. It will calculate the states based on the current 
	 * time and the states before and after the current state of interest.
	 * @param timeArray Array of times from the integrator.
	 * @param currentTime the current user defined output time.
	 * @param store vector to store the location in the timeArray where the currentTime is located.
	 * @param neqns is the number of states.
	 * @param newStates double array containing the states from the integrator.
	 * @return
	 */
	public static double[] LinearInterpolator(double[] timeArray, double currentTime, int neqns, double[][] newStates)
	{
		double earlierState;
		double laterState;
		double nextTime;
		double earlierTime;
		double slope;
		double stateDifference=0.0;
		double[] interpolatedState = new double[neqns];
		Vector<Integer> store = new Vector<Integer>();

		for(int n=0; n<timeArray.length; n++)
		{
			if (currentTime<timeArray[n])
			{
				store.add(new Integer(n));
			}
			else{
				if (n==(timeArray.length-1))
				{
					store.add(new Integer(timeArray.length-1));
				}

			}
		}
		Integer timePosition = (Integer) store.get(0);
		store.clear();
		nextTime = timeArray[timePosition.intValue()];
		earlierTime = timeArray[timePosition.intValue()-1];

		for (int g =0; g<neqns; g++)
		{
			earlierState = newStates[timePosition.intValue()-1][g];
			laterState = newStates[timePosition.intValue()][g];
			slope = (laterState-earlierState)/(nextTime - earlierTime);
			double earlierDistance = currentTime-earlierTime;
			double laterDistance = nextTime - currentTime;

			if (laterDistance<earlierDistance)
			{
				stateDifference = laterState - (slope*(nextTime-currentTime));
			}
			else
			{
				stateDifference = earlierState + (slope*(currentTime-earlierTime));
			}

			interpolatedState[g] = stateDifference;
		}
		return interpolatedState;
	}

	/** 7th Order Lagrangian Interpolator adapted from Steve Hughes Matlab Code. 
	 * @author Kate Bradley
	 * @date 20 March 2006
	 * @param timeArray is the array of times from the interpolator.
	 * @param currentTime is the current user specified output time.
	 * @param neqns is the number of states.
	 * @param newStates is the double array of states from the interpolator.
	 * @return
	 */
	public static double[] LagrangianInterpolator(double[] timeArray, double currentTime, int neqns, double[][] newStates)
	{
//		int n = timeArray.length;
		double product;
		double [] interpolatedState = new double[neqns];
		int position=0;
		int start;
		int end;
		
		if (timeArray.length < 2) 
		{
			throw new java.lang.IndexOutOfBoundsException("Number of integrator steps = " 
					+ timeArray.length + "\nThe Lagrangian Interpolator requires at least 2 integration time steps.");
		} 
		
		if (timeArray.length < 8)
		{
			interpolatedState = LinearInterpolator(timeArray, currentTime, neqns, newStates);
			return interpolatedState;
		}
		
		for (int t=0; t<timeArray.length-1; t++)
		{
			if((currentTime > timeArray[t]) && (currentTime <timeArray[t+1]))
			{
				if (currentTime > timeArray[timeArray.length-1])
				{
					position = t;
				}
				else 
				{
					position=t+1;
				}
				break;
			}
		}
		if (currentTime > timeArray[timeArray.length-1])
		{
			// Use Linear Extrapolation
			for (int g=0; g<neqns; g++)
			{
				double state = newStates[timeArray.length-1][g];
				double yi = (state/timeArray[timeArray.length-1])*currentTime;
				interpolatedState[g] = yi;
			}
		}
		else
		{
			if (position<5)
			{
				/*Means we are at the beginning of the array and we need to take values 
				 * from the beginning of the array.
				 */
				start = 0;
				end = 8;
			}
			else if (timeArray.length <position+4)
			{
				/*Means we are at the end of the array and we need to take values from the 
				 * end of the array.
				 */
				start = timeArray.length-8;
				end = timeArray.length;
			}
			else
			{
				/*Means we are somewhere in the middle of the array.
				 * 
				 */
				start = position-4;
				end = position+3;
			}

			for (int g=0; g<neqns; g++)
			{
				double yi=0;
				for (int i=start; i<end; i++)
				{
					product=newStates[i][g];

					for (int j=start; j<end; j++)
					{
						if(i!=j)
						{
							product=product*(currentTime-timeArray[j])/(timeArray[i]-timeArray[j]);
						}
					}
					yi=yi+product;
				}
				interpolatedState[g] = yi;
			}
		}
		return interpolatedState;
	}

	/** ddHermite adapted from Sun Hur-Diaz Matlab Code. 
	 * Forms the Hermite interpolating polynominal by modifying the divided 
	 * difference basis of order neqns and interpolates for the given xo.
	 * @author Stephen Metcalfe
	 * @date September 2010
	 * @param	xo		current user specified output time.			(scalar)
	 * @param	x		array of times from the interpolator.		(N*1)
	 * @param	f		array of states from the interpolator.		(N*1)
	 * @param	fd		array of derivative data.					(N*1)
	 * @param	neqns	number of states.							(scalar)
	 * @return	interpolated state									(scalar)
	 */
	public static double[] ddHermite(double xo, double[] x, double[] f, double[] fd, int neqns)
	{
		int i, j, k, l, di, diMax;
		int iMax = x.length*2;
		int jMax = (neqns*2)+1;
		double p, s, fo, fdo;
		double[][] tab = new double[iMax][jMax];
		double[] result = new double[2];

		// construct the divided difference table
		for (i=0; i<iMax-1; i+=2)
		{
			k = i>>1;
			l = i+1;
			tab[l][0] = tab[i][0] = x[k];
			tab[l][1] = tab[i][1] = f[k];
						tab[i][2] = fd[k];
			for (j=3; j<jMax; j++)
			{
				tab[l][j] = tab[i][j] = 0.0;
			}
		}
		// derive the missing states
		for (i=1; i<iMax-1; i+=2)
		{
			tab[i][2] = (tab[i+1][1] - tab[i-1][1])/(tab[i+1][0] - tab[i-1][0]);
		}
		// derive the coefficients
		for (j=3; j<jMax; j++)
		{
			di = j-1;
			diMax = iMax-j+1;
			for (i=0; i<diMax; i++)
			{
				tab[i][j] = (tab[i+1][j-1]-tab[i][j-1])/(tab[i+di][0]-tab[i][0]);
			}
		}
		// Interpolate on the function
		double v[] = new double[jMax];
		fo = tab[0][1];
		p  = 1;
		for (i=2; i<jMax; i++)
		{
			l = i-2;
			v[l] = (xo - tab[l][0]);
			p = p * v[l];
			fo = fo + (tab[0][i] * p);
		}
		// Interpolate on the function derivative
		fdo = tab[0][2];
		for (i=3; i<jMax; i++)
		{
			s = 0;
			for (j=0; j<=i-2; j++)
			{
				p = 1;
				for (k=0; k<=i-2; k++)
				{
					if(k!=j)
					{
						p = p*v[k];
					}
				}
				s = s+p;
			}
			fdo = fdo+tab[0][i]*s;
		}
		result[0] = fo;
		result[1] = fdo;
		return result;
	}

	/** HermiteInterpolator
	 * Calculate states based on less than 8 newStates
	 * Uses ddHermite to interpolate based on modified divided difference
	 * @author	Stephen Metcalfe
	 * @date	September 2010
	 * @param	simTime		array of times for the interpolator.			(N*1)
	 * @param	currentTime current user specified output time.				(scalar)
	 * @param	neqns		number of states.								(scalar)
	 * @param	simStates	double array of states (x,y,z,xdot,ydot,zdot).	(N*neqns)
	 * @return	array of interpolated states								(neqns*1)
	 *
	 * NOTE: HermiteInterpolator expects paired state and derivative data elements in the form x,y,z,xdot,ydot,zdot
	 * The number of dimensions must match half the number of states ie [x,y,zdot,ydot] gives an neqns of 4.
	 * 
	 * NOTE: Be sure that the units of the simTime, currentTime, and the 
	 * simStates values and derivatives are consistent.
	 */
	public static double[] HermiteInterpolator(double[] simTime, double currentTime, int neqns, double[][] simStates)
	{
		int g, i, l, d;
		int deriv = neqns / 2;	// index to derivatives, xdot, ydot, zdot
		int entries = simTime.length;
		int start = 0;
		double[] interpolatedState = new double[neqns];

		if (neqns != 2 && neqns != 4 && neqns != 6)
		{
			throw new java.lang.IndexOutOfBoundsException("Not enough equations passed " 
				+ neqns + "\nThe Hermite Interpolator expects paired state and derivative data.");
		}
		else if (simTime.length != simStates.length)
		{
			throw new java.lang.IndexOutOfBoundsException("Length of simTime ("
				+ simTime.length + ") does not match length of simStates ("
				+ simStates.length + ")");
		}

		if (entries > 7)							// if more than 7 time entries
		{	
			for (start=0; start<entries; start++)	// find the first time AFTER currentTime
			{	
				if(currentTime<simTime[start])		// if time in array is after currentTime
				{	
					break;							// leave start as this index
				}
			}

			if (start < 4)							// if currentTime is near the beginning
			{
				start = 0;							//		use the first 7 entries
			}
			else if (entries-start < 3)				// if currentTime is near the end
			{
				start = entries-7;					//		use the last 7 entries
			}
			else									// if currentTime is in the middle
			{
				start = start-4;					//		bracket currentTime
			}
			entries = 7;							// use only 7 entries
		}

		// allocate arrays to hold data once count is known
		double[]  x  = new double[entries];			// time array
		double[]  f  = new double[entries];			// position array
		double[]  fd = new double[entries];			// velocity array

		// extract the subset of times 
		for (l=start, i=0; i<entries; i++, l++)
		{
			x[i] = simTime[l];						// extract time values
		}

		for (d=deriv, g=0; g<deriv; g++, d++)
		{											// extract the subset of states
			for (l=start, i=0; i<entries; i++, l++)
			{
				f[i]  = simStates[l][g];			// extract state value
				fd[i] = simStates[l][d];			// extract derivative value
			}

			// interpolate using divided difference
			double[] partial = ddHermite(currentTime, x, f, fd, entries);
			
			// check for limit conditions being returned
			if (Double.isNaN(partial[0]) || Double.isInfinite(partial[0]) ||
				Double.isNaN(partial[1]) || Double.isInfinite(partial[1]))
			{
				throw new java.lang.IndexOutOfBoundsException("ddHermite returned Not A Number or Infinite");
			}
			interpolatedState[g]       = partial[0];// interpolated state value
			interpolatedState[g+deriv] = partial[1];// interpolated derivative value
		}
		return interpolatedState;
	}
	
}
