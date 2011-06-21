/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2010 United States Government as represented by the
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
 *  Stephen D Metcalfe  06/13/2010   	Original
 *  Stephen D Metcalfe  06/29/2010   	Removed unused imports causing warnings. Added testing 
 *  									at each time entry to verify changes to start selection 
 *  									in HermiteInterpolator per Kevin Berry. Modified test 
 *  									conditions to avoid errors due to double data type limits. 
 *  									Improved instrumentation to identify NaN and Infinity being 
 *  									returned from ddHermite. Refactored to create a function 
 *  									to be called from Matlab. 
 *  Stephen D Metcalfe  06/29/2010   	Removed HermiteInterpolator limit on current time < start time 
 */

package jat.matlabInterface.unittest; 

import junit.framework.TestCase;

import org.junit.Test;

/** The JATIntegratorsTest provides a JUnit TestCase for the functions in the JATIntegrators class.
 * This used to be the ODToolboxIntegratorsTest but the integrators were moved to the JATIntegrators class.
 * 
 * Currently this is incomplete and does not test all of the integrators.
 *  
 * @author Stephen D Metcalfe, Emergent Space Technologies
 */
public class JATIntegratorsTest extends TestCase {
	
	/** TODO Generate a set of generic test points defined by parameters */
	
	/* the test point class contains the parameters used for testing interpolators */
	/* 
	private class testpoint
	{
		public int neqns;		// number of equations (2*number of dimensions)
		public int entries;		// number of entries
		public int entry;		// entry under test
		public int offset;		// offset - 0 = before, 1 = at time, 2 = after time
		public double start;	// start time
		public double step;		// step interval size
		public double time;		// entry time
		public double current; 	// current time
		public double[] A;		// acceleration
		public double[] x;		// time array
		public double[][] f;	// position & velocity array
		public double[] ex;		// expected results
		
		// default constructor sets values for an undefined test
		public testpoint()
		{
			neqns = entries = entry = offset = 0;
			start = step = time = current = 0;
			A = null;
			x = null;
			f = null;
			ex = null;
		}
	}
	*/
	
	/** Generate a generic set of test points defined by the parameters
	 * @param	DIMENSIONS  - array of dimensions
	 * @param	ENTRIES		- array of numbers of entries 
	 * @param	START		- array of start times
	 * @param	STEPS		- array of interval step sizes
	 * @param	P			- initial position
	 * @param	V			- initial velocity
	 * @param	AN			- acceleration when more than 1 entry
	 * 
	 * NOTE no acceleration is used when a test point only has 1 entry
	 */
	/*
	private static ArrayList < testpoint > gentests(int[] DIMENSIONS, int[] ENTRIES, double[] START, double[] STEPS,
										double[] P, double[] V, double[] AN) {

		int d, dims, neqns, e, entries, s, st, te, p, o, numtests = 0, testnum = 0;
		double start, step, t, offset, c;

		// validate the acceleration parameter - must be enough for the number of dimensions		
		for (d = 0; d < DIMENSIONS.length; d++)
		{
			if (DIMENSIONS[d] > AN.length)
			{
				throw new java.lang.IndexOutOfBoundsException("Insufficient accelleration parameters for dimensions.");
			}
		}
		
		// find how many test points will be generated
		d = DIMENSIONS.length * ENTRIES.length * START.length * STEPS.length * 3;
		for(e = 0; e < ENTRIES.length; numtests += d * ENTRIES[e++]);

		// allocate space for the tests
		testpoint   test;
		ArrayList < testpoint > testset = new ArrayList;
		
		for (d = 0; d < DIMENSIONS.length; d++)								// iterate the dimensions
		{
			dims = DIMENSIONS[d];											// select the number of dimensions
			neqns = dims*2; 												// set the number of equations
			double[] ex = new double[neqns];								// expected results
			double[] A  = new double[dims];
			test.neqns	= neqns;											// number of equations

			for (e = 0; e < ENTRIES.length; e++)							// iterate the number of entries
			{
				entries = ENTRIES[e];										// select the number of entries
				if (entries == 1)
					for (p = 0; p < dims; A[p++] = 0.0);					// use no acceleration if only 1 entry
				else
					A = AN;													// otherwise use acceleration defined

				double[] x = new double[entries];							// allocate time array
				double[][] f = new double[entries][neqns];					// allocate data array
				
				for (s = 0; s < START.length; s++)							// iterate start times
				{
					start = START[s];										// select starting time
					
					for (st = 0; st < STEPS.length; st++)					// iterate time interval sizes
					{
						step = STEPS[st];									// select time interval size

						double OFFSETS[] = {-(step/4), 0, step/4};			// offsets before, at and after entry times
						
						for (te = 0; te < entries; te++)					// for each entry in the test data
						{
							t = start + (te * step);						// calculate time at each interval
							x[te] = t;										// set the entry time
							
							for (p = 0; p < dims; p++)						// for each dimension
							{	
								f[te][p] = P[p]+(V[p]*t)+(0.5*A[p]*(t*t));	// calculate position
								f[te][p+dims] = V[p]+(A[p]*t);				// calculate velocity 
							}
						}
						
						
						
						for (te = 0; te < entries; te++)					// for each entry in the test data
						{
							for (o = 0; o < 3; o++)							// for offsets before, at and after
							{
								t = x[te]; 
								offset = OFFSETS[o];						// get offset size 
								c = t+offset;								// calculate time with offset

								for (p = 0; p < dims; p++)					// for each dimension
								{
									ex[p] = P[p]+(V[p]*c)+(0.5*A[p]*(c*c));	// calculate expected position
									ex[p+dims] = V[p]+(A[p]*c);				// calculate expected velocity 
								}
								
								testset[testnum++].neqns	= neqns;		// number of equations
								testset[testnum++].entries 	= entries;		// number of entries
								testset[testnum++].entry	= te;			// entry under test
								testset[testnum++].start	= start;		// start time
								testset[testnum++].step		= step;			// step interval size
								testset[testnum++].time     = t;			// entry time
								testset[testnum++].offset	= o;			// before, at time, after
								testset[testnum++].current	= c; 			// current time
								testset[testnum++].A		= A;			// acceleration
								testset[testnum++].x		= x;			// time array
								testset[testnum++].f		= f;			// position/velocity array
								testset[testnum++].ex		= ex;			// expected results
							} 
						}
					}
				}
			}
		}
		return testset;
	}
	*/
	
	/** LinearInterpolator_regression
	 *  TODO Test the LinearInterpolator using a range of conditions
	 * @return	double array of results 								(N*N)
	 * 			each row consist of
	 *				number of tests expected
	 *				number of tests performed
	 *				number of dimensions
	 *				number of entries
	 *				start time
	 *				interval step
	 *				current time
	 *				results of test
	 *
	public static double[][] LinearInterpolator_regression() {
	}
	*/
	
	/** TODO Perform the LinearInterpolator_regression test and check for errors returned.
	 *
	@Test
	public void testLinearInterpolator() {
	}
	*/
	
	/** LagrangianInterpolator_regression
	 *  TODO Test the LagrangianInterpolator using a range of conditions
	 * @return	double array of results 								(N*N)
	 * 			each row consist of
	 *				number of tests expected
	 *				number of tests performed
	 *				number of dimensions
	 *				number of entries
	 *				start time
	 *				interval step
	 *				current time
	 *				results of test
	 *
	public static double[][] LinearInterpolator_regression() {
	}
	*/
	
	/** TODO Perform the LagrangianInterpolator_regression test and check for errors returned.
	 *
	@Test
	public void testLagrangianInterpolator() {
	}
	*/
	
	@Test
	/** Perform the HermiteInterpolator_regression test and check for errors returned.
	 */
	public void testHermiteInterpolator() {
		
		final int TODATA  = 7;	// index to first result element
	
		double[][] rtn = new double[1][1];
		try
		{
			rtn = HermiteInterpolator_regression.run_regression();
		}
		finally
		{
			// test for exception
			assertFalse(rtn.length == 0);
			
			// check each result test
			for (int t = 0; t < rtn.length; t++)
			{
				// first value contains the number of tests expected
				int max = (int)rtn[t][0];

				// first value contains the number of tests performed
				int done = (int)rtn[t][1];

				// test for exception causing missed tests
				assertFalse(max != done);
				
	 			// check each test result for failure
				for (int r = TODATA; r < TODATA+max; r++)
				{
					assertTrue(rtn[t][r] == 0);
				}
			}
		}
	};
	
};