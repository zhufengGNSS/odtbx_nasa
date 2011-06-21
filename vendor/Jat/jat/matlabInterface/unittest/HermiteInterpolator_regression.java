/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2011 United States Government as represented by the
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
 */

package jat.matlabInterface.unittest;

import jat.matlabInterface.JATIntegrators;

/**
 * A standard Java class that contains the significant data for testing the 
 * HermiteInterpolator.  This class is extracted from the JUnit TestCase,
 * JATIntegratorsTest, so that it can be run without the JUnit framework.
 * 
 * @author abrown
 */
public class HermiteInterpolator_regression {

	/** 
	 *  Test the HermiteInterpolator over 1, 2 and 3 dimensions using a range of data set sizes with different 
	 *  starting times, interval sizes. Sample at the beginning, mid and end points of each data set with
	 *  offsets to before, at and after each sample point by 1/4 the time interval.  
	 *  Results are compared as a percentage of the expected result against an allowed percentage error. 
	 *  NOTE: This test limits acceleration when there is only 1 entry in the data set.
	 *  TODO remove the loops and use a generic set of test points   
	 * @author	Stephen Metcalfe
	 * @date	9/29/2010
	 * @return	double array of results 								(N*N)
	 * 			each row consist of
	 *				number of tests expected
	 *				number of tests performed
	 *				third value is number of dimensions
	 *				number of entries
	 *				start time
	 *				interval step
	 *				current time
	 *				results of test
	 */
	public static double[][] run_regression() {
		/* configure test */ 
		final int[]		DIMENSIONS	= { 1, 2, 3};									/** Dimensions to test - min 1 max 3 */
		final int[] 	ENTRIES 	= {	20, 15, 10, 8, 7, 6, 5, 4, 3, 2, 1 }; 		/** Number of entries sequence */
		final double[] 	START 		= {	99999.0, 999.0, 0.0, -999.0, -99999.0 };	/** Starting times sequence */
		final double[] 	STEPS 		= {	1, 999.0, 9999.0 };							/** Time interval sequence - NO NEGATIVES */
		final double[] 	P 			= {	1.0, 2000.0, -3000000.0 };					/** Initial positions */
		final double[] 	V 			= {	4.0, 5000.0, 6000000.0 };					/** Initial velocities */
		final double[] 	AN 			= {	7.0, 800.0, 90000.0 };						/** Acceleration with more that 1 entry */
		final double[]	A1 			= { 0.0, 0.0,    0.0 };							/** Acceleration with only 1 entry */
		final double	POSERR		= 0.000000001;									/** Position error allowed */
		final double	VELERR		= 0.000000001;									/** Velocity error allowed */

		final String VALUE[] = {"x   ","y   ","z   ","xdot","ydot","zdot"}; 		// dimension names
		final String BEFOREATAFTER[] = { "Before", "At", "After" };					// offset names

		/* local variables */
		double err, percent, max;													// error values
		int rtn_t, rtn_r, d, dims, neqns, e, entries, s, st, te, p, o, r, test = 1, failcount = 0;
		double start, step, time, offset, t, c, sum[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		String status, condition;
		boolean failed = false;
		double[] result;

		// return array
		final int TODATA  = 7;	// index to first result element
		final int MAXTEST = 6;	// maximum number of tests in a set
		double[][] rtn = new double[DIMENSIONS.length*ENTRIES.length*START.length*STEPS.length*20*3][TODATA+MAXTEST];
		for (rtn_t = 0; rtn_t < rtn.length; rtn_t++)
		{
			for (r = 0; r < MAXTEST+TODATA; rtn[rtn_t][r++] = 0);
		}

		for (rtn_t = 0, d = 0; d < DIMENSIONS.length; d++)							// iterate the dimensions
		{
			dims = DIMENSIONS[d];													// select the number of dimensions

			neqns = dims*2; 														// set the number of equations
			double[] ex = new double[neqns];										// expected results

			for (e = 0; e < ENTRIES.length; e++)									// iterate the number of entries
			{
				entries = ENTRIES[e];												// select the number of entries

				double[] A = (entries == 1) ? A1 : AN;								// acceleration for number of entries 
				double[] x = new double[entries];									// allocate time array
				double[][] f = new double[entries][neqns];							// allocate data array

				for (s = 0; s < START.length; s++)									// iterate start times
				{
					start = START[s];												// select starting time
					for (r = 0; r < 6; sum[r++] = 0.0);								// clear error stats

					for (st = 0; st < STEPS.length; st++)							// iterate time interval sizes
					{
						step = STEPS[st];											// select time interval size

						double OFFSETS[] = {-(step/4), 0, step/4};					// offsets before, at and after entry times

						for (te = 0; te < entries; te++)							// generate the test data 
						{
							t = start + (te * step);								// calculate time at each interval
							x[te] = t;												// set the entry time

							for (p = 0; p < dims; p++)								// for each dimension
							{	
								f[te][p] = P[p]+(V[p]*t)+(0.5*A[p]*(t*t));			// calculate position
								f[te][p+dims] = V[p]+(A[p]*t);						// calculate velocity 
							}
						}

						for (te = 0; te < entries; te++)							// for each entry in the test data
						{
							for (o = 0; o < 3; o++, rtn_t++)						// for offsets before, at and after
							{
								time = x[te]; 
								offset = OFFSETS[o];								// get offset size 
								c = time+offset;									// calculate time with offset

								rtn_r = TODATA;										// start the next test set
								rtn[rtn_t][0] = neqns;								// number of tests expected
								rtn[rtn_t][1] = 0;									// number of tests performed
								rtn[rtn_t][2] = dims;								// number of dimensions
								rtn[rtn_t][3] = entries;							// number of entries
								rtn[rtn_t][4] = start;								// start time
								rtn[rtn_t][5] = step;								// interval step
								rtn[rtn_t][6] = c;									// current time

								for (p = 0; p < dims; p++)							// for each dimension
								{
									ex[p] = P[p]+(V[p]*c)+(0.5*A[p]*(c*c));			// calculate expected position
									ex[p+dims] = V[p]+(A[p]*c);						// calculate expected velocity 
								}

								status	  =	"Test   " + test;
								condition =	  "\tDimensions    " + dims +
								"\n\tStart time    " + start + 
								"\n\tInterval size " + step +
								"\n\tTest entry    " + te +  " of " + entries +
								"\n\tEntry time    " + time +
								"\n\tOffset        " + BEFOREATAFTER[o] + " time by " + offset +
								"\n\tTest time     " + c + "\n";

								for (r = 0; r < neqns; rtn[rtn_t][rtn_r+r++] = 1);	// assume test will fail
								failed = true;
								try
								{													// run the HermiteInterpolator
									result = JATIntegrators.HermiteInterpolator(x,c,neqns,f);
									failed = false;
									for (r = 0; r < neqns; r++)						// examine each result
									{
										status = "Test   " + (test + r);
										max = ( r < dims ) ? POSERR : VELERR;		// get maximum error limit
										err = Math.abs(ex[r]-result[r]);			// calculate absolute error
										percent = (err/ex[r])*100;					// calculate percentage error
										// check for errors
										if (Double.isNaN(result[r]) || Double.isInfinite(result[r]) || (percent > max))
										{
											failed = true;							// record the failure
											failcount++;
										}
										else
										{
											rtn[rtn_t][rtn_r+r] = 0;				// record the success
										}

										status += "\t" + VALUE[( r < dims ) ? r : (r - dims) + 3];
										status += (failed ? " failed" : " passed");
										condition += "\t" + VALUE[( r < dims ) ? r : (r - dims) + 3] + 
										(failed ? " failed" : " passed") +
										"\n\t\tExpected      " + ex[r] +
										"\n\t\tReceived      " + result[r] +
										"\n\t\tError         " + err +
										"\n\t\tPercent error " + percent +
										"\n\t\tPercent limit " + max +
										"\n\t\tAverage error " + (sum[r]/test) +
										"\n";
										sum[r] += percent;								// sum percentage errors
									}
									rtn_r += neqns;
									rtn[rtn_t][1] = rtn_r-TODATA;
								}
								finally
								{
									//									System.out.print(status + "\n");
									if (failed == true)
									{
										System.out.print(status + "\n" + condition);
									}
									test += neqns;
								}
							}
						}
					}
				}
			}
		}
		return rtn;
	}

}
