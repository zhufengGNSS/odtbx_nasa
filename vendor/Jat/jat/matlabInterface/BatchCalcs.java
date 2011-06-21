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

package jat.matlabInterface;

import jat.matvec.data.Matrix;
import jat.spacetime.EarthRef;
import jat.spacetime.Time;

/**
 * The Matlab-Java interface does impose some overhead, and Matlab code 
 * with a tight iteration over calls to JAT functions can be very
 * inefficient.  This class provides a single JAT method call to perform
 * multiple iterations of common functions.
 * 
 * @author rantonucci
 *
 */
public class BatchCalcs {

	/**
	 * This determines the tranformation matrices for doing ECI to ECEF transformation 
	 * at multiple times
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @param todExpireDays for efficiency, the true of date transformation (precession
	 * and nutation) does not need to be recomputed every time.  todExpireDays indicates 
	 * for how long a single precession and nutation can be reused.  Values are in days.
	 * 0 indicates to always recompute the precession and nutation.
	 * @return a 3x3xN array where each 3x3 array is the transformation matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] eciToEcefXform(double[] mjd_utc, double todExpireDays) {
		return eciEcefXform(mjd_utc, todExpireDays, false);
	}		
	
	/**
	 * This determines the tranformation matrices for doing ECEF to ECI transformation 
	 * at multiple times
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @param todExpireDays for efficiency, the true of date transformation (precession
	 * and nutation) does not need to be recomputed every time.  todExpireDays indicates 
	 * for how long a single precession and nutation can be reused.  Values are in days.
	 * 0 indicates to always recompute the precession and nutation.
	 * @return a 3x3xN array where each 3x3 array is the transformation matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] ecefToEciXform(double[] mjd_utc, double todExpireDays) {
		return eciEcefXform(mjd_utc, todExpireDays, true);
	}		

	/**
	 * Since the ECI-to-ECEF transform is just the transpose of the ECEF-to-ECI transform, both
	 * methods are handled by a single method with a argument of whether to transpose the result
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @param todExpireDays for how long a single precession and nutation can be reused (in days).
	 * @param transpose true for ECEF-to-ECI, false for ECI-to-ECEF
	 * @return a 3x3xN array of transformation matrices
	 */
	private static double[][][] eciEcefXform(double[] mjd_utc, double todExpireDays, boolean transpose) {
			
		// Preallocate the array
		int numSamples = mjd_utc.length;
		double[][][] out = new double[3][][];
		for(int rowctr=0; rowctr<3; ++rowctr) {
			out[rowctr] = new double[3][];
			for(int colctr=0; colctr<3; ++colctr) {
				out[rowctr][colctr] = new double[numSamples];
			}
		}
		
		
		double refTime = -Integer.MIN_VALUE;
		EarthRef ref = null;
		Matrix xform = null;
		for(int timeCtr=0; timeCtr<numSamples; ++timeCtr) {
			// Get the transformation matrix.
			// Either create a new Earth reference frame which computes a new 
			// precession and nutation, or use the existing Earth reference frame
			// which reuses precession and nutation.
			double nextTime = mjd_utc[timeCtr];
			if (nextTime - refTime > todExpireDays) {
				ref = new EarthRef(nextTime);
				refTime = nextTime;
				xform = ref.ECI2ECEF();
			}
			else {
				Time t = new Time(nextTime);
				xform = ref.eci2ecef(t.mjd_ut1(), t.mjd_tt(), false);
			}
			
			// Now copy the matrix into the array.
			double[][] A = xform.getArray();
			for(int rowctr=0; rowctr<3; ++rowctr) {
				for(int colctr=0; colctr<3; ++colctr) {
					out[rowctr][colctr][timeCtr] = (transpose ? A[colctr][rowctr] : A[rowctr][colctr]);
				}
			}			
		}
		
		return out;
	}
	
	/**
	 * This determines the tranformation matrices for doing J2000 to True of Date transformation 
	 * at multiple times
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the transformation matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] j2000ToTODXform(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				return ref.trueOfDate(new Time(t_mjd_utc));
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		
	
	/**
	 * This gets Earth's precession matrices at multiple time samples
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] getPrecession(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				return ref.PrecMatrix(new Time(t_mjd_utc).mjd_tt());
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		
	
	/**
	 * This gets Earth's nutation matrices at multiple time samples
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] getNutation(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				return ref.NutMatrix(new Time(t_mjd_utc).mjd_tt());
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		
	
	/**
	 * This gets Earth's Greenwich Hour Angle matrices at multiple time samples
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] getGHAMatrix(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				Time t = new Time(t_mjd_utc);
				return ref.GHAMatrix(t.mjd_ut1(), t.mjd_tt());
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		
	
	/**
	 * This gets Earth's pole matrices at multiple time samples
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] getPoleMatrix(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				EarthRef newref = new EarthRef(t_mjd_utc);
				return newref.PoleMatrix();
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		

	/**
	 * This gets Earth's ecliptic matrices at multiple time samples
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	public static double[][][] getEcliptic(double[] mjd_utc) {
		GetMatrixFunc func = new GetMatrixFunc() {
			public Matrix get(EarthRef ref, double t_mjd_utc) {
				return ref.EclMatrix(new Time(t_mjd_utc).mjd_tt());
			}
		};
		return getEarthRefMatrices(mjd_utc, func);
	}		

	/**
	 * This is the interface of a function object, which we use so we can combine
	 * the code for all the different matrix request objects and reuse a single
	 * piece of code.
	 * Anonymous subclasses of this function object are created to call the specific
	 * method to return the desired matrix.
	 */
	private static interface GetMatrixFunc {
		public Matrix get(EarthRef ref, double mjd_utc);
	}

	/**
	 * This is is the main code for all the "getMatrix" functions (except for the ECI2ECEF and
	 * ECEF2ECI which need to do special caching).  This common code preallocates the array and
	 * copies the desired matrix for each time sample into the array, but it relies on a function
	 * object to actually request the desired matrix from the EarthRef.
	 * @param mjd_utc a 1xN array of UTC times in MJD format
	 * @param func a function object that defines a get() method that calls the right
	 * method on EarthRef to get the desired matrix
	 * @return a 3x3xN array where each 3x3 array is the matrix for the
	 * corresponding time in mjd_utc.
	 */
	private static double[][][] getEarthRefMatrices(double[] mjd_utc, GetMatrixFunc func) {
		
		// Preallocate the array
		int numSamples = mjd_utc.length;
		double[][][] out = new double[3][][];
		for(int rowctr=0; rowctr<3; ++rowctr) {
			out[rowctr] = new double[3][];
			for(int colctr=0; colctr<3; ++colctr) {
				out[rowctr][colctr] = new double[numSamples];
			}
		}
		
		// Now get the matrices.  Use a function object to actually
		// invoke the right method to return the desired matrix.
		EarthRef ref = new EarthRef(mjd_utc[0]);
		Matrix m = null;
		for(int timeCtr=0; timeCtr<numSamples; ++timeCtr) {
			m = func.get(ref, mjd_utc[timeCtr]);
			
			// Now copy the matrix into the array.
			double[][] A = m.getArray();
			for(int rowctr=0; rowctr<3; ++rowctr) {
				for(int colctr=0; colctr<3; ++colctr) {
					out[rowctr][colctr][timeCtr] = A[rowctr][colctr];
				}
			}			
		}
		
		return out;
	}
}
