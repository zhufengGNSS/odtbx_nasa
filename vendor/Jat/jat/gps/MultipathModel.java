package jat.gps;

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

import jat.matvec.data.*;
import jat.math.*;
import jat.cm.*;
import jat.alg.*;

/**
* <P>
* The MultipathModel.java Class provides a statistic model of the 
* multipath effects on GPS carrier phase and C/A code measurements 
* near the ISS.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class MultipathModel implements ScalarFunction {

	/** Radar cross-sectional area of the ISS (or reflector) */
	private double a_rcs = 5000.0;
	
	/** Number of multipath rays */
	private int n = 10;
	
	/** GPS L1 frequency */
	private static final double freq = 1575.42E+06;
	
	/** GPS L1 wavelength in meters */
	private static final double lambda = Constants.c / freq;
	
	/** GPS L1 wavelength squared */
	private static final double lamsq = lambda * lambda;
		
	/** 2 PI */
	private static final double twopi = 2.0 * MathUtils.PI;
	
	/** 4 PI */
	private static final double fourpi = 4.0 * MathUtils.PI;
	
	/** (4 PI)^2 */
	private static final double fourpi2 = fourpi * fourpi;
	
	/** (4 PI)^3 */
	private static final double fourpi3 = fourpi2 * fourpi;

	/** C/A code chipping rate in Hz */
	private static final double chipRate = 1.023E+06;
	
	/** C/A code chipping period in seconds */
	private static final double T = 1.0 / chipRate;
	
	/** Time advance/delay of early/late code */
	private static final double tau_d = 0.5 * T;
	
	private long seed = -1;

	/** Random numbers */
	private RandomNumber tk = new RandomNumber(seed);
	private RandomNumber bk = new RandomNumber(2 * seed);
	private RandomNumber pk = new RandomNumber(3 * seed);

	/** ZeroFinder object */
	private ZeroFinder rf = new ZeroFinder(this, 500, 1.0E-12, 1.0E-12);

	/** Array containing multipath relative time delays */
	private double[] tau;
	
	/** Array containing multipath amplitudes */
	public double[] beta;
	
	/** Array containing the the multipath amplitudes squared */
	private double[] betasq;
	
	/** Array containing the multipath relative phase angles in radians */
	private double[] psi;
	
	/** Direct signal amplitude */
	public double beta0 = 0.0;

	/** Constructor
	 * @param area radar cross-sectional area of ISS (or reflector)
	 * @param nrays number of multipath rays to simulate
	 */
	public MultipathModel(double area, int nrays) {
		this.a_rcs = area;
		this.n = nrays;
		tau = new double[this.n];
		beta = new double[this.n];
		betasq = new double[this.n];
		psi = new double[this.n];
	}
	
	/**
	 * Set the seed
	 * @param sd long containing seed to be use
	 */
	public void setSeed(long sd){
		this.seed = sd;
	}

	/** Returns the antenna gain. Currently using a cardiod pattern.
	 * Any antenna gain pattern function could be substituted here.
	 * @param theta angle between the direct line of sight 
	 * and the antenna boresight in radians.
	 * @return the antenna gain
	 */
	private double antennaGain(double theta) {
		double gain = Math.abs(Math.cos(theta));
		return gain;
	}

	/** PRN code correlation function
	 * @param tau time delay
	 * @return the output of the correlator
	 */
	private double correlation(double tau) {
		double R = 0.0;
		double abstau = Math.abs(tau);
		if (abstau <= T) {
			R = 1.0 - abstau / T;
		}
		return R;
	}
	
	/** Compute the multipath environment based on the input value.
	 * Fills the tau, beta and psi arrays.
	 * @param prn GPS SV PRN number
	 * @param dr distance from the ISS (or reflector)
	 * @param theta angle between direct line of sight and antenna boresight
	 * @param rGPS vector containing the ECI position vector of the GPS SV
	 * @param rISS vector containing the ECI position vector of the ISS
	 * @param rSTS vector containing the ECI position vector of the spacecraft
	 */
	public void environment(int prn, double dr, double theta, VectorN rGPS,
		VectorN rISS, VectorN rSTS) {
			
		// compute the distance from the GPS satellite to the ISS
		VectorN rGPSISS = rGPS.minus(rISS);
		double rgi = rGPSISS.mag();
		double rgi2 = rgi * rgi;
		double dr2 = dr * dr;

		// generate some seeds for random number generators
		long seed1 = (long) prn;
		long seed2 = (long) prn - 1;
		long seed3 = (long) prn - 2;

		// compute the distance from the GPS satellite to the receiver
		VectorN rGPSSTS = rGPS.minus(rSTS);
		double rgs = rGPSSTS.mag();
		double rgs2 = rgs * rgs;

		// compute beta0
		double antgain = antennaGain(theta);
		this.beta0 = Math.sqrt((lamsq * antgain) / (fourpi2 * rgs2));
		//		System.out.println(beta0+"\t0.0");

		// compute tau_bar
		double tau_bar = Math.abs(dr) / Constants.c;

		// loop thru the n multipath rays
		for (int k = 0; k < this.n; k++) {
			this.tau[k] = tk.exponential(tau_bar);

			// Compute bbarsq
			double decay = Math.exp(-tau[k] / tau_bar);
			double numerator = a_rcs * lamsq * decay;
			double denominator = fourpi3 * dr2 * rgi2;
			double bbarsq = numerator / denominator;

			// get beta[k]
			this.betasq[k] = bk.exponential(bbarsq);
			this.beta[k] = Math.sqrt(betasq[k]);

			// get psi[k]
			this.psi[k] = pk.uniform(0.0, twopi);
		}
	}

	/** Compute the range errors in carrier phase measurements.
	 * Must be called only after environment() has been called.
	 * @return range error in carrier phase measurement in meters.
	 */
	public double carrierPhaseError() {
		double sinSum = 0.0;
		double cosSum = 0.0;

		// loop thru the n multipath rays
		for (int k = 0; k < this.n; k++) {

			sinSum = sinSum + beta[k] * Math.sin(psi[k]);
			cosSum = cosSum + beta[k] * Math.cos(psi[k]);
		}

		double dphi = Math.atan2(sinSum, (beta0 + cosSum));
		double error = lambda * dphi / twopi;
		return error;
	}

	/** Compute the range errors in C/A code measurements.
	 * Must be called only after environment() has been called.
	 * @return range error in C/A code measurement in meters.
	 */
	public double pseudorangeError() {
		
		// solve D(tau) = 0 for tau using the secant method
		double tau = this.rf.secant(0.0, 0.01 * T);
//		double tau = this.rf.ridder(-1.499*this.T, 1.499*this.T);

		// convert to range error in meters
		double out = -1.0 * tau * Constants.c;
		return out;
	}

	/** Discriminator function to be solved for tau.
	 * This function is called by the ZeroFinder (secant method) to
	 * determine the tracking loop error, which is converted to 
	 * a range error for C/A code measurements.
	 * @param tau
	 * @return the output of the discriminator function
	 */
	public double evaluate(double tau) {
		//		double tau = tau_in/constants.c;
		double d = 0.0;
		double rplus = this.correlation(tau + tau_d);
		double rminus = this.correlation(tau - tau_d);
		double term1 = rplus * rplus - rminus * rminus;

		double term2 = 0.0;
		double term3 = 0.0;

		// form the sum over i
		for (int i = 0; i < this.n; i++) {
			double rplus_i = this.correlation(tau + tau_d + this.tau[i]);
			double rminus_i = this.correlation(tau - tau_d + this.tau[i]);
			double cospsi = Math.cos(this.psi[i]);
			double beta_i = this.beta[i] / this.beta0;
			term2 = term2+2.0*beta_i*cospsi*(rplus*rplus_i-rminus*rminus_i);

			// form the sum over j
			for (int j = 0; j < this.n; j++) {
				double rplus_j =
					this.correlation(tau + tau_d + this.tau[j]);
				double rminus_j =
					this.correlation(tau - tau_d + this.tau[j]);
				double cosdpsi = Math.cos(this.psi[i] - this.psi[j]);
				double beta_j = this.beta[j] / this.beta0;
				term3 = term3+beta_i*beta_j*cosdpsi*(rplus_i*rplus_j-rminus_i*rminus_j);
			}
		}

		// sum up the pieces
		d = term1 + term2 + term3;
		return d;
	}
}
