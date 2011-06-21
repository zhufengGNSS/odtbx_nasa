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

package jat.matvec.data;
import java.io.*;

/**
 * <P>
 * The RandomNumber Class is a random number generator used to create samples of uniform
 * and Gaussian distributions. The algorithms come from Chapter 7 of Numerical Recipes.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class RandomNumber {

	protected long idum;

	/** Creates a new instance of RandomNumber */
	public RandomNumber() {
		idum = -1;
	}

	/** Create a RandomNumber generator given the dimension and seed.
	 * @param seed Seed for the random number generator.
	 */
	public RandomNumber(long seed) {
		if (seed > 0)
			seed = -seed;
		idum = seed;
	}

	private static final int IM1 = 2147483563;
	private static final int IM2 = 2147483399;
	private static final int IMM1 = IM1 - 1;
	private static final int IA1 = 40014;
	private static final int IA2 = 40692;
	private static final int IQ1 = 53668;
	private static final int IQ2 = 52774;
	private static final int IR1 = 12211;
	private static final int IR2 = 3791;
	private static final int NTAB = 32;
	private static final long NDIV = 1 + IMM1 / NTAB;
	private static final double AM = 1.0 / IM1;
	private static final double EPS = 1.0E-07;
	private static final double RNMX = 1.0 - EPS;

	private long[] iv = new long[NTAB];
	private long iy = 0;
	private long idum2 = 123456789;

	/** Returns a uniform random deviate between 0.0 and 1.0, exclusive of the endpoints. 
	 * Call with a negative integer to initialize and do not alter iseed between successive calls. 
	 * Based on ran2.c from Numerical Recipes.
	 * @param iseed Seed for random number generator, initialize with iseed < 0.
	 * @return returns a uniform random deviate.
	 */
	public double uniformDeviate() {
		double temp;
		int j;
		long k;
		if (idum <= 0) {
			idum = Math.max(-idum, 1);
			idum2 = idum;
			for (j = NTAB + 7; j >= 0; j--) {
				k = idum / IQ1;
				idum = IA1 * (idum - k * IQ1) - k * IR1;
				if (idum < 0)
					idum = idum + IM1;
				if (j < NTAB)
					iv[j] = idum;
			}
			iy = iv[0];
		}
		k = idum / IQ1;
		idum = IA1 * (idum - k * IQ1) - k * IR1;
		if (idum < 0)
			idum = idum + IM1;
		k = idum2 / IQ2;
		idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
		if (idum2 < 0)
			idum2 = idum2 + IM2;
		j = (int) (iy / NDIV);
		iy = iv[j] - idum2;
		iv[j] = idum;
		if (iy < 1)
			iy = iy + IMM1;
		temp = Math.min(AM * iy, RNMX);
		return temp;
	}

	/** Generate a random number from a uniform random variable.
	 @param min    Min of the random variable.
	 @param max    Max of the random variable.
	 @return      A double.
	 */

	public double uniform(double min, double max) {
		double x = min + (max - min) * uniformDeviate();
		return x;
	}

	private int iset = 0;
	private double gset;

	/** Returns a normally distributed deviate with zero mean and unit variance using uniform(iseed) as the source of uniform deviates. 
	 * Based on gasdev.c from Numerical Recipes. Uses a Box-Muller transformation.
	 * @param iseed Seed for random number generator.
	 * @return returns a sample from normally distribution.
	 */
	public double normal() {
		double fac;
		double rsq;
		double v1;
		double v2;

		if (idum < 0)
			iset = 0;
		if (iset == 0) {
			do {
				v1 = 2.0 * uniformDeviate() - 1.0;
				v2 = 2.0 * uniformDeviate() - 1.0;
				rsq = v1 * v1 + v2 * v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac = Math.sqrt(-2.0 * Math.log(rsq) / rsq);
			gset = v1 * fac;
			iset = 1;
			return v2 * fac;
		} else {
			iset = 0;
			return gset;
		}
	}

	/** Generate a random number from an normal random variable.
	@param mean mean.
	@param sigma sigma.
	@return      A double.
	*/
	public double normal(double mean, double sigma) {
		double x = this.normal() * sigma + mean;
		return x;
	}

	/** Generate a random number from an exponantial random variable (Mean = 1/lambda, variance = 1/lambda^2).
	@param lambda    Parmaeter of the exponential random variable.
	@return      A double.
	*/
	public double exponential(double beta) {
		double x = -beta * Math.log(uniformDeviate());
		return x;
	}

	public static void main(java.lang.String args[]) throws IOException {

		FileOutputStream outf = new FileOutputStream("c:\\temp\\random.txt");
		PrintWriter pw = new PrintWriter(outf);
		RandomNumber x = new RandomNumber();
		double beta = 3.30E-07;
		for (int i = 0; i < 3000; i++) {
			pw.println(
				x.uniformDeviate()
					+ "\t"
					+ x.normal()
					+ "\t"
					+ x.exponential(beta));
			System.out.println(
				x.uniformDeviate()
					+ "\t"
					+ x.normal()
					+ "\t"
					+ x.exponential(beta));
		}
		pw.close();
		outf.close();
	}

}
