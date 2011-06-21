/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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
* 
*/
package jat.ground_tracking;
import jat.constants.*;
import jat.gps.*;
import jat.math.*;
import jat.matvec.data.*;

/**
 * Tropospheric Refraction Model
 * Ref: Montenbruck, pp. 221 - 224.
 * "A general and accurate tropospheric refraction model is the 
 * Hopfield model, modified by Goad to use the Saastamoinen zenith
 * range correction. It is applicable both to radar data as well as
 * to optical observations."
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 *
 */
public class TropoModel {
	
	private double Ce;
	
	private double Crho;
	
	private static final double term2 = 78.8828 / 77.624;
	private static final double Re = IERS_1996.R_Earth;
	
/**
 * Constructor
 * @param lambda double containing wavelength in meters
 */
	public TropoModel(double lambda){
		double lp2_inv = 1.0 / ((lambda * 1.0E+06)*(lambda * 1.0E+06));
		double denom = (173.3 - lp2_inv);
		double term1 = 170.2649 / denom;
		Ce = term1 * term2;
		double term3 = (173.3 + lp2_inv)/denom;
		Crho = Ce * term3;		
	}

	/**
	 * Compute refraction corrections.
	 * @param p double containing pressure in hPa
	 * @param T double containing temperature in deg K
	 * @param fh double containing relative humidity (0 <= fh <= 1)
	 * @param E double containing elevation angle in radians
	 * @param rho double containing range in m
	 * @return double[] containing tropospheric refraction corrections for range (m)
	 * and elevation (arcsec) measurements
	 */
	public double [] corrections(double p, double T, double fh, double E, double rho){
		// refractivities
		double[] N = new double[2];
		// compute dry component refractivity
		N[0] = 77.624 * p / T;
		// compute wet component refractivity
		double Tc = T - 273.15;
		double e = 6.10 * fh * Math.exp(17.15 * Tc /(234.7 + Tc));
		N[1] = 371900.0 * e / (T*T) - 12.92 * e/T;
		// troposphere heights
		double[] h = new double[2];
		// compute dry troposphere height
		h[0] = 5.0 * 0.002277 * p / (N[0] * 1.0E-06);
		// compute wet troposphere height
		h[1] = 5.0 * 0.002277 * e * (1255.0/T + 0.05) / (N[1] * 1.0E-06);
		// distances to top of the troposphere
		double [] r = new double [2];
		double [][] alpha = new double[9][2];
		double [][] beta = new double[7][2];
		double cosE = Math.cos(E);
		double cosE2 = cosE * cosE;
		double sinE = Math.sin(E);
		for (int j = 0; j < 2; j++){
			r[j] = Math.sqrt((Re + h[j])*(Re + h[j]) - (Re*Re*cosE2)) - Re*sinE;
			double aj = -1.0 * sinE/h[j];
			double bj = - 1.0 * cosE2/(2.0 * h[j] * Re);
			alpha[0][j] = 1.0;
			alpha[1][j] = 4.0*aj;
			alpha[2][j] = 6.0*aj*aj + 4.0*bj;
			alpha[3][j] = 4.0*aj*(aj*aj + 3.0*bj);
			alpha[4][j] = Math.pow(aj, 4) + 12.0*aj*aj*bj + 6.0*bj*bj;
			alpha[5][j] = 4.0*aj*bj*(aj*aj + 3.0*bj);
			alpha[6][j] = bj*bj*(6.0*aj*aj + 4.0*bj);
			alpha[7][j] = 4.0 * aj * bj*bj*bj;
			alpha[8][j] = Math.pow(bj,4);
			beta[0][j] = 1.0;
			beta[1][j] = 3.0*aj;
			beta[2][j] = 3.0*(aj*aj + bj);
			beta[3][j] = aj*(aj*aj + 6.0*bj);
			beta[4][j] = 3.0*bj*(aj*aj + bj);
			beta[5][j] = 3.0 * aj * bj*bj;
			beta[6][j] = Math.pow(bj,3);			
		}
		
		double drho = 0.0;
		double dE = 0.0;
		for (int j = 0; j < 2; j++){
			double sum1 = 0.0;
			for (int i = 0; i < 9; i++){
				double ii = (double)i;
				double temp1 = alpha[i][j]*Math.pow(r[j],(i+1))/(ii+1.0);
				sum1 = sum1 + temp1;
			}
			double sum2 = 0.0;			
			for (int k = 0; k < 7; k++){
				double kk = (double)k;
				double temp2 = beta[k][j]*Math.pow(r[j],(k+2))/((kk+1.0)*(kk+2.0)) + beta[k][j]*Math.pow(r[j],(k+1))*(rho - r[j])/(kk+1);
				sum2 = sum2 + temp2;				
			}
			drho = drho + N[j] * 1.0E-06 * sum1;
			dE = dE + N[j] * 1.0E-06 * sum2 / h[j];
		}
		drho = Crho * drho;
		dE = Ce * 4.0 * cosE * dE/ rho;
		dE = dE / MathUtils.ARCSEC2RAD;
		double[] out = new double[2];
		out[0] = drho;
		out[1] = dE;
		return out;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TropoModel x = new TropoModel(GPS_Utils.lambda);
		double [] els = {1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0};
		VectorN E = new VectorN(els).times(MathUtils.DEG2RAD);
		for (int i = 0; i < E.length; i++) {
			double [] out = x.corrections(938.0, 286.0, 0.73, E.get(i), 1000000.0);
			System.out.println(els[i]+" "+out[0]+" "+out[1]);
		}

	}

}
