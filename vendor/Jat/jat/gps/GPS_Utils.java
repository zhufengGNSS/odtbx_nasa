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
 * 
 * File Created on Jun 20, 2003
 */
 
import jat.matvec.data.*;

/**
* The GPSutils.java Class provides general utility functions needed for GPS
* data processing.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/

public class GPS_Utils {
	
    /** Speed of Light in m/s from IAU 1976.
     */
    public static final double c = 299792458.0;

	/** GPS L1 frequency */
	public static final double L1_freq = 1575.42E+06;
	
	/** GPS L1 wavelength in meters */
	public static final double lambda = c / L1_freq;
	
	
	
	/** compute the line of sight vector from a spacecraft to a GPS satellite
	 * @param r the spacecraft position vector
	 * @param rGPS the GPS SV position vector
	 * @return the line of sight vector
	 */	
	public static VectorN lineOfSight(VectorN r, VectorN rGPS) {
		VectorN los = rGPS.minus(r);
		return los;
	}
	
	/** Iterative solution for time of transmission
	 * @param t_mjd time of reception in MJD
	 * @param sv GPS_SV object
	 * @return time of transmission in MJD
	 * 
	 */	
	public static double transmitTime(double t_mjd, GPS_SV sv, VectorN r){
		int maxit = 500;
		int i = 0;
		double ts_new = 0.0;
		double ts_old = t_mjd;
		double eps = 1.0E-16;
		double diff = Math.abs(ts_new - ts_old);
		while ((diff > eps)&&(i < maxit)) {
			VectorN rGPS = sv.rECI(ts_old);
			VectorN los = GPS_Utils.lineOfSight(r, rGPS);
			double rho = los.mag();
			ts_new = t_mjd - rho/(86400.0*c);
			diff = Math.abs(ts_new - ts_old);
			ts_old = ts_new;
			i = i + 1;
			if (i >= maxit){
				System.out.println("transmitTime too many iterations, diff = "+ diff);
				if (diff > 1.0E-10) {
					throw new RuntimeException("GPS_Utils.transmitTime too many iterations at t_mjd = "+t_mjd+", diff = "+ diff);
				}
			}
		}
		return ts_new;		
	}	
	
	/** compute the range rate
	 * @param rho the range vector
	 * @param v the spacecraft ECI velocity vector at reception time
	 * @param vGPS the GPS SV ECI velocity vector at transmission time
	 */	
	public static double rangeRate(VectorN rho, VectorN v, VectorN vGPS) {
		
		// compute the relative velocity
		VectorN vrel = vGPS.minus(v);
		double num = rho.dotProduct(vrel);
		
		double timeCorrection = rho.dotProduct(vGPS)/c;
		double denom = rho.mag() + timeCorrection;
		
//		double range_rate = num / rho.mag();  //DMS this was original
		double range_rate = num / denom;
			
		return range_rate;
	}
	
	/** OD Toolbox interface to compute the range rate
	 * @param rho the range vector
	 * @param v the spacecraft ECI velocity vector at reception time
	 * @param vGPS the transmitting SV ECI velocity vector at transmission time
	 */	
	public static double rangeRate(double[] range, double[] vR, double[] vT) {
		
		VectorN rho  = new VectorN(range);
		VectorN v    = new VectorN(vR);
		VectorN vGPS = new VectorN(vT);
		
		return rangeRate(rho,v,vGPS);
	}
	
	/** Compute the line of sight declination angle (angle from zenith).
	 * @param r spacecraft position vector
	 * @param rGPS GPS SV position vector
	 * @return the angle between the GPS line of sight and zenith in radians.
	 */ 	
	public static double declination(VectorN r, VectorN rGPS) {
		double dot = r.dotProduct(rGPS);
		double rmag = r.mag();
		double rGPSmag = rGPS.mag();
		double dec = Math.acos(dot/(rmag*rGPSmag));
		return dec;
	}

	/** Compute the line of sight elevation angle (angle from local horizontal).
	 * @param r spacecraft position vector
	 * @param rGPS GPS SV position vector
	 * @return the angle between the GPS line of sight and horizon in radians.
	 */ 		
	public static double elevation(VectorN r, VectorN rGPS) {
		double dec = GPS_Utils.declination(r, rGPS);
		double pi2 = Math.PI/2.0;
		double el = pi2 - dec;
		return el;
	}
	/** ODToolbox interface to elevation 
	 * Compute the line of sight elevation angle (angle from local horizontal).
	 * @param r spacecraft position vector
	 * @param rGPS GPS SV position vector
	 * @return the angle between the GPS line of sight and horizon in radians.
	 */ 		
	public static double elevation(double[] r, double[] rGPS) {
		VectorN r1 = new VectorN(r);
		VectorN r2 = new VectorN(rGPS);
		double  el = elevation(r1,r2);
	return el;
	}
	
}
