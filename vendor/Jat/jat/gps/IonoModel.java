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
 * File Created on Jun 19, 2003
 */

import jat.matvec.data.*;
//import jat.math.*;
import jat.timeRef.*;
//import jat.gps.*;


/**
 * The IonoModel.java Class provides a simple model for range errors
 * due to ionospheric delay.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class IonoModel {

	/** GPS L1 frequency */
	private static final double freq = 1575.42E+06;

	/** Reference TVEC in electrons/m^2 */
	private double tvec_r = 2.0E+17;

	public double correlationTime = 2000.0;

	public double sigma = 0.82;

	public double Ivbar = 5.2;

	private double q;

	/** Default constructor */
	public IonoModel(){
		this.q = 2.0 * this.sigma*this.sigma /this.correlationTime;

	}

	/** Construct an IonoModel with a given reference TVEC
	 * @param ref_tvec reference TVEC in electrons/m^2
	 */
	public IonoModel(double ref_tvec){
		this.tvec_r = ref_tvec;
		this.q = 2.0 * this.sigma*this.sigma /this.correlationTime;
	}

	/** Compute the range error due to iono
	 * @param t_mjd the current GPS time in MJD format
	 * @param r the current spacecraft ECI position vector
	 * @param rGPS the current GPS SV ECI position vector
	 * @return the range error in meters for code measurement
	 * (multiply by -1.0 for carrier phase measurement)
	 */
	public double error(double t_mjd, VectorN r, VectorN rGPS, double iv) {

		// compute signal elevation angle
		double elevation = GPS_Utils.elevation(r, rGPS);
		double out = del_iono(iv, elevation);
		return out;
	}

	/** OD Toolbox interface to compute the range error due to iono
	 * @param t_mjd the current GPS time in MJD format
	 * @param r the current spacecraft ECI position vector
	 * @param rGPS the current GPS SV ECI position vector
	 * @return the range error in meters for code measurement
	 * (multiply by -1.0 for carrier phase measurement)
	 */
	public double error(double t_mjd, double[] r, double[] rGPS, double iv) {
		VectorN r1 = new VectorN(r);
		VectorN r2 = new VectorN(rGPS);

		return error(t_mjd,r1,r2,iv);
	}

	/**
	 * Compute Iv
	 * @param t_mjd time in MJD
	 * @param r VectorN containing ECI position vector
	 */
	public double Iv(double t_mjd, VectorN r){
		// get the sun unit vector
		GPSTimeFormat gpstime = new GPSTimeFormat(t_mjd);
		CalDate utc = gpstime.GPS2UTC();
		EarthRef er = new EarthRef(utc);
		VectorN sun = er.sunVector();
		VectorN usun = sun.unitVector();

		// get the spacecraft unit vector
		VectorN ur = r.unitVector();

		// compute tvec
		double sdotr = usun.dotProduct(ur);
		double tvec = this.tvec_r * Math.pow((1.0 + 0.143*sdotr),8.0);

		double out = 40.3 * tvec / (freq * freq);
		return out;
	}

	/**
	 * OD Toolbox interface to compute Iv
	 * @param t_mjd time in MJD
	 * @param r VectorN containing ECI position vector
	 */
	public double Iv(double t_mjd, double[] r){
		VectorN r1 = new VectorN(r);
		return Iv(t_mjd,r1);
	}
	
	/** compute the range error due to Iono
	 * @param Iv range error along the vertical direction
	 * @param elev line of sight elevation angle in radians
	 */	
	public static double del_iono(double Iv, double elev){
		double sinE = Math.sin(elev);
		double sin2E = sinE * sinE;		
		double out = 2.04 * Iv / (sinE + Math.sqrt(sin2E + 0.076));
		return out;
	}

	public double ionoProcess(double div){
		double out = -1.0 * div/this.correlationTime;
		return out;
	}

	public double ionoQ(double dt){
		double sig2 = this.sigma*this.sigma;
		double exp = Math.exp(-2.0*dt/this.correlationTime);
		double out = sig2 * (1.0 - exp);
		return out;
	}

	public double Q(){
		return this.q;
	}



}
