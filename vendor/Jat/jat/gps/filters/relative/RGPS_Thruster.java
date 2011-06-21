package jat.gps.filters.relative;

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
 * File Created on May 9, 2003
 */
import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.gps.*;
//import jat.gps_ins.*;
//import jat.gps_ins.*;

/**
* The RGPS_Thruster.java Class runs the RGPS-only EKF with thruster model.
* 
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RGPS_Thruster {

	/**
	 * main - runs the demo.
	 * @params args none.
	 */
	public static void main(String[] args) {
		String dir = "C:\\Jat\\jat\\input\\";
		String gpsmeasfile = "gps\\rgpsmeas_rbar_geom4_mp.jat";
		String rinexfile = "gps\\rinex.n";
		String burnfile = "burns\\rbar_burns.jat";
		RGPS_MeasurementList gps = RGPS_MeasurementList.recover(dir+gpsmeasfile);
		GPS_Constellation constellation = new GPS_Constellation(dir+rinexfile);
		LinePrinter lp1 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_rgps_thr_rbar_geom4_mp_new_bt_1.txt");
		LinePrinter lp2 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_rgps_thr_rbar_geom4_mp_new_bt_2.txt");
		LinePrinter lp3 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_rgps_thr_rbar_geom4_mp_resid_bt.txt");
		ProcessModel process = new RGPS_Thruster_ProcessModel(constellation, lp1, lp2, lp3, dir+burnfile);
		MeasurementModel meas = new RGPS_MeasurementModel(gps, constellation);
		ExtendedKalmanFilter ekf = new ExtendedKalmanFilter(meas, gps, process);
		long ts = System.currentTimeMillis();

		ekf.process();
		
		long tf = System.currentTimeMillis();
		double dtf = (tf - ts)/(60.0 * 1000.0);
		System.out.println("EKF Time Elapsed: "+dtf);
		

	}
}
