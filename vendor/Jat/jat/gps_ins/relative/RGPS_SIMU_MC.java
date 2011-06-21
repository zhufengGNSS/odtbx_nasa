package jat.gps_ins.relative;

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
import jat.ins.*;

/**
* The ODdemo.java Class is a demonstration of Orbit Determination.
* It processes range data in the file OBSDATA2 with an observation model
* defined by the ObsData class and dynamics model defined by J2DragProcss
* using an Extended Kalman Filter.
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RGPS_SIMU_MC {

	/**
	 * main - runs the demo.
	 * @params args none.
	 */
	public static void main(String[] args) {
		
		long counter = 11;
		int j = 11;
		for (int i = 0; i < 20; i++) {
			
			String num = Integer.toString(j);
			long seed = -1 * counter;
			
			
			String dotJat = ".jat";
			String dotTxt = ".txt";
		
			String dir = "C:\\Jat\\jat\\input\\";
			String gpsmeasfile = "monte\\rgpsmeas_rbar_geom3_mp_" + num + dotJat;
			String insmeasfile = "monte\\simu_rbar_" + num + dotJat;
			String rinexfile = "gps\\rinex.n";
			RGPS_MeasurementList gps = RGPS_MeasurementList.recover(dir+gpsmeasfile);
			INS_MeasurementList ins = INS_MeasurementList.recover(dir+insmeasfile);
			GPS_Constellation constellation = new GPS_Constellation(dir+rinexfile);
			LinePrinter lp1 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_simu_rbar_geom1_mp_mc_1_" + num + dotTxt);
			LinePrinter lp2 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_simu_rbar_geom1_mp_mc_2_" + num + dotTxt);
			LinePrinter lp3 = new LinePrinter("C:\\Jat\\jat\\output\\ekf_simu_rbar_geom1_mp_mc_resid_" + num + dotTxt);
			ProcessModel process = new RGPS_SIMU_ProcessModel(ins, constellation, lp1, lp2, lp3, seed);
			MeasurementModel meas = new RGPS_INS_MeasurementModel(gps, constellation);
			ExtendedKalmanFilter ekf = new ExtendedKalmanFilter(meas, gps, process);
			long ts = System.currentTimeMillis();
			
			System.out.println("Run Number: "+num); 
	
			ekf.process();
			
			long tf = System.currentTimeMillis();
			double dtf = (tf - ts)/(60.0 * 1000.0);
			System.out.println("EKF Time Elapsed: "+dtf);

			counter = counter + 1;
			j = j + 1;
			
		}
		

	}
}
