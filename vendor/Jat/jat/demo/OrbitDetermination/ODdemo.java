package jat.demo.OrbitDetermination;

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

/**
* The ODdemo.java Class is a demonstration of Orbit Determination.
* It processes range data in the file OBSDATA2 with an observation model
* defined by the ObsData class and dynamics model defined by J2DragProcss
* using an Extended Kalman Filter.
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class ODdemo {

	/**
	 * main - runs the demo.
	 * @params args none.
	 */
	public static void main(String[] args) {
		ObsData obs = new ObsData();
		LinePrinter lp = new LinePrinter();
		LinePrinter lp2 = new LinePrinter("C:\\Jat\\output\\resid.txt");
		ProcessModel process = new J2DragProcess(lp, lp2);		
		ExtendedKalmanFilter ekf = new ExtendedKalmanFilter(obs, obs, process);
		System.out.println("Processing..");
		ekf.process();
		System.out.println("Processing completed");

	}
}
