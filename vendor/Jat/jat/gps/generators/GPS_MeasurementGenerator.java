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
 * File Created on May 22, 2003
 */

package jat.gps.generators;

import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.math.*;
//import jat.cm.*;
import jat.traj.*;
import jat.gps.*;
import java.io.*;

/**
* The GPS_MeasurementGenerator.java Class generates GPS measurements
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class GPS_MeasurementGenerator {

	private Trajectory truth;
	
	private Trajectory iss;

	private static final double t_mjd0 = 51969.0;

	private ReceiverModel rcvr;

	private VectorN r;
	private VectorN v;
	private VectorN rISS;

	private GPS_Constellation constell;
	
	private IonoModel iono = new IonoModel();
	
	private URE_Model ure;

	private RungeKutta8 rk8 = new RungeKutta8(1.0);

	private String outfile;
	
	private double [] integerAmbiguity;
	
	private LinePrinter clockLP;
	
	private LinePrinter measLP;
	
	private Visible vis;

	/**
	 * Constructor
	 * @param t Spacecraft Trajectory
	 * @param i ISS Trajectory
	 * @param c GPS_Constellation
	 * @param v Visible checker
	 * @param file String containing the output file name
	 * @param clp LinePrinter for clock data output
	 * @param mlp LinePrinter for measurement data output
	 */
	public GPS_MeasurementGenerator(
		Trajectory t, Trajectory i,
		GPS_Constellation c, Visible v,
		String file, LinePrinter clp, LinePrinter mlp) {
		
		// generator setup	
		this.truth = t;
		this.iss = i;
		this.constell = c;
		this.outfile = file;
		this.clockLP = clp;
		this.measLP = mlp;
		this.vis = v;
		this.rcvr = new ReceiverModel();
		
		int nsv = this.constell.size();
		this.ure = new URE_Model(nsv);
		
		// initialize integer ambiguity
		this.integerAmbiguity = new double[nsv];
		RandomNumber rn = new RandomNumber();		
		for (int j = 0; j < nsv; j++) {
			double number = rn.normal(0.0, 1.0E+06);
			double num = MathUtils.round(number);
			this.integerAmbiguity[j] = GPS_Utils.lambda*num;
		}
			
	}

	/**
	 * Constructor
	 * @param t Spacecraft Trajectory
	 * @param i ISS Trajectory
	 * @param c GPS_Constellation
	 * @param v Visible checker
	 * @param file String containing the output file name
	 * @param clp LinePrinter for clock data output
	 * @param mlp LinePrinter for measurement data output
	 * @param seed long containing the random number seed to be used
	 */
	public GPS_MeasurementGenerator(
		Trajectory t, Trajectory i,
		GPS_Constellation c, Visible v,
		String file, LinePrinter clp, LinePrinter mlp, long seed) {
		
		// generator setup	
		this.truth = t;
		this.iss = i;
		this.constell = c;
		this.outfile = file;
		this.clockLP = clp;
		this.measLP = mlp;
		this.vis = v;
		this.rcvr = new ReceiverModel(seed);
		
		int nsv = this.constell.size();
		this.ure = new URE_Model(nsv);
		
		// initialize integer ambiguity
		this.integerAmbiguity = new double[nsv];
		RandomNumber rn = new RandomNumber(seed);		
		for (int j = 0; j < nsv; j++) {
			double number = rn.normal(0.0, 1.0E+06);
			double num = MathUtils.round(number);
			this.integerAmbiguity[j] = GPS_Utils.lambda*num;
		}
			
	}
	
	private double[] initClock() {
		double[] out = new double[2];
		out[0] = 0.0;
		out[1] = 0.0;
		return out;
	}
	
	/**
	 * Generate the measurements
	 * @throws IOException
	 */
	public void generate() throws IOException {

		System.out.println("Generating GPS Measurements");
		GPS_MeasurementList list = new GPS_MeasurementList();

		// initialize
		double tprev = 0.0;
		double t = 0.0;
		double[] oldclock = this.initClock();
		double[] newclock = this.initClock();
		double t_mjd = t_mjd0; // MJD at beginning of the sim

		System.out.println("initialized, number of points " + truth.npts());

		// grab the first data
		double[] data = truth.next();
		tprev = data[0];
		this.r = new VectorN(data[1], data[2], data[3]);
		this.v = new VectorN(data[4], data[5], data[6]);
		
		double[] issdata = iss.next();
		this.rISS = new VectorN(issdata[1], issdata[2], issdata[3]);
		
		double iv = iono.Iv(t_mjd, this.r);
		
		int nvis = 0;
		
		// loop
		//		for (int i = 0; i<2; i++){
		while (truth.hasNext()) {
			clockLP.println(t + "\t" + oldclock[0] + "\t" + oldclock[1]+"\t"+nvis+"\t"+iv);
			
			nvis = 0;
			
			// process each GPS satellite
			for (int isv = 0; isv < constell.size(); isv++) {

				// retrieve the SV
				GPS_SV sv = constell.getSV(isv);

				int prn = sv.prn();
				
				// compute time of transmission
				double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);

				// compute the GPS position vector
				VectorN rvGPS = sv.rvECI(ts_mjd);
				VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
				VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);

				// compute the line of sight vector
				VectorN los = GPS_Utils.lineOfSight(this.r, rGPS);
				VectorN losu = los.unitVector();

				// check visibility
				boolean visible = vis.visible(losu, this.r, this.rISS);

				if (visible) {
					
					nvis = nvis + 1;

					// compute measured range
					double truerho = los.mag();
					double clock_err = rcvr.clockError(los, this.v, vGPS, oldclock[0]);
					double code_noise = rcvr.codeNoise();
					
					double cp_noise = rcvr.cpNoise();
					
					double ionoDelay = iono.error(t_mjd, this.r, rGPS, iv);
					
					double uree = ure.ure(isv, los, rGPS, vGPS);
//					double rho = truerho;
					double rho = truerho + clock_err + uree + code_noise;
//					double rho = truerho + clock_err + ionoDelay + uree + code_noise;
					
					double ia = integerAmbiguity[isv];					
//					double carrier_rho = truerho + clock_err - ionoDelay + uree + ia + cp_noise;
					double carrier_rho = truerho + clock_err + uree + ia + cp_noise;
//					double carrier_rho = truerho + clock_err - ionoDelay + uree + cp_noise;
					
//					double diff = carrier_rho - rho;
//					double twoiono = -2.0*ionoDelay;

					// output data
					GPS_Measurement meas =
						new GPS_Measurement(t, t_mjd, rho, carrier_rho, prn);
					list.add(meas);
					double dr = rho - truerho;

					measLP.println(meas + "\t" + isv + "\t" + iv+ "\t" + uree + "\t" + ia);
				}

			}

			// get the next truth data
			data = truth.next();
			t = data[0];
			issdata = iss.next();

			double dt = t - tprev;
			if (dt < 0.0) {
				System.out.println(
					"backwards time jump in GPS_MeasurementGenerator");
				System.exit(99);
			}
			if (dt > 0.0) {
				//				newclock = rk8.integrate(tprev, oldclock, t, this);
				newclock = rcvr.propClock(dt, oldclock);
				t_mjd = t_mjd + dt / 86400.0;
			}
			// update stuff

			this.r = new VectorN(data[1], data[2], data[3]);
			this.v = new VectorN(data[4], data[5], data[6]);
			this.rISS = new VectorN(issdata[1], issdata[2], issdata[3]);
			iv = iono.Iv(t_mjd, this.r);

			oldclock[0] = newclock[0];
			oldclock[1] = newclock[1];
			tprev = t;

		}
		list.serialize(this.outfile);
		clockLP.close();
		measLP.close();
		System.out.println("finished generating GPS measurements");
	}

	public static void main(String[] args) throws IOException {

		// get the true trajectory data
		System.out.println("recovering truth data");
		String in_directory = "C:\\Jat\\jat\\traj\\";
		String trajfile = "rvtraj.jat";
		String file = in_directory + trajfile;
		Trajectory truetraj = Trajectory.recover(file);
				
		// get the ISS data
		String issfile = "isstraj.jat";
		file = in_directory + issfile;
		Trajectory isstraj = Trajectory.recover(file);

		// get the GPS Constellation
		String rinexfile = "C:\\Jat\\jat\\input\\rinex.n";
		GPS_Constellation constellation = new GPS_Constellation(rinexfile);

		// set the output file
		String out_directory = "C:\\Jat\\jat\\output\\";
		String measfile = "C:\\Jat\\jat\\input\\gpsmeas_if_rv.jat";
		String outfile = measfile;
		
		// set the clock output
		String clockfile = "gpsclock_if_rv.txt";
		LinePrinter clp = new LinePrinter(out_directory+clockfile);
		
		// set the text output file
		String msfile = "gpsmeas_if_rv.txt";
		LinePrinter mlp = new LinePrinter(out_directory+msfile);
		
		// visibility checker
		ISS_Blockage block = new ISS_Blockage();
//		ElevationMask block = new ElevationMask();
		
		GPS_MeasurementGenerator x =
			new GPS_MeasurementGenerator(truetraj, isstraj, constellation, block, outfile, clp, mlp);

		x.generate();

	}
}
