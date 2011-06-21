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
//import jat.gps_ins.*;

/**
* The RGPS_MP_EarthBlock_MeasurementGenerator.java Class generates relative GPS measurements
* including the effects of multipath and blockage due to a spherical earth.

* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RGPS_MP_EarthBlock_MeasurementGenerator {

	private Trajectory truth;

	private Trajectory iss;

	private static final double t_mjd0 = 51969.0;

	private ReceiverModel rcvr1;
	private ReceiverModel rcvr2;

	private VectorN r;
	private VectorN v;
	private VectorN rISS;
	private VectorN vISS;

	private GPS_Constellation constell;
	
	private double arcs = 500.0;
	private int nrays = 5;
	
	private MultipathModel mp = new MultipathModel(arcs, nrays);

	private IonoModel iono = new IonoModel();

	private URE_Model ure;

	private RungeKutta8 rk8 = new RungeKutta8(1.0);

	private String outfile;

	private double[] iaChaser;

	private double[] iaISS;

	private LinePrinter clockLP;

	private LinePrinter measLP;

	private Visible vis1;

	private Visible vis2;

	/**
	 * Constructor
	 * @param t Spacecraft Trajectory
	 * @param i ISS Trajectory
	 * @param c GPS_Constellation
	 * @param v Visible checker
	 * @param file String containing the output file name
	 * @param clp LinePrinter for clock data output
	 * @param mlp LinePrinter for measurement data output
	 * @param seed random number generator seed to be used
	 */
	public RGPS_MP_EarthBlock_MeasurementGenerator(
		Trajectory t,
		Trajectory i,
		GPS_Constellation c,
		Visible v1,
		Visible v2,
		String file,
		LinePrinter clp,
		LinePrinter mlp,
		long seed) {

		// generator setup	
		this.truth = t;
		this.iss = i;
		this.constell = c;
		this.outfile = file;
		this.clockLP = clp;
		this.measLP = mlp;
		this.vis1 = v1;
		this.vis2 = v2;
		this.rcvr1 = new ReceiverModel(seed);
		this.rcvr2 = new ReceiverModel(2 * seed);

		int nsv = this.constell.size();
		this.ure = new URE_Model(nsv);

		// initialize integer ambiguity
		this.iaChaser = new double[nsv];
		this.iaISS = new double[nsv];
		RandomNumber rn1 = new RandomNumber(seed);
		RandomNumber rn2 = new RandomNumber(2 * seed);
		for (int j = 0; j < nsv; j++) {
			double number1 = rn1.normal(0.0, 1.0E+06);
			double num1 = MathUtils.round(number1);
			this.iaChaser[j] = GPS_Utils.lambda * num1;
			double number2 = rn2.normal(0.0, 1.0E+06);
			double num2 = MathUtils.round(number2);
			this.iaISS[j] = GPS_Utils.lambda * num2;
		}

	}

	private double[] initClock1() {
		double[] out = new double[2];
		out[0] = 0.0;
		out[1] = 0.0;
		return out;
	}

	private double[] initClock2() {
		double[] out = new double[2];
		out[0] = 0.0;
		out[1] = 0.0;
		return out;
	}

	private double[] compute(double t, double t_mjd, ReceiverModel rcv, VectorN los, 
	VectorN r, VectorN v, VectorN rGPS, VectorN vGPS, double clock0, double iv, double ia, int isv, int prn, int type) {
			
		double[] out = new double[2];
		double truerho = los.mag();
		double clock_err = rcv.clockError(los, v, vGPS, clock0);
		double code_noise = rcv.codeNoise();
		double ionoDelay = iono.error(t_mjd, r, rGPS, iv);
		double uree = ure.ure(isv, los, rGPS, vGPS);
		
		double rho = truerho + clock_err + uree + ionoDelay + code_noise;

		double cp_noise = rcv.cpNoise();
		double carrier_rho = truerho + clock_err + uree + ia - ionoDelay + cp_noise;
		measLP.println(t+"\t"+t_mjd+"\t"+rho+"\t"+carrier_rho+"\t"+type+"\t"+prn +"\t"+ isv + "\t"+ iv + "\t" + uree+ "\t"+ ia);

		out[0] = rho;
		out[1] = carrier_rho;
		return out;
	}


	private double[] compute_mp(double t, double t_mjd, ReceiverModel rcv, VectorN los, 
	VectorN r, VectorN v, VectorN rGPS, VectorN vGPS, VectorN rISS, double clock0, double iv, double ia, int isv, int prn, int type) {
			
		double[] out = new double[2];
		double truerho = los.mag();
		double clock_err = rcv.clockError(los, v, vGPS, clock0);
		double code_noise = rcv.codeNoise();
		double ionoDelay = iono.error(t_mjd, r, rGPS, iv);
		double uree = ure.ure(isv, los, rGPS, vGPS);
		
		// multipath
		double theta = GPS_Utils.declination(r, rGPS);
		VectorN dr_vec = rISS.minus(r);
		double dr = dr_vec.mag();		
		mp.environment(prn, dr, theta, rGPS, rISS, r);
		double mperr = mp.pseudorangeError();
		double cperr = mp.carrierPhaseError();
		
		double rho = truerho + clock_err + uree + ionoDelay + code_noise + mperr;

		double cp_noise = rcv.cpNoise();
		double carrier_rho = truerho + clock_err + uree + ia - ionoDelay + cp_noise + cperr;
		measLP.println(t+"\t"+t_mjd+"\t"+rho+"\t"+carrier_rho+"\t"+type+"\t"+prn +"\t"+ isv + "\t"+ iv + "\t" + uree+ "\t"+ ia + "\t"+ mperr+ "\t"+ cperr);

		out[0] = rho;
		out[1] = carrier_rho;
		return out;
	}

	/**
	 * Method generate.
	 * @throws IOException
	 */
	public void generate() throws IOException {

		System.out.println("Generating RGPS Measurements");
		RGPS_MeasurementList list = new RGPS_MeasurementList();

		// initialize
		double tprev = 0.0;
		double t = 0.0;
		double[] oldclock1 = this.initClock1();
		double[] newclock1 = this.initClock1();
		double[] oldclock2 = this.initClock2();
		double[] newclock2 = this.initClock2();

		double t_mjd = t_mjd0; // MJD at beginning of the sim

		System.out.println("initialized, number of points " + truth.npts());

		// grab the first data
		double[] data = truth.next();
		tprev = data[0];
		this.r = new VectorN(data[1], data[2], data[3]);
		this.v = new VectorN(data[4], data[5], data[6]);

		double[] issdata = iss.next();
		this.rISS = new VectorN(issdata[1], issdata[2], issdata[3]);
		this.vISS = new VectorN(issdata[4], issdata[5], issdata[6]);

		double iv1 = iono.Iv(t_mjd, this.r);
		double iv2 = iono.Iv(t_mjd, this.rISS);

		int nvis1 = 0;
		int nvis2 = 0;
		int ncommon = 0;

		// loop
		//		for (int i = 0; i<2; i++){
		while (truth.hasNext()) {
			clockLP.println(t + "\t" + oldclock1[0] + "\t"+ oldclock1[1]+"\t"+ oldclock2[0]+ "\t"+ oldclock2[1] +"\t" + nvis1 + "\t"+ nvis2+ "\t" + ncommon+ "\t"+ iv1+ "\t" + iv2);
//			System.out.println(t);
			nvis1 = 0;
			nvis2 = 0;
			ncommon = 0;

			// process each GPS satellite
			for (int isv = 0; isv < constell.size(); isv++) {

				// retrieve the SV
				GPS_SV sv = constell.getSV(isv);

				int prn = sv.prn();

				// compute time of transmission
				double ts_mjd1 = GPS_Utils.transmitTime(t_mjd, sv, r);
				double ts_mjd2 = GPS_Utils.transmitTime(t_mjd, sv, rISS);

				// compute the GPS position vector
				VectorN rvGPS1 = sv.rvECI(ts_mjd1);
				VectorN rGPS1 = new VectorN(rvGPS1.x[0], rvGPS1.x[1], rvGPS1.x[2]);
				VectorN vGPS1 = new VectorN(rvGPS1.x[3], rvGPS1.x[4], rvGPS1.x[5]);
				
				VectorN rvGPS2 = sv.rvECI(ts_mjd2);
				VectorN rGPS2 = new VectorN(rvGPS2.x[0], rvGPS2.x[1], rvGPS2.x[2]);
				VectorN vGPS2 = new VectorN(rvGPS2.x[3], rvGPS2.x[4], rvGPS2.x[5]);

				// compute the line of sight vector
				VectorN los1 = GPS_Utils.lineOfSight(this.r, rGPS1);
				VectorN losu1 = los1.unitVector();
				VectorN los2 = GPS_Utils.lineOfSight(this.rISS, rGPS2);
				VectorN losu2 = los2.unitVector();

				// check visibility
				boolean visible1 = vis1.visible(losu1, this.r, this.rISS);
				boolean visible2 = vis2.visible(losu2, this.rISS, this.rISS);

				double [] meas1 = new double[2];
				double [] meas2 = new double[2];
				
				// if visible by chaser
				if (visible1) {

					nvis1 = nvis1 + 1;
					
					meas1 = this.compute_mp(t, t_mjd, this.rcvr1, los1, this.r, this.v, rGPS1, vGPS1, rISS, oldclock1[0], iv1, this.iaChaser[isv], isv, prn, 0);

					// output data
					RGPS_Measurement meas =
						new RGPS_Measurement(t, t_mjd, meas1[0], 0, prn);
					list.add(meas);

				}
				
				// if visible by ISS
				if (visible2) {

					nvis2 = nvis2 + 1;
					
					meas2 = this.compute(t, t_mjd, this.rcvr2, los2, this.rISS, this.vISS, rGPS2, vGPS2, oldclock2[0], iv2, this.iaISS[isv], isv, prn, 1);

					// output data
					RGPS_Measurement meas =
						new RGPS_Measurement(t, t_mjd, meas2[0], 1, prn);
					list.add(meas);

				}
				
				// if visible by both, then create a single diff cp measurement
				if (visible1 && visible2) {
					ncommon = ncommon + 1;
					double diff = meas1[1] - meas2[1];
					RGPS_Measurement meas =
						new RGPS_Measurement(t, t_mjd, diff, 2, prn);
					list.add(meas);					
				}
								
			}

			// get the next truth data
			data = truth.next();
			t = data[0];
			issdata = iss.next();

			// move time
			double dt = t - tprev;
			if (dt < 0.0) {
				System.out.println(
					"backwards time jump in GPS_MeasurementGenerator");
				System.exit(99);
			}
			if (dt > 0.0) {
				newclock1 = rcvr1.propClock(dt, oldclock1);
				newclock2 = rcvr2.propClock(dt, oldclock2);
				
				t_mjd = t_mjd + dt / 86400.0;
			}

			// update stuff
			this.r = new VectorN(data[1], data[2], data[3]);
			this.v = new VectorN(data[4], data[5], data[6]);
			this.rISS = new VectorN(issdata[1], issdata[2], issdata[3]);
			this.vISS = new VectorN(issdata[4], issdata[5], issdata[6]);
			iv1 = iono.Iv(t_mjd, this.r);
			iv2 = iono.Iv(t_mjd, this.rISS);

			oldclock1[0] = newclock1[0];
			oldclock1[1] = newclock1[1];
			oldclock2[0] = newclock2[0];
			oldclock2[1] = newclock2[1];
			
			tprev = t;

		}
		list.serialize(this.outfile);
		clockLP.close();
		measLP.close();
		System.out.println("finished generating RGPS measurements");
	}

	public static void main(String[] args) throws IOException {

		// get the true trajectory data
		System.out.println("recovering truth data");
		String in_directory = "C:\\Jat\\jat\\traj\\reference\\";
		String trajfile = "rvtraj_vbar.jat";
		String file = in_directory + trajfile;
		Trajectory truetraj = Trajectory.recover(file);

		// get the ISS data
		String issfile = "isstraj_vbar.jat";
		file = in_directory + issfile;
		Trajectory isstraj = Trajectory.recover(file);

		// get the GPS Constellation
		String rinexfile = "C:\\Jat\\jat\\input\\gps\\rinex.n";
		GPS_Constellation constellation = new GPS_Constellation(rinexfile);

		// set the output file
		String out_directory = "C:\\Jat\\jat\\output\\";
		String measfile = "C:\\Jat\\jat\\input\\gps\\rgpsmeas_vbar_geom1_mp_eb2.jat";
		String outfile = measfile;

		// set the clock output
		String clockfile = "gpsref\\rgpsclock_vbar_geom1_mp_eb2.txt";
		LinePrinter clp = new LinePrinter(out_directory + clockfile);

		// set the text output file
//		String msfile = "rgpsmeas_vbar.txt";
		String msfile = "gpsref\\rgpsmeas_vbar_geom1_mp_eb2.txt";
		LinePrinter mlp = new LinePrinter(out_directory + msfile);

		// visibility checker
		ISS_Earth_Blockage block = new ISS_Earth_Blockage();
//		ElevationMask mask = new ElevationMask(10.0);
		Earth_Blockage mask = new Earth_Blockage();
		
		long seed = -1;

		RGPS_MP_EarthBlock_MeasurementGenerator x = new RGPS_MP_EarthBlock_MeasurementGenerator(truetraj,isstraj, constellation, block, mask, outfile, clp, mlp, seed);
		x.generate();

	}
}
