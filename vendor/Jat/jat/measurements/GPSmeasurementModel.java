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

package jat.measurements;

import jat.alg.estimators.EKF;
import jat.alg.estimators.MeasurementFileModel;
import jat.alg.estimators.MeasurementModel;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.math.*;
import jat.sim.*;
import jat.spacetime.EarthRef;
import jat.spacetime.FitIERS;
import jat.spacetime.GPSTimeFormat;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.traj.*;
import jat.util.FileUtil;
import java.util.HashMap;
import jat.gps.*;

import java.io.*;


import ptolemy.plot.*;

import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.JFrame;

/**
* The GPS_MeasurementGenerator.java Class generates GPS measurements
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
* 
* NOTE:  This model assumes that the the satellite will have 
* position and velocity in states 0-5.  The clock bias and
* drift states are defined in the input file
*/
public class GPSmeasurementModel implements MeasurementFileModel,MeasurementModel{


	private boolean obsfromfile = false;
	public static LinePrinter gps_rv;
	
	private ReceiverModel rcvr;
	public GPS_Constellation constell;
	private IonoModel iono = new IonoModel();
	private URE_Model ure;	
	private static GEO_Blockage_Models block;
	private double [] integerAmbiguity;
	static public  LinePrinter noiseOut;
	static public  LinePrinter azOut;
	static public  LinePrinter elOut;
	private Visible vis;
	public static int CurrentSv;
	public double numvis;
	public static HashMap hm;
	static boolean [] AcquisitionFlag = new boolean[33];
	double AcquisitionThreshold = 25;// (dB-Hz)
	double TrackingThreshold    = 21;// (dB-Hz)
	double[] oldclock = this.initClock();
	double[] newclock = oldclock;
	double currentTime = 0;
	double lastTime = 0;
	double MJD0;
	boolean firstTime = true;
	int clockState, biasState;
	public static VectorN H;
	static public boolean visible;
	public static double oldtime, newitme;
	public static VectorN Cn0_out;
	private static FileOutputStream visableSats;
	
	private int FILTER_states;
	
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
	public GPSmeasurementModel(HashMap h) {
		obsfromfile = initializer.parseBool(h,"init.fromfile");
		hm = h;  //closedLoopSim.hm; //* *NOTE* added argument rather than static var
		this.FILTER_states = initializer.parseInt(hm,"FILTER.states");
		block = new GEO_Blockage_Models();
		
		//Set up the GPS constellation
		String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"input"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        String fileName = initializer.parseString(hm,"GPS.const");
        constell = new GPS_Constellation(dir_in+fileName);
		
		// Set up the Receiver Model
		this.rcvr = new ReceiverModel();
		
		int nsv = this.constell.size();
		this.ure = new URE_Model(nsv);
		
		//Determine which state is the clock state
		clockState = initializer.parseInt(hm,"FILTER.clock");
		biasState  = initializer.parseInt(hm,"FILTER.bias");
		
		MJD0 = initializer.parseDouble(hm,"init.MJD0") + initializer.parseDouble(hm,"init.T0")/86400.0;
		Cn0_out = new VectorN(33);
		dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
		String fileName5 = dir_in+"Visible.txt";
		try {
			visableSats = new FileOutputStream(fileName5);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
			
			
		// initialize integer ambiguity if it is used
		this.integerAmbiguity = new double[nsv];
		RandomNumber rn = new RandomNumber();		
		for (int j = 0; j < nsv; j++) {
			double number = rn.normal(0.0, 1.0E+06);
			double num = MathUtils.round(number);
			this.integerAmbiguity[j] = GPS_Utils.lambda*num;
		}
			
		try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        gps_rv = new LinePrinter(dir_in+"gps_rv.txt");
	}

	
	private double[] initClock() {
		double[] out = new double[2];
		out[0] = (Math.random() - 0.5)*10;
		out[1] = .1;
			
		
		firstTime = true;
		return out;
	}

	

	public double R(ObservationMeasurement om)
	{
		int measNum    = 0;//EKF.measNum;
		String tmp = "MEAS."+measNum+".R";
		double sigma = initializer.parseDouble(hm,tmp);
		double R = sigma*sigma;
		return R;
	}
	public double R()
	{
		int measNum    = EKF.measNum;
		String tmp = "MEAS."+measNum+".R";
		double sigma = initializer.parseDouble(hm,tmp);
		double R = sigma*sigma;
		return R;
	}
	//* Cheating *
	//* NOTE * H is actually computed in the call to "predictedMeasurement"
	public VectorN H(ObservationMeasurement om,VectorN state)
	{
		//int numStates = initializer.parseInt(hm,"FILTER.states");
		//VectorN H = new VectorN(numStates);
		return H;
	}
	//* Cheating *
	//* NOTE * H is actually computed in the call to "predictedMeasurement"
	public VectorN H(VectorN state)
	{
		//int numStates = initializer.parseInt(hm,"FILTER.states");
		//VectorN H = new VectorN(numStates);
		return H;
	}

	
	public double zPred(int isv, double t_mjd, VectorN state){
		double out,obs;
		
		
		
		obs = observedMeasurement(isv,t_mjd,state);
		double pred = predictedMeasurement(isv, t_mjd, state);
		
		if(obs == 0)
		    out = 0.0;
		else
			out = obs-pred;	
		return out;
	}
	
	public double zPred(ObservationMeasurement om, int isv, double t_mjd, VectorN state){
		
		double out,obs;
		if(obsfromfile){
			obs = om.get_range();
		}else{
			obs  = observedMeasurement(isv, t_mjd, state);
		}
		//*TODO watch - the following line is a random attempt at hard-tweaking data
		//t_mjd = t_mjd+30*60/86400;
		//double pred = predictedMeasurement(isv, t_mjd, state, om.get_ECF2ECI());
		double pred = predictedMeasurement(isv, t_mjd, state);
	    //if(Double.isNaN(pred)) return Double.NaN;
		if(obs == 0)
		    out = 0.0;
		else
			out = obs-pred;
		
		//*TODO watch should this be absolute value?
		//return Math.abs(out);
		return out;
	}
	public double observedMeasurement(int isv, double t, VectorN state)
	{
		double[] Cn0 = new double[33];
		double out,dt;
		
		
		currentTime = t/(double)86400 + MJD0;
		if(firstTime == true)
		{
			double filterStart = initializer.parseDouble(hm,"init.T0");
			lastTime = currentTime;
			dt = t-filterStart;
//			Propagate the clock forward
			newclock = rcvr.propClock(dt, oldclock);

			oldtime = 99;
			firstTime = false;
		}
		else if(currentTime > lastTime && firstTime == false)
		{
			dt = (currentTime - lastTime)*86400;
			double[] tmpClock = newclock;
			newclock = rcvr.propClock(dt, oldclock);
			oldclock = tmpClock;
			lastTime = currentTime;
		}
		
		//state= new VectorN (closedLoopSim.truth[0].sc.get_spacecraft().toStateVector());
//* TODO Cheating - static reference to EstimatorSimModel		
		//state= new VectorN (CEVSim.truth_traj[0].getState(0));
		state = new VectorN(EstimatorSimModel.truth[0].get_spacecraft().toStateVector());
		
		
		
		// get the SVID and MJD time of measurement
		GPS_SV sv = constell.getSV(isv);
		int prn = sv.prn();
		
		double t_mjd = t/(double)86400 + MJD0; 
		
//		compute time of transmission
		VectorN r = new VectorN(state.get(0,3));
		VectorN v = new VectorN(state.get(3,3));
		double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		
		
//		compute the GPS position vector
		
		/*Follow the method used in the Predicted Measurement method*/
	    Time time = new Time(ts_mjd);
		FitIERS iers = new FitIERS();
		double[] params = iers.search(time.mjd_utc());
		time.set_UT1_UTC(params[2]);
		EarthRef earth = new EarthRef(time);
		earth.setIERS(params[0], params[1]);
		
		//Previous Implimentation
		//Time time = new Time(ts_mjd);
		//EarthRef earth = new EarthRef(time);
		
		Matrix pole = earth.PoleMatrix();
		Matrix gha = earth.GHAMatrix(time.mjd_ut1(),time.mjd_tt());
		Matrix tod = earth.TOD();
		VectorN rvGPS = sv.rvECI(ts_mjd,pole,gha,tod);
		VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
		VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);
		
		// compute the line of sight vector
		VectorN los = GPS_Utils.lineOfSight(r, rGPS);
		VectorN losu = los.unitVector();
		
		
		boolean visibleBlockage = block.visible(losu, r, rGPS);

		visible = acquisitionAndTracking(visibleBlockage, prn, t);
		
		if (visible) {
			
			// compute measured range
			double truerho = los.mag();
			double clock_err = rcvr.clockError(los, v, vGPS, oldclock[0]);
			double code_noise = rcvr.codeNoise();
			
			double cp_noise = rcvr.cpNoise();
			
			//double ionoDelay = iono.error(t_mjd, this.r, rGPS, iv);
			double ionoDelay = 0;
			
			double uree = ure.ure(isv, los, rGPS, vGPS);
			
			//Use 2 meters of random noise for URE
			
			/*Added in the bias correction from the ephemeris here, because it is in the 
			 * predicted measurement method.  Not sure why it was omitted here*/
		    double rho = truerho + clock_err + code_noise + 4*(Math.random()-0.5) - GPS_Utils.c*sv.biasCorrection(new GPSTimeFormat(ts_mjd).gps_sow());//+ code_noise;// + 4*(Math.random()-0.5);
			
		    //double rho = truerho + 4*(Math.random()-0.5) - GPS_Utils.c*sv.biasCorrection(new GPSTimeFormat(ts_mjd).gps_sow());
		    
		    //double rho = truerho + clock_err + uree + code_noise;
			//double rho = truerho + clock_err + ionoDelay + uree + code_noise;
			
			double ia = integerAmbiguity[isv];					
//			double carrier_rho = truerho + clock_err - ionoDelay + uree + ia + cp_noise;
			double carrier_rho = truerho + clock_err + uree + ia + cp_noise;
//			double carrier_rho = truerho + clock_err - ionoDelay + uree + cp_noise;
			
//			double diff = carrier_rho - rho;
//			double twoiono = -2.0*ionoDelay;
			
			out = rho;
		}
		else
			out = 0;
		
		return out;
	}
	
	
	public boolean acquisitionAndTracking(boolean visIn, int prn, double t)
	{
		/*Todo:  Need to be able to get the Cn0 without
		 * directly referencing the visibility file!!
		 * sort of stuck with the way the visible interface
		 * is set up
		 */
		double Cn0 = GEO_Blockage_Models.Cn0;
		
		boolean visible = false;
		if(visIn == false)
			AcquisitionFlag[prn] = false;
		
		//System.out.println("Cn0 for PRN : " + prn + " is: "+Cn0);
		//If it is visible and strong enough, "Acquire" satellite
		if(visIn && AcquisitionFlag[prn] == false && Cn0 > AcquisitionThreshold)
			AcquisitionFlag[prn] = true;
		
		//If we have "acquired" a satellite and if the Cn0 is greater than
		//the tracking threshold, then it is visible
		if(AcquisitionFlag[prn])
		{
			if(Cn0 > TrackingThreshold)
			{
				visible = true;
			}
		}
		else
		{
			visible = false;
		}
		//System.out.println("Cn0 Value: " + Cn0 + " PRN: " + prn + " Visible: " + visible + " Time :" + t);
	
		if(visible)
		{
			if(t != oldtime)
			{	
				//The time has been updated, so end out the information
				new PrintStream(visableSats).println (Cn0_out.toString());
				//Reset the Cn/0 vector
				Cn0_out = new VectorN(33);
				Cn0_out.set(0,t);
				Cn0_out.set(prn,Cn0);
				oldtime = t;
			}	
			else
			{
				Cn0_out.set(prn,Cn0);
			}
		}
		
	return visible;
	}	
	
	public double predictedMeasurement(int isv, double t, VectorN state) 
	{
	//private double rangePred(double t_mjd, VectorN xref, int prn, int type ){
		
		// get the SVID of measurement
		GPS_SV sv;
		//int index = isv;//constell.getIndex(isv);
		//int index = constell.getIndex(isv);
		sv = constell.getSV(isv);
		int prn = sv.prn();
//		if(prn!=13)
//			return Double.NaN;
		int clockIndex = clockState;

		// extract the current spacecraft position vector		
		VectorN r = new VectorN(state.get(0,3));
		VectorN v = new VectorN(state.get(3,3));
			
		//* The following is used when feeding in the sim_time
		//* currently when using files, the mjd_utc time is fed in
		double t_mjd = t/(double)86400 + MJD0;
		//double t_mjd = t; 
		double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		
		//* TODO watch this
		// compute the GPS SV position vector at transmission time
//		* TODO watch this - should this be transmit time?
		//Time time = new Time(t_mjd);
		Time time = new Time(ts_mjd);
		FitIERS iers = new FitIERS();
		double[] params = iers.search(time.mjd_utc());
		time.set_UT1_UTC(params[2]);
		EarthRef earth = new EarthRef(time);
		earth.setIERS(params[0], params[1]);
		
		Matrix pole = earth.PoleMatrix();
		Matrix gha = earth.GHAMatrix(time.mjd_ut1(),time.mjd_tt());
		Matrix tod = earth.TOD();
		VectorN rvGPS_ECEF = sv.rvECEF(ts_mjd);
		VectorN rvGPS = sv.rvECI(ts_mjd,pole,gha,tod);
		VectorN rvGPS_ECEF_meas = sv.rvECEF(t_mjd);
		gps_rv.println("t:\t"+(new Time(ts_mjd).secOfDay())+"\t prn:\t"+sv.prn()+"\t rveci:\t"+rvGPS.toString()+"\t rvecef:\t"+rvGPS_ECEF.toString()+"\t rvecef_meas:\t"+rvGPS_ECEF_meas.toString());

		VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
		VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);
		
		// compute the LOS vector
		VectorN los = rGPS.minus(r);
		VectorN losu = los.unitVector();
		
		// compute the expected range
		double range = los.mag();
				
		// compute the range rate
		double range_rate = GPS_Utils.rangeRate(los, v, vGPS);
		
		// add clock bias contribution to range
		double bc = state.x[clockIndex];
		double omrroc = 1.0 - (range_rate/GPS_Utils.c);
		double bcpart = omrroc*bc;
		range = range + bcpart - GPS_Utils.c*sv.biasCorrection(new GPSTimeFormat(ts_mjd).gps_sow());
		//range = range - GPS_Utils.c*sv.biasCorrection(new GPSTimeFormat(ts_mjd).gps_sow());
		
		// add iono contribution to range
		//int iono_index = 9;
		//double diono = xref.x[iono_index];
		//double iv = (1.0 + diono)*iono.Ivbar;
		//double elev = GPS_Utils.elevation(r, rGPS);
		//double ionopart = IonoModel.del_iono(iv, elev);
		//Assume there is no Ionosphere Delay
		//range = range + ionopart;
		
		// add ure contribution to range
		//int ure_index = 9 + svindex;
		//double dure = xref.x[ure_index];
		//range = range + dure;
		
		//Create and zero out the H vector
		int numStates = this.FILTER_states;//initializer.parseInt(hm,"FILTER.states");
		H = new VectorN(numStates);
		H.set(0.0);
		
		// compute the H vector
		VectorN drrdr = drrdr(los, vGPS, range_rate);
		double term1 = -1.0 * bc / GPS_Utils.c;
		drrdr = drrdr.times(term1);
		VectorN drrdv = this.drrdv(los, vGPS);
		VectorN drhodv = drrdv.times(term1);
		VectorN drhodr = drrdr.minus(losu);
		H.set(0, drhodr);
		H.set((3), drhodv);
		double drdbc = 1.0 - range_rate/GPS_Utils.c;
		
		H.set(clockIndex, drdbc); 
		
		
		
		//double drdiono = IonoModel.del_iono(iono.Ivbar, elev);
		//htilde.set(iono_index, drdiono);
		//htilde.set(ure_index, 1.0);
				
		return range;
	}
	public double predictedMeasurement(int isv, Time t, VectorN state) 
	{
	//private double rangePred(double t_mjd, VectorN xref, int prn, int type ){
		
		// get the SVID of measurement
		GPS_SV sv;
		int index = constell.getIndex(isv);
		sv = constell.getSV(index);
		//int prn = sv.prn();
		
		int clockIndex = clockState;
		
		
		// extract the current spacecraft position vector		
		VectorN r = new VectorN(state.get(0,3));
		VectorN v = new VectorN(state.get(3,3));
			
		double t_mjd = t.mjd_utc(); 
		double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		t.updateTo(ts_mjd);
		
		//* TODO watch this
		// compute the GPS SV position vector at transmission time		
		EarthRef earth = new EarthRef(t);
		
		Matrix pole = earth.PoleMatrix();
		Matrix gha = earth.GHAMatrix(t.mjd_ut1(),t.mjd_tt());
		Matrix tod = earth.TOD();
		VectorN rvGPS = sv.rvECI(ts_mjd,pole,gha,tod);
//		VectorN rvGPSecef = sv.rvECI(ts_mjd);
//		VectorN recef = new VectorN(rvGPSecef.get(0,3));
//		VectorN vecef = new VectorN(rvGPSecef.get(3,3));
//		VectorN rvGPS = earth.ecf2eci(recef,vecef,time);
//		VectorN rvGPS = sv.rvECI(ts_mjd);
		VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
		VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);
		
		// compute the LOS vector
		VectorN los = rGPS.minus(r);
		VectorN losu = los.unitVector();
		
		// compute the expected range
		double range = los.mag();
				
		// compute the range rate
//		VectorN vrel = vGPS.minus(v);
		double range_rate = GPS_Utils.rangeRate(los, v, vGPS);
		
		// add clock bias contribution to range
		double bc = state.x[clockIndex];
		double omrroc = 1.0 - (range_rate/GPS_Utils.c);
		double bcpart = omrroc*bc;
		
		range = range + bcpart;
		
		// add iono contribution to range
		//int iono_index = 9;
		//double diono = xref.x[iono_index];
		//double iv = (1.0 + diono)*iono.Ivbar;
		//double elev = GPS_Utils.elevation(r, rGPS);
		//double ionopart = IonoModel.del_iono(iv, elev);
		//Assume there is no Ionosphere Delay
		//range = range + ionopart;
		
		// add ure contribution to range
		//int ure_index = 9 + svindex;
		//double dure = xref.x[ure_index];
		//range = range + dure;
		
		//Create and zero out the H vector
		int numStates = this.FILTER_states;//initializer.parseInt(hm,"FILTER.states");
		H = new VectorN(numStates);
		H.set(0.0);
		
		// compute the H vector
		VectorN drrdr = drrdr(los, vGPS, range_rate);
		double term1 = -1.0 * bc / GPS_Utils.c;
		drrdr = drrdr.times(term1);
		VectorN drrdv = this.drrdv(los, vGPS);
		VectorN drhodv = drrdv.times(term1);
		VectorN drhodr = drrdr.minus(losu);
		H.set(0, drhodr);
		H.set((3), drhodv);
		double drdbc = 1.0 - range_rate/GPS_Utils.c;
		
		H.set(clockIndex, drdbc); 
		
		
		
		//double drdiono = IonoModel.del_iono(iono.Ivbar, elev);
		//htilde.set(iono_index, drdiono);
		//htilde.set(ure_index, 1.0);
				
		return range;
	}
	public double predictedMeasurement(int isv, double t, VectorN state, RotationMatrix ECF2ECI) 
	{
	//private double rangePred(double t_mjd, VectorN xref, int prn, int type ){
		
		// get the SVID of measurement
		GPS_SV sv;
		int index = constell.getIndex(isv);
		sv = constell.getSV(index);
		//int prn = sv.prn();
		
		int clockIndex = clockState;

		// extract the current spacecraft position vector		
		VectorN r = new VectorN(state.get(0,3));
		VectorN v = new VectorN(state.get(3,3));
			
		double t_mjd = t;//t/(double)86400 + MJD0;
		//* TODO watch this 
		double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		//double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		
		// compute the GPS SV position vector at transmission time
		VectorN rvGPS = sv.rvECI(ts_mjd);
		
		VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
		VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);
//		* TODO watch this for accuracy
		//* sv.rvECI really seems to output as rvWGS84 (EarthFixed)
		//* the velocity should be inertial but still unsure
		rGPS = ECF2ECI.transform(rGPS);
		vGPS = ECF2ECI.transform(vGPS);
		
		// compute the LOS vector
		VectorN los = rGPS.minus(r);
		VectorN losu = los.unitVector();
		
		// compute the expected range
		double range = los.mag();
				
		// compute the range rate
//		VectorN vrel = vGPS.minus(v);
		double range_rate = GPS_Utils.rangeRate(los, v, vGPS);
		
		// add clock bias contribution to range
		double bc = state.x[clockIndex];
		double omrroc = 1.0 - (range_rate/GPS_Utils.c);
		double bcpart = omrroc*bc;
		range = range + bcpart;
		
		// add iono contribution to range
		//int iono_index = 9;
		//double diono = xref.x[iono_index];
		//double iv = (1.0 + diono)*iono.Ivbar;
		//double elev = GPS_Utils.elevation(r, rGPS);
		//double ionopart = IonoModel.del_iono(iv, elev);
		//Assume there is no Ionosphere Delay
		//range = range + ionopart;
		
		// add ure contribution to range
		//int ure_index = 9 + svindex;
		//double dure = xref.x[ure_index];
		//range = range + dure;
		
		//Create and zero out the H vector
		int numStates = this.FILTER_states;//initializer.parseInt(hm,"FILTER.states");
		H = new VectorN(numStates);
		H.set(0.0);
		
		// compute the H vector
		VectorN drrdr = drrdr(los, vGPS, range_rate);
		double term1 = -1.0 * bc / GPS_Utils.c;
		drrdr = drrdr.times(term1);
		VectorN drrdv = this.drrdv(los, vGPS);
		VectorN drhodv = drrdv.times(term1);
		VectorN drhodr = drrdr.minus(losu);
		H.set(0, drhodr);
		H.set((3), drhodv);
		double drdbc = 1.0 - range_rate/GPS_Utils.c;
		
		H.set(clockIndex, drdbc); 
		
		
		
		//double drdiono = IonoModel.del_iono(iono.Ivbar, elev);
		//htilde.set(iono_index, drdiono);
		//htilde.set(ure_index, 1.0);
				
		return range;
	}

	private VectorN drrdr (VectorN los, VectorN vGPS, double rr){
		VectorN vc = vGPS.divide(GPS_Utils.c);
		double range = los.mag();
		double corr = los.dotProduct(vc);
		double factor = 1.0 / (range + corr);
		VectorN losu = los.unitVector();
		VectorN part2 = losu.plus(vc);
		part2 = part2.times(rr);
		VectorN numerator = part2.minus(los);
		VectorN out = numerator.times(factor);
		return out;
	}

	private VectorN drrdv (VectorN los, VectorN vGPS){
		VectorN vc = vGPS.divide(GPS_Utils.c);
		double range = los.mag();
		double corr = los.dotProduct(vc);
		double factor = -1.0 / (range + corr);
		VectorN out = los.times(factor);
		return out;
	}

}
