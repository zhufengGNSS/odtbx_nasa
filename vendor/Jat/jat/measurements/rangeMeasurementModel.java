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
 */

package jat.measurements;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

//import jat.sim.closedLoopSim;
import jat.matvec.data.*;

import java.util.Random;
import jat.alg.estimators.*;
import jat.cm.Constants;
import jat.sim.*;
import jat.spacetime.EarthRef;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.traj.Trajectory;
import jat.util.FileUtil;

public class rangeMeasurementModel implements MeasurementFileModel,MeasurementModel{
	
	private boolean obsfromfile = false;
	private double pseudorange;
	private int prn;
	private PreciseEphemeris[] ephem;
	
	public static VectorN R;
	public static int numStates;
	HashMap hm;// = closedLoopSim.hm;
	//* *NOTE* added argument to default constructor instead of call to static
	//*        at instantiation
	Random generator;
	
	public rangeMeasurementModel(ObservationMeasurement om, HashMap h){
		hm = h;
		obsfromfile = true;
		pseudorange = om.get_obs_data(ObservationMeasurement.DATA_PSEUDORANGE);
		prn = om.get_PRN();
		
		/*Add a sleep in here to insure that the Random Number
		 * Seeds don't allign with any other random number generator
		 */
		try {
			Thread.sleep(20);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		generator = new Random();
	}
	
	public rangeMeasurementModel(HashMap h) {
		obsfromfile = initializer.parseBool(h,"init.fromfile");
		hm = h;
		
		if(obsfromfile){
			int n = initializer.parseInt(hm,"RANGE.numsc");
			String[] fileName = new String[n];
			ephem = new PreciseEphemeris[n];
			for(int i=0; i<n; i++){
				ephem[i] = new PreciseEphemeris(initializer.parseString(hm,"RANGE.file."+i));
			}
		}
		/*Add a sleep in here to insure that the Random Number
		 * Seeds don't allign with any other random number generator
		 */
		try {
			Thread.sleep(20);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		generator = new Random();
	}
	
	public double getMeasurement()
	{
		String tmp;
//		double[] truth0 = closedLoopSim.truth[0].sc.get_spacecraft().toStateVector();
//		double[] truth1 = closedLoopSim.truth[1].sc.get_spacecraft().toStateVector();
		double[] truth0 = EstimatorSimModel.truth[0].get_spacecraft().toStateVector();
		double[] truth1 = EstimatorSimModel.truth[1].get_spacecraft().toStateVector();

		double x2 = (truth0[0] - truth1[0])*(truth0[0] - truth1[0]);
		double y2 = (truth0[1] - truth1[1])*(truth0[1] - truth1[1]);
		double z2 = (truth0[2] - truth1[2])*(truth0[2] - truth1[2]);
		
		double range = Math.sqrt(x2 + y2 + z2);
		
		//Add in the measurement noise read out of the file
		tmp = "MEAS."+EKF.measNum+".R.0";
		double R = initializer.parseDouble(hm,tmp);
			
		/*Scale the error as Gaussian noise times the 
		 square of the measurement noise*/
		range += generator.nextGaussian()*R*R; 
			
		
		return range;
	}
	
	public double predictMeasurement(VectorN state){
		
		double x2 = (state.get(0) - state.get(6))*(state.get(0) - state.get(6));
		double y2 = (state.get(1) - state.get(7))*(state.get(1) - state.get(7));
		double z2 = (state.get(2) - state.get(8))*(state.get(2) - state.get(8));
		
		double range = Math.sqrt(x2 + y2 + z2);
		
		return range;
		
	}
	/**
	 * For use with crosslink precise ephemeris file data
	 * @return predicted range [meters]
	 */
	public double predictMeasurement(VectorN state, ObservationMeasurement om){
		int id = om.get_PRN();
		VectorN obs_r = get_obs_pos(id,om.get_mjd());
		VectorN rel = state.get(0,3).minus(obs_r);
		double range = rel.mag();
		//* TODO check that the range is being calculated correctly
		return range;
	}
	
	public double zPred(int i, double time, VectorN state){
		double oMinusC;
		double pred = predictMeasurement(state);
		double obs= getMeasurement();
		oMinusC      = obs-pred;
		return oMinusC;
	}
	public double zPred(ObservationMeasurement om, int i, double time, VectorN state){
		double oMinusC;
		double pred = predictMeasurement(state,om);
		double obs = om.get_range();
		oMinusC      = obs-pred;
		return oMinusC;
	}
	
	public double R()
	{	
		String tmp = "MEAS."+EKF.measNum+".R."+0;
		double R = initializer.parseDouble(hm,tmp);
		return R;
	}
	public double R(ObservationMeasurement om)
	{
		return R();
	}
	
	public VectorN H(VectorN state)
	{
		/*NOTE:  Relative state are computed by differencing
		 		 the host (states 6 - 11) minus the local 
		 		 (states 0-5)  */
			
		
		/*Determine the number of states*/
		int numStates = initializer.parseInt(hm,"FILTER.states");
		VectorN H = new VectorN(numStates);
		
		double x2 = (state.get(0) - state.get(6))*(state.get(0) - state.get(6));
		double y2 = (state.get(1) - state.get(7))*(state.get(1) - state.get(7));
		double z2 = (state.get(2) - state.get(8))*(state.get(2) - state.get(8));
		
		double range = Math.sqrt(x2 + y2 + z2);
		
		H.set(0,(state.get(0)-state.get(6))/range);
		H.set(1,(state.get(1)-state.get(7))/range);
		H.set(2,(state.get(2)-state.get(8))/range);
		H.set(3,0.0);
		H.set(4,0.0);
		H.set(5,0.0);
		H.set(6,(state.get(6)-state.get(0))/range);
		H.set(7,(state.get(7)-state.get(1))/range);
		H.set(8,(state.get(8)-state.get(2))/range);
		H.set(9,0.0);
		H.set(10,0.0);
		H.set(11,0.0);
		H.set(12,0.0);
		H.set(13,0.0);
		H.set(14,0.0);
		return H;
	}
	public VectorN H(ObservationMeasurement om, VectorN state){
		/*NOTE:  Relative state are computed by differencing
		 the host (states 6 - 11) minus the local 
		 (states 0-5)  */
		int id = om.get_PRN();
		VectorN obs_r = get_obs_pos(id,om.get_mjd());
		VectorN rel = state.get(0,3).minus(obs_r);
		
		/*Determine the number of states*/
		int numStates = initializer.parseInt(hm,"FILTER.states");
		VectorN H = new VectorN(numStates);
		
		//double x2 = (state.get(0) - state.get(6))*(state.get(0) - state.get(6));
		//double y2 = (state.get(1) - state.get(7))*(state.get(1) - state.get(7));
		//double z2 = (state.get(2) - state.get(8))*(state.get(2) - state.get(8));
		
		double range = rel.mag();//Math.sqrt(x2 + y2 + z2);
		
		H.set(0,(state.get(0)-obs_r.get(0))/range);
		H.set(1,(state.get(1)-obs_r.get(1))/range);
		H.set(2,(state.get(2)-obs_r.get(2))/range);
		H.set(3,0.0);
		H.set(4,0.0);
		H.set(5,0.0);
		//* TODO Should these be set to zero?
		H.set(6,0.0);
		H.set(7,0.0);
		H.set(8,0.0);
		
		//H.set(6,(obs_r.get(0)-state.get(0))/range);
		//H.set(7,(obs_r.get(1)-state.get(1))/range);
		//H.set(8,(obs_r.get(2)-state.get(2))/range);
		//H.set(9,0.0);
		//H.set(10,0.0);
		//H.set(11,0.0);
		//H.set(12,0.0);
		//H.set(13,0.0);
		//H.set(14,0.0);
		return H;
	}
	
	private VectorN get_obs_pos(int id,double mjd){
		for(int i=0; i<ephem.length; i++){
			if(ephem[i].id == id){
				return ephem[i].traj.getPositionAt(mjd);
			}
		}
		//* TODO watch for errors
		return new VectorN(3);
	}
	
	private static class PreciseEphemeris {
		
		public int id;
		public Trajectory traj;
		private EarthRef earth;
		
		public PreciseEphemeris(String fileName){
			traj = new Trajectory();
			try {
				String path = FileUtil.getClassFilePath("jat.measurements","rangeMeasurementModel");
				String fs = FileUtil.file_separator();
				BufferedReader in = new BufferedReader(new FileReader(new File(path+fs+fileName)));
				String line = "start";
				for(int i=0; i<23; i++){
					line = in.readLine();
				}
				RotationMatrix rot;
				earth = new EarthRef(new Time(53000));
				boolean first = true;
				while(!line.equalsIgnoreCase("EOF")){
					StringTokenizer tok = new StringTokenizer(line, " ");
					int year,month,day,hour,min;
					double sec,mjd;
					double[] x = new double[3];
					tok.nextToken();
					year = Integer.parseInt(tok.nextToken());
					month = Integer.parseInt(tok.nextToken());
					day = Integer.parseInt(tok.nextToken());
					hour = Integer.parseInt(tok.nextToken());
					min = Integer.parseInt(tok.nextToken());
					sec = Double.parseDouble(tok.nextToken());
					//* TODO watch this!!!
					sec = 0;
					Time time = new Time(year,month,day,hour,min,sec);
					mjd = time.mjd_utc();
					line = in.readLine();
					tok = new StringTokenizer(line, " ");
					tok.nextToken();
					id = Integer.parseInt(tok.nextToken());
					x[0] = Double.parseDouble(tok.nextToken());
					x[1] = Double.parseDouble(tok.nextToken());
					x[2] = Double.parseDouble(tok.nextToken());
					line = in.readLine();
					tok = new StringTokenizer(line, " ");
					tok.nextToken();
					tok.nextToken();
					VectorN r = new VectorN(x);
					earth.update(time);
					//rot = new RotationMatrix((earth.eci2ecef(time)).transpose());
					//r = rot.transform(r);
					x[0] = Double.parseDouble(tok.nextToken());
					x[1] = Double.parseDouble(tok.nextToken());
					x[2] = Double.parseDouble(tok.nextToken());
					VectorN v = new VectorN(x);
					VectorN data = ecf2eci(r,v,time);
					traj.add(mjd,data.x);
					line = in.readLine();					
				}
				
				
				
			} catch (FileNotFoundException e) {
				System.err.println("Precise Ephemeris file: "+fileName+" not found.");
				e.printStackTrace();
				System.exit(0);
			} catch (IOException ioe){
				System.err.println("Error reading from Precise Ephem File: "+fileName);
				ioe.printStackTrace();
				System.exit(0);
			}
		}
		
		private Matrix computePole(Time t){
			double a1 = -0.097873978689492;
			double a2 = -0.024595601235973;
			double a3 = -0.65797118746339;
			double a4 =  0.095555969865679;
			double a5 =  0.94278239745940;
			double a6 =  0.13268750786818;
			double a7 = -0.29741426364758;
			double a8 = -0.38409384587563;
			double a9 =  0.62823543875476;
			double a10 = 0.53300376467789;
			double Tp = 51013.0;
			double A = 2*Constants.pi/365.25*(t.mjd_utc()-Tp);
			double C = 2*Constants.pi/435*(t.mjd_utc()-Tp);
			double xp = a1 + a2 *Math.cos(A) + a3* Math.sin(A) + a4*Math.cos(C) + a5*Math.sin(C);
			double yp = a6 + a7 *Math.cos(A) + a8* Math.sin(A) + a9*Math.cos(C) + a10*Math.sin(C);
			xp = xp*Constants.arcsec2rad;
			yp = xp*Constants.arcsec2rad;
			Matrix out = new Matrix(3);
			out.A[0][2] = xp;
			out.A[2][0] = -xp;
			out.A[1][2] = -yp;
			out.A[2][1] = yp;
			return out;
		}
		private VectorN ecf2eci(VectorN recf, VectorN vecf, Time t){
//	  	Compute derivative of GHA Matrix (S) and its transpose
	    	double omega = Constants.WE_WGS84;
	    	Matrix todMatrix = earth.trueOfDate(t.mjd_tt());
	        Matrix ghaMatrix = earth.GHAMatrix(t.mjd_ut1(), t.mjd_tt());	   	        
	        Matrix poleMatrix = computePole(t);
	        
	        Matrix A = poleMatrix.times(ghaMatrix);
	        Matrix E = A.times(todMatrix);
	    	VectorN omegaE = new VectorN(0,0,omega);
	    	
//	    	% ---- perform transformations
//	        thetasa= 7.29211514670698e-05 * (1.0  - lod/86400.0 );
//	        omegaearth = [0; 0; thetasa;];

//	        rpef = pm'*recef;
//	        reci = prec'*nut'*st'*rpef;
	    	VectorN rpef = poleMatrix.transpose().times(recf);
	    	VectorN reci = E.transpose().times(recf);
//	        vpef = pm'*vecef;
//	        veci = prec'*nut'*st'*(vpef + cross(omegaearth,rpef));
	    	VectorN vpef = poleMatrix.transpose().times(vecf);
	    	VectorN veci = todMatrix.transpose().times(ghaMatrix.transpose()).times(vpef.plus(omegaE.crossProduct(rpef)));
	    	VectorN out = new VectorN(reci,veci);
	    	return out;
	    }
	}
	
	public static void main(String[] args0){
		String file = "C:/Code/Misc/PreciseEphem/test3a.auro2.prec";
		rangeMeasurementModel.PreciseEphemeris ephem = new rangeMeasurementModel.PreciseEphemeris(file);
		int end = 0;
		end++;
	}
}