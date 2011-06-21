package jat.alg.estimators;

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
 * File Created on May 7, 2003
 */
import java.util.HashMap;
import jat.sim.*;
import jat.matvec.data.*;
import jat.alg.integrators.LinePrinter;
import jat.eph.DE405;
import jat.measurements.*;
import jat.spacecraft.SpacecraftModel;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.util.FileUtil;

//import jat.audio.*;

/**
 * The ExtendedKalmanFilter Class processes measurements using an EKF algorithm,
 * given the measurements, measurement model and dynamics (or process) model.
* Assumes a scalar measurement update.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class EKF {

	DE405 ephem = new DE405();
	double mjd_epoch;

	/** Dynamics or Process Model */
	public ProcessModel  process;

	private createMeasurements measurements;
	
	//** Next observation */
	//private ObservationMeasurement next_obs;
	
	/** Number of states or unknowns */
	public int n;

	/** Nominal time step */
	private double dtNominal;
	
	/** Time in the filter*/
	 private double filterTime;
	 
	 private double finalTime;
	  
	 /**Filter State**/
	 private double[] xprev;
	 
	 /**Previous Time**/
	 private double tprev;
	 
	 /**State Transition Matrix*/
	 public EstSTM xref;
	 
	 /**Need to keep track of the new and old covariances*/
//	 public static Matrix pold;
//	 public static Matrix pnew;
	 public Matrix pold;
	 public Matrix pnew;
	 
	 public LinePrinter residuals; 
	 public HashMap hm;

	private boolean verbose=false;
	public static boolean visible;
	public static int measNum, stateNum;
	
	private boolean runMonteCarlo = false;
	 
	/**
	 * Constructor.
	 * @param h HashMap from parsing input
	 */
	public EKF(HashMap h) {
		hm = h;
		measurements = new createMeasurements(hm);
		
        String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        try{
        	runMonteCarlo = initializer.parseBool(hm, "init.runMonteCarlo");
        }catch(NullPointerException ne){ runMonteCarlo = false;}
        if(!runMonteCarlo){
        	residuals = new LinePrinter(dir_in+"Residuals.txt");
        }
		this.n = initializer.parseInt(hm,"FILTER.states");
		String stringPm = initializer.parseString(hm,"FILTER.pm");
		dtNominal = initializer.parseInt(hm,"FILTER.dt");
		mjd_epoch = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0")/86400.0;
		
		filterTime = 0;//initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0");
		System.out.println(stringPm);
		if(stringPm.equals("JGM4x4SRPProcess15state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess15state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4SRPProcess9state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess9state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("Simple") || stringPm.equals("Lunar"))
		{
			this.process = new SimpleProcessModel(hm);
		}else{
			System.out.println("Process model not recognized.  Aborting");
			System.exit(1);
		}
		//double[] X = new double[n];
		filterInitialize();
		
	}
	
	/**
	 * Constructor.
	 * @param h HashMap from parsing input
	 */
	public EKF(HashMap h, int JAT_Case) {
		hm = h;
		measurements = new createMeasurements(hm);
		
        String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        try{
        runMonteCarlo = initializer.parseBool(hm, "init.runMonteCarlo");        
        }catch(Exception e){
        	runMonteCarlo = false;
        }
        if(!runMonteCarlo){
        	residuals = new LinePrinter(dir_in+"Residuals.txt");
        }
		this.n = initializer.parseInt(hm,"FILTER.states");
		String stringPm = initializer.parseString(hm,"FILTER.pm");
		dtNominal = initializer.parseInt(hm,"FILTER.dt");
		mjd_epoch = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0")/86400.0;
		
		filterTime = 0;//initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0");
		System.out.println(stringPm);
		if(stringPm.equals("JGM4x4SRPProcess15state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess15state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4SRPProcess9state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess9state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4DragProcess9state")){
			this.process= new JGM4x4DragProcess9state(hm);
		}
		else if(stringPm.equals("Simple") || stringPm.equals("Lunar"))
		{
			this.process = new SimpleProcessModel(hm);
		}else{
			System.out.println("Process model not recognized.  Aborting");
			System.exit(1);
		}
		//double[] X = new double[n];
		filterInitialize();
		
	}
	
	/**
	 * Constructor.
	 * @param ol ObservationMeasurementList initialized from a file
	 */
	public EKF(ObservationMeasurementList ol, HashMap input) {
		//obs_list = ol;
		//next_obs = obs_list.getFirst();
		hm = input;
        String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        try{
        runMonteCarlo = initializer.parseBool(hm, "init.runMonteCarlo");
        }catch(Exception e){
        	runMonteCarlo = false;
        }
        if(!runMonteCarlo){
        	residuals = new LinePrinter(dir_in+"Residuals.txt");
        }
		this.n = initializer.parseInt(hm,"FILTER.states");
        //this.n=6;
		String stringPm = initializer.parseString(hm,"FILTER.pm");
		dtNominal = initializer.parseInt(hm,"FILTER.dt");
		mjd_epoch = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0")/86400.0;
		
		filterTime = 0;//initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0");
		System.out.println(stringPm);
		if(stringPm.equals("JGM4x4SRPProcess15state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess15state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4SRPProcess9state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess9state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4DragProcess9state")){
			this.process= new JGM4x4DragProcess9state(hm);
		}
		else
		{
			System.out.println("Process model not recognized.  Aborting");
			System.exit(1);
		}
		//double[] X = new double[n];
		filterInitialize();
		
	}
	
	public EKF(HashMap input, int jat_case, createMeasurements created_meas) {
		hm = input;
		measurements = created_meas;
		
        String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
        }catch(Exception e){
            dir_in = "";
        }
        runMonteCarlo = initializer.parseBool(hm, "init.runMonteCarlo");
        if(!runMonteCarlo){
        	residuals = new LinePrinter(dir_in+"Residuals.txt");
        }
		this.n = initializer.parseInt(hm,"FILTER.states");
		String stringPm = initializer.parseString(hm,"FILTER.pm");
		dtNominal = initializer.parseInt(hm,"FILTER.dt");
		mjd_epoch = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0")/86400.0;
		
		filterTime = 0;//initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0");
		System.out.println(stringPm);
		if(stringPm.equals("JGM4x4SRPProcess15state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess15state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("JGM4x4SRPProcess9state"))
		{
			LinePrinter lp1 = new LinePrinter(dir_in+"geom1_1.txt");
	 		LinePrinter lp2 = new LinePrinter(dir_in+"geom1_2.txt");
			this.process= new JGM4x4SRPProcess9state(lp1, lp2,hm);
	
		}
		else if(stringPm.equals("Simple") || stringPm.equals("Lunar"))
		{
			this.process = new SimpleProcessModel(hm);
		}else{
			System.out.println("Process model not recognized.  Aborting");
			System.exit(1);
		}
		//double[] X = new double[n];
		filterInitialize();
	}

	public int get_numStates(){
		return this.xref.numberOfStates();
	}
	private Matrix updateCov(VectorN k, VectorN h, Matrix p) {
		Matrix eye = new Matrix(this.n);
		Matrix kh = k.outerProduct(h);
		Matrix i_kh = eye.minus(kh);
		Matrix out = i_kh.times(p);
		return out;
	}

	private Matrix updateCov(Matrix k, Matrix h, Matrix p){
		Matrix eye = new Matrix(this.n);
		Matrix kh = k.times(h);
		Matrix i_kh = eye.minus(kh);
		Matrix out = i_kh.times(p);
		return out;
	}
	
	private Matrix updateCovJoseph(VectorN k, VectorN h, Matrix p, double measurementNoise){
		Matrix eye = new Matrix(this.n);
		Matrix kh = k.outerProduct(h);
		Matrix i_kh = eye.minus(kh);
		Matrix i_khT = i_kh.transpose();
		double r = measurementNoise;
		Matrix kkT = k.outerProduct(k);
		Matrix krkT = kkT.times(r);

		Matrix part1 = i_kh.times(p);
		part1 = part1.times(i_khT);
		Matrix out = part1.plus(krkT);
		return out;
	}

	private VectorN updateState(VectorN k, double z) {
		VectorN xhat = k.times(z);
		return xhat;
	}

	private Matrix propCov(Matrix p, Matrix phi, Matrix q) {
		Matrix phitrans = phi.transpose();
		Matrix temp1 = p.times(phitrans);
		Matrix temp2 = phi.times(temp1);
		Matrix out = temp2.plus(q);
		return out;
	}

	private VectorN kalmanGain(Matrix p, VectorN h, double r) {
		VectorN ph = p.times(h);
		VectorN hp = h.times(p);
		double hph = hp.dotProduct(h);
		double hphr_inv = 1.0 / (hph + r);
		VectorN out = ph.times(hphr_inv);
		
		return out;
	}
	
	private Matrix kalmanGain(Matrix p, Matrix h, Matrix r) {
		Matrix ph = p.times(h);
		Matrix hp = h.times(p);
		Matrix hph = hp.times(h);
		Matrix hphTranspose = hph.transpose();
		//Matrix hphr = (hphTranspose.plus(r));
		Matrix hphr_inv = hphTranspose.inverse();
		Matrix out = ph.times(hphr_inv);
		return out;
	}
	
	
	public void filterInitialize() {
		
		System.out.println("Initializing the Filter .  .  . ");

		// initialize
		tprev = 0.0;
		pold = process.P0();
		pnew = pold.copy();
		xref = new EstSTM(process.xref0());
		xprev = xref.longarray();
		//VectorN k = new VectorN(process.numberOfStates());
				
	}

	public void closeFiles(){
		if(!runMonteCarlo){
			this.residuals.close();
		}
	}
	public double get_filterTime(){
		return this.filterTime;
	}
	
	private void propagate(double stoptime){
		while(filterTime<stoptime){
			double dt;
			if(stoptime-filterTime > this.dtNominal)
				dt = dtNominal;
			else
				dt = stoptime-filterTime;
			
			double tnext = filterTime + dt;
			
//			if (tnext > simTime.get_sim_time()) {
//			tnext = simTime.get_sim_time();
//			System.out.println(
//			"measurement gap not an integer number of time steps");
//			}
			
			//Propagate the state forward using the process model
			double[] xnew = process.propagate(filterTime, xprev, tnext);
			
			//Get the state transition matrix for the current state
			xref = new EstSTM(xnew, this.n);
			Matrix phi = xref.phi();
			
			//Calculate the process noise matrix
			Matrix q = process.Q(tnext, this.dtNominal, xref);
			
			//Propagate the covariance matrix forward
			pnew = this.propCov(pold, phi, q);
			
			//Update the filter time and reset required variables
			filterTime = tnext;
			xref.resetPhi();
			xprev = xref.longarray();
			pold = pnew.copy();
		}		
		if(this.verbose) System.out.println("Running... "+filterTime+" / "+finalTime);//+"  range: "+new VectorN(xprev).get(0,3).mag());
	}
	
	private void process(SpacecraftModel sc, ObservationMeasurement obs){
		double y = obs.get_residual(xref.get(0,n));
		//double y = createMeasurements.mm[measNum].zPred(whichMeas,simTime,xref.get(0,n));
		
		/*Catch the case where the measurement doesn't occur*/
		//*TODO Watch this
		if( !Double.isNaN(y)) //Math.abs(y) > 0)
		{
			double r = obs.get_noise(sc);
			//double r = createMeasurements.mm[measNum].R();
			
			if(!runMonteCarlo){
//			String residualsOut = 
//				"Time:  "+simTime+"  Residual:  "+y+" Measurement Type:  "+
//				createMeasurements.measurementTypes[measNum] + " State "+whichMeas;
//			residuals.println(residualsOut);
			String residualsOut = "Time: "+filterTime+" Residual: "+y+"    Measurement Type: "+
				obs.get_measurementType()+" State "+obs.get_PRN();
			residuals.println(residualsOut);
			}
			//Use the current reference trajectory to form the H matrix
			VectorN H = obs.get_H(new VectorN(6));
			//VectorN  H = createMeasurements.mm[measNum].H(new VectorN(6));

			// compute the Kalman gain
			VectorN k = this.kalmanGain(pnew, H, r);
			
			// compute new best estimate
			VectorN xhat = k.times(y);
			

			// update state and covariance
			xref.update(xhat); 
//			if(Double.isNaN(xref.state().x[0])){
//				int donothing = 0;
//			}
			pold = this.updateCov(k, H, pnew);
		}else{
			System.err.println("Error: negative residual!"); //else visible = false;
		}
	}
	
	/** Process the measurements (using measurements from a file)
	 *  If the next measurement comes within the next timestep, propagate to the measurement
	 */
	public VectorN estimate(Time sim_Time,  SpacecraftModel sc, ObservationMeasurementList obslist, boolean measFlag) {
		
		double simTime = sim_Time.get_sim_time();
		//* catch the end of the set of measurements and return null
		//if(obs==null) measFlag=false;
		double dt = simTime-filterTime;
		
		if(measFlag){
			ObservationMeasurement obs = obslist.getCurrent();
			double measTime = Math.round(TimeUtils.days2sec*(obs.time_mjd()-sim_Time.get_epoch_mjd_utc()));
			while(obs.time_mjd()<sim_Time.get_epoch_mjd_utc() && obs!=null){ 
				obs = obslist.getNext();
				measTime=Math.round(TimeUtils.days2sec*(obs.time_mjd()-sim_Time.get_epoch_mjd_utc()));
			}
			while(measTime<filterTime && measTime>=0 && obs!=null){ 
				obs = obslist.getNext(); 
				measTime=Math.round(TimeUtils.days2sec*(obs.time_mjd()-sim_Time.get_epoch_mjd_utc()));
			}			
			while(measTime<simTime && measTime>=0){
				
				/*If necessary move to  a new time*/
				//measTime = TimeUtils.days2sec*(obs.time_mjd()-sim_Time.get_epoch_mjd_utc());
				double dt_obs = measTime-filterTime;
				//double dt_sim = simTime.get_sim_time()- filterTime;
				// detect backwards time jump
				if (dt_obs < 0.0) {
					System.out.println("backwards time jump");
					System.exit(1);
				}
				
				
				// propagate state and covariance to new time
				if (dt_obs > 0.0) {
					//while (filterTime < simTime.get_sim_time()) {
					propagate(measTime);				
				} 
				else {			
					//The time is the same so don't move the covariance
					pnew = pold.copy(); // dt = 0, no change to covariance
				}
				
				/* perform the measurement update  Currently we feed in the position that the measurement
				 * is in the state.  This is probably not used by truly scalar measurements
				 * and can be safely set to zero in those cases.
				 */	
				
				//if(initializer.parseInt(hm,"MEAS.types")!=0 && measFlag==true){
				if(measFlag==true){
					process(sc,obs);
				}
				
				// check the update
//				double zafter = obs.get_residual(xref.get(0,n));
//				//double zafter = measModel.zPred(i, t, xref.state());
//				double yafter = z - zafter;
//				process.printResiduals(simTime, y, yafter);
				
				xref.resetPhi(); // re-linearize
				//filterTime = measTime; //simTime;
				xprev = xref.longarray();
				
				obs = obslist.getNext();
				if(obs!=null)
					measTime = Math.round(TimeUtils.days2sec*(obs.time_mjd()-sim_Time.get_epoch_mjd_utc()));
				else
					measTime = -1;
			}
		}
		
		//* Propagate to simTime
		while(filterTime < simTime){
			propagate(simTime);
		}
		
		VectorN out = new VectorN(xref.get(0,n));
		
//		if(Double.isNaN(out.x[0])){
//			int donothing = 0;
//		}
		return out;
		
	}	
	
	/** Process the measurements
	 * 
	 */
	public VectorN estimate(double simTime, int measurementNum, int whichMeas, boolean measFlag) {
		
		measNum = measurementNum;
		stateNum = whichMeas;
		
		/*If necessary move to  a new time*/
		double dt = simTime- filterTime;
		
		// detect backwards time jump
		if (dt < 0.0) {
			System.out.println("backwards time jump");
			System.exit(1);
		}
		
		// propagate state and covariance to new time
		if (dt > 0.0) {
			
			while (filterTime < simTime) {
				
				double tnext = filterTime + this.dtNominal;
				
				if (tnext > simTime) {
					tnext = simTime;
					System.out.println(
					"measurement gap not an integer number of time steps");
				}
				
				//Propagate the state forward using the process model
				double[] xnew = process.propagate(filterTime, xprev, tnext);
				if(this.verbose) System.out.println("Running... "+filterTime+" / "+finalTime);
				//if(this.verbose) System.out.println("Running... "+(int)(100*filterTime/finalTime)+" %");
				
				//Get the state transition matrix for the current state
				xref = new EstSTM(xnew, this.n);
				Matrix phi = xref.phi();
				
				//Calculate the process noise matrix
				Matrix q = process.Q(tnext, this.dtNominal, xref);

				//Propagate the covariance matrix forward
				pnew = this.propCov(pold, phi, q);
				
				//Update the filter time and reset required variables
				filterTime = tnext;
				xref.resetPhi();
				xprev = xref.longarray();
				pold = pnew.copy();
				
			}
		} 
		else {
			
			//The time is the same so don't move the covariance
			pnew = pold.copy(); // dt = 0, no change to covariance
		}
	
		
		/* perform the measurement update  Currently we feed in the position that the measurement
		 * is in the state.  This is probably not used by truly scalar measurements
		 * and can be safely set to zero in those cases.
		 */
		
		
		if(initializer.parseInt(hm,"MEAS.types")!=0 && measFlag==true)
		{
			if(simTime>6680){
				int donothing = 0;
				donothing++;
			}
			double y = measurements.mm[measNum].zPred(whichMeas,simTime,xref.get(0,n));
			
			/*Catch the case where the measurement doesn't occur*/
			if( Math.abs(y) > 0)
			{
			//* TODO watch this 
				//VectorN moon = ephem.get_Geocentric_Moon_pos(TimeUtils.MJDtoJD(Time.TTtoTDB(Time.UTC2TT(mjd_epoch+filterTime/86400.0)))).times(1000);
				//VectorN r_eci = xref.get(0,3);
				//* Gets earth range
				//double dist = r_eci.mag();
				//* Gets Moon range
				double dist = 0;//(moon.minus(r_eci)).mag(); 
				
				if(!runMonteCarlo){
				if(measurements.measurementTypes[measNum].equalsIgnoreCase("y_angle_los")){
					//String residualsOut = "Time:  " + simTime +
					String residualsOut = "Time:  " + simTime +"  Dist:  " + dist +
					"  Residual:  " + (jat.math.MathUtils.RAD2DEG*y) + "  deg    Measurement Type:  " + 
					measurements.measurementTypes[measNum] + " State " + whichMeas;
					residuals.println(residualsOut);
				}else if(measurements.measurementTypes[measNum].equalsIgnoreCase("range")){
					String residualsOut = "Time:  " + simTime +"  Dist:  " + dist + 
					"  Residual:  " + (jat.math.MathUtils.RAD2DEG*y) + "  deg    Measurement Type:  " + 
					measurements.measurementTypes[measNum] + " State " + whichMeas;
					residuals.println(residualsOut);
				}else{
					String residualsOut = "Time:  " + simTime + 
					"  Residual:  " + y + " Measurement Type:  " + 
					measurements.measurementTypes[measNum] + " State " + whichMeas;
					residuals.println(residualsOut);
				}
				}
				
				//Use the current reference trajectory to form the H matrix
				VectorN  H = measurements.mm[measNum].H(new VectorN(6));

				if(!H.equals(new VectorN(n))){
				double r = measurements.mm[measNum].R();
				// compute the Kalman gain
				VectorN k = this.kalmanGain(pnew, H, r);
				
				// compute new best estimate
				VectorN xhat = k.times(y);
				
//				if(xhat.get(6,3).mag() > 15){
//					int s= 1;
//					s++;
//				}
				// update state and covariance
				
				
				xref.update(xhat); 
				
				
				y = measurements.mm[measNum].zPred(whichMeas,simTime,xref.get(0,n));

				pold = this.updateCov(k, H, pnew);
				}
			} //else visible = false;
		}
		
		// check the update
		//double zafter = measModel.zPred(i, t, xref.state());
		//double yafter = z - zafter;
		//process.printResiduals(simTime, y, yafter);
		
		xref.resetPhi(); // re-linearize
		filterTime = simTime;
		xprev = xref.longarray();
		VectorN out = new VectorN(xref.get(0,n));
		return out;
		
	}
	public void set_verbose(boolean b, double tf) {
		this.verbose = b;		
		this.finalTime = tf;
	}
	public Matrix get_pold() {
		return this.pold;
	}	

}

