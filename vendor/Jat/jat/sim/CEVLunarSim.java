/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
 * Emergent Space Technologies
 * File created by Richard C. Page III 
 **/
package jat.sim;

import jat.alg.estimators.EKF;
import jat.alg.estimators.LunarEOM;
import jat.alg.integrators.LinePrinter;
import jat.alg.integrators.RungeKutta8;
import jat.cm.Constants;
import jat.cm.TwoBody;
import jat.eph.DE405;
import jat.forces.GravitationalBody;
import jat.forces.gravity.*;
import jat.forces.density.earth.*;
import jat.forces.Moon;
import jat.forces.SolarRadiationPressure;
import jat.forces.Sun;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.measurements.ObservationMeasurementList;
import jat.measurements.OpticalMeasurementModel;
import jat.measurements.createMeasurements;
import jat.spacecraft.Spacecraft;
import jat.spacecraft.SpacecraftModel;
import jat.spacetime.EarthRef;
import jat.spacetime.LunaFixedRef;
import jat.spacetime.LunaRef;
import jat.spacetime.RSW_Frame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.spacetime.UniverseModel;
import jat.traj.RelativeTraj;
import jat.traj.Trajectory;
import jat.util.Celestia;
import jat.util.FileUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.StringTokenizer;

public class CEVLunarSim {

	//** Static Variables **//

	//* TODO Cheating
	//** Externally referenced Static Variables **//
	public static SpacecraftModel truth[];
	public static String InputFile;
	public static int JAT_case = 0;
	public static boolean JAT_runtruth = true;
	public static String MEAS_GPSSTATE,MEAS_GPS;
	public static String Truth,GEONS_Ref;
	public static boolean PlotJAT = false;
	public static boolean PlotTruth = false;
	public static boolean PlotGEONSRef = false;
	public static boolean PlotGEONSBoth = false;
	public static boolean PlotMeasurements = false;
	public static boolean Flag_GPS = false;
	public static boolean Flag_GPSState = false;
	public static boolean Flag_Cross = false;
	public static boolean COV_printoffdiag = false;

	//** Object Variables **//

	protected boolean gravityModel;
	protected HashMap input;
	protected createMeasurements created_meas;

	//private SpacecraftModel truth[];
	protected SpacecraftModel ref[];
	protected Trajectory truth_traj[];
	protected Trajectory ref_traj[];
	protected Trajectory sim_truth;
	protected Trajectory sim_cov;
	protected static int numSpacecraft;
	//public SimModel[] truth = null;
	//public SimModel[] ref   = null;
	protected FileOutputStream[] trajectories;
	protected FileOutputStream[] truths;
	protected FileOutputStream[] ECIError;
	protected FileOutputStream[] covariance;

	protected int simStep;
	public EKF filter;
	protected double dt;
	private VectorN newState;
	private int numStates;
	//public double simTime;
	protected Time simTime;

	protected boolean verbose_estimation=false;
	protected boolean useMeas=true;
	private boolean obsFromFile=false;
	private boolean verbose_timestep = false;
	private double tf;
	private RungeKutta8 rk8;
	private boolean runMonteCarlo = false;

	public static String JAT_name;
	private static boolean PlotBoeing;
	private static boolean PlotCov;
	private static boolean MONTE_Plot_override = false;

	public CEVLunarSim(boolean useFilter) {
		//super(useFilter);
		this.useMeas = useFilter;
		initializeConst();
	}

//	** Object Methods **//
	protected void initializeConst(){		
		String fs, dir_in;
		fs = FileUtil.file_separator();
		try{
			dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"input"+fs;
		}catch(Exception e){
			dir_in = "";
		}
		this.input = initializer.parse_file(dir_in+InputFile);
		runMonteCarlo = initializer.parseBool(input, "init.runMonteCarlo");
		double MJD0 =  initializer.parseDouble(input,"init.MJD0");
		double T0 = initializer.parseDouble(input, "init.T0");
		double MJDF =  initializer.parseDouble(input, "init.MJDF");
		double TF = initializer.parseDouble(input, "init.TF");
		simTime = new Time(MJD0+T0/86400);

		//geons_truth = parseGEONSTruth(simTime.get_epoch_mjd_utc(),MJDF+TF/86400);

		numSpacecraft = initializer.parseInt(input,"prop.NumSpacecraft");
		numStates = initializer.parseInt(input,"FILTER.states");
		dt = initializer.parseInt(input,"init.dt");
		if(JAT_runtruth){
			truth = new SpacecraftModel[numSpacecraft];
			truth_traj = new Trajectory[2];
		}
		ref   = new SpacecraftModel[numSpacecraft];
		ref_traj = new Trajectory[2];
		for(int i=0; i<2; i++){
			if(JAT_runtruth) truth_traj[i] = new Trajectory();
			ref_traj[i] = new Trajectory();
		}
		sim_cov = new Trajectory();
		int i=0;
		double[] sigmas = new double[6];
		String tmp = "P0."+i+".X";
		sigmas[6*i + 0] = initializer.parseDouble(input,tmp);
		tmp = "P0."+i+".Y";
		sigmas[6*i + 1] = initializer.parseDouble(input,tmp);
		tmp = "P0."+i+".Z";
		sigmas[6*i + 2] = initializer.parseDouble(input,tmp);
		tmp = "P0."+i+".VX";
		sigmas[6*i + 3] = initializer.parseDouble(input,tmp);
		tmp = "P0."+i+".VY";
		sigmas[6*i + 4] = initializer.parseDouble(input,tmp);
		tmp = "P0."+i+".VZ";
		sigmas[6*i + 5] = initializer.parseDouble(input,tmp);
		sim_cov.add(simTime.get_sim_time(), sigmas);

		created_meas = new createMeasurements(input);
		filter = new EKF(input,CEVLunarSim.JAT_case,created_meas);
		rk8 = new RungeKutta8(dt);

	}

	protected void initialize()
	{
		double[] r = new double[3];
		double[] tr = new double[3];  //variable for true trajectory
		double[] v = new double[3];
		double[] tv = new double[3];  //variable for true trajectory
		double cr,cd,area,mass,dt;

		double MJD0 =  initializer.parseDouble(input,"init.MJD0") + initializer.parseDouble(input,"init.T0")/86400.0;
		simTime = new Time(MJD0);
		double MJDF =  initializer.parseDouble(input,"init.MJDF") + initializer.parseDouble(input,"init.TF")/86400.0;
		this.tf = (MJDF-MJD0)*86400.0;
		//For each spacecraft extract the initial vector and force information
		//Use this information to create a sim model
		System.out.println("Propagating "+numSpacecraft+" Spacecraft");
		for(int i = 0;i < numSpacecraft; i++)
		{	
			/*Position*/
			String refs = "REF_STATE.";
			String tru = "TRUE_STATE.";

			String str  = refs+i+".X";
			String strt = tru+i+".X";
			r[0] = initializer.parseDouble(this.input,str);
			tr[0] = initializer.parseDouble(this.input,strt);

			str  = refs+i+".Y";
			strt = tru+i+".Y";
			r[1] = initializer.parseDouble(this.input,str);
			tr[1] = initializer.parseDouble(this.input,strt);

			str  = refs+i+".Z";
			strt = tru+i+".Z";
			r[2] = initializer.parseDouble(this.input,str);
			tr[2] = initializer.parseDouble(this.input,strt);

			/*Velocity*/
			str  = refs+i+".VX";
			strt = tru+i+".VX";
			v[0] = initializer.parseDouble(this.input,str);
			tv[0] = initializer.parseDouble(this.input,strt);

			str  = refs+i+".VY";
			strt = tru+i+".VY";
			v[1] = initializer.parseDouble(this.input,str);
			tv[1] = initializer.parseDouble(this.input,strt);

			str  = refs+i+".VZ";
			strt = tru+i+".VZ";
			v[2] = initializer.parseDouble(this.input,str);
			tv[2] = initializer.parseDouble(this.input,strt);


			/*Solar Radiation Pressure Coefficient*/
			str = "jat."+i+".Cr";
			cr   = initializer.parseDouble(this.input,str);

			/*Drag Coefficient*/
			str = "jat."+i+".Cd";
			cd   = initializer.parseDouble(this.input,str);

			/*Initial Mass*/
			str = "jat."+i+".mass";
			mass = initializer.parseDouble(this.input,str);

			/*Initial Area*/
			str = "jat."+i+".area";
			area = initializer.parseDouble(this.input,str);
			
			LunaFixedRef lfr = new LunaFixedRef();
			LunaRef lref = new LunaRef();
			ReferenceFrameTranslater trans = new ReferenceFrameTranslater(lfr,lref,simTime);
			/*Read in the appropriate model flags*/			
			VectorN rr = new VectorN(r);
			VectorN vv = new VectorN(v);
//			DE405 ephem = new DE405();
//			VectorN e2m = ephem.get_Geocentric_Moon_pos(TimeUtils.MJDtoJD(MJD0));
			//rr = (rr.plus(e2m));
			rr = rr.times(1000);
			vv = vv.times(1000);			
			vv = trans.translateVelocity(vv,rr);
			rr = trans.translatePoint(rr);

//			VectorN out = new VectorN(6);
//			out.x[0] = 0.01721537446918e7;
//			out.x[1] = 0.16362908568758e7;
//			out.x[2] = -0.08151111642821e7;
//			out.x[3] = 0.13822834993145e7;
//			out.x[4] = 0.71321929326561e7;
//			out.x[5] = 1.46094273856495e7;
//			rr = out.get(0,3);
//			vv = out.get(3,3);
			
			Spacecraft s = new Spacecraft(rr,vv,cr,cd,area,mass);
			s.set_use_params_in_state(false);
			int ncoef = initializer.parseInt(input,"jat.0.lunar_n.jat");
			LunarEOM eom = new LunarEOM(input,ncoef,false);
			ref[i] = new SpacecraftModel(s,eom,MJD0);

			rr = new VectorN(tr);
			vv = new VectorN(tv);
			//rr = (rr.plus(e2m));
			rr = rr.times(1000);
			vv = vv.times(1000);		
			TwoBody orbit2 = new TwoBody(LunaRef.GM_Luna,rr,vv);
			orbit2.printElements("Lunar LCF");

			vv = trans.translateVelocity(vv,rr);
			rr = trans.translatePoint(rr);
			s = new Spacecraft(rr,vv,cr,cd,area,mass);
			s.set_use_params_in_state(false);
			truth[i] = new SpacecraftModel(s,eom,MJD0);			

			TwoBody orbit = new TwoBody(LunaRef.GM_Luna,rr,vv);
			orbit.printElements("Lunar");


			/*Set the step size for the trajectory generation*/
			/*Set the integrator Step size*/
			dt = initializer.parseInt(this.input,"init.dt");
			//ref[i].set_sc_dt(dt);
			//truth[i].set_sc_dt(dt);
			int numSats = i;
			VectorN newState = new VectorN(filter.xref.state().get(0, 6));
			double tmpState[] = new double[newState.length];
				//Extract the state of the current satellite
				for(int k = 0;k < newState.length; k++)
				{
					tmpState[k]=newState.x[numSats*6 + k];
				}
				ref[numSats].get_spacecraft().updateMotion(tmpState);	
				ref[numSats].update(simTime.get_sim_time());
//				Write out the current True States
				double [] true_state = truth[numSats].get_spacecraft().toStateVector();

				//Print out simStep + 1 because we didn't output the initial state
				VectorN vecTime =  new VectorN(1,(simStep)*dt);
				VectorN trueState = new VectorN(true_state);
				VectorN truthOut = new VectorN(vecTime,trueState);
				if(!runMonteCarlo)
					new PrintStream(truths[numSats]).println (truthOut.toString());

				//Write out the current State estimates
				double [] ref_state  = tmpState;//ref[numSats].get_spacecraft().toStateVector();
				VectorN vecState = new VectorN(ref_state);
				VectorN stateOut = new VectorN(vecTime,vecState);
				if(!runMonteCarlo)
					new PrintStream(trajectories[numSats]).println (stateOut.toString());

				//Output the current ECI error
				VectorN error_out = new VectorN(6);
				for(int k = 0; k < 6; k++)
				{
					error_out.x[k] = true_state[k] - ref_state[k];
				}
				VectorN ErrState = new VectorN(error_out);
				stateOut = new VectorN(vecTime,ErrState);
				new PrintStream(ECIError[numSats]).println (stateOut.toString());

//				Output the current Covariances
//				Matrix Covariance = EKF.pold;
				Matrix Covariance = filter.get_pold();
				VectorN var = new VectorN(6);
				if(!runMonteCarlo || MONTE_Plot_override){
					RSW_Frame rictrans = new RSW_Frame(trueState.get(0, 3),trueState.get(3, 3));
					Matrix covtrans = rictrans.ECI2RSW();
					Matrix covr =   covtrans.times((Covariance.getMatrix(0, 2, 0, 2)).times(covtrans.transpose()));
					Matrix covv =   covtrans.times((Covariance.getMatrix(3, 5, 3, 5)).times(covtrans.transpose()));
					var = new VectorN(covr.diagonal(),covv.diagonal());
				}
				int numStates = filter.get_numStates();
				//double[] tmp = new double[numStates*numStates];
				double[] tmp = new double[numStates];
				double[] tmp2 = new double[numStates];
				int m = 0;
				for(int k = 0; k < numStates; k++)
				{
					//for(int j = 0; j < numStates; j++)
					//{
					tmp[m] = Covariance.get(k,k);
					try{
						if(!runMonteCarlo || MONTE_Plot_override)
							tmp2[m] = Math.sqrt(var.get(k));
					}catch(Exception e){
						//out of range
					}
					m++;
					//}
				}			
				if(!runMonteCarlo || MONTE_Plot_override)
					sim_cov.add(simTime.get_sim_time()/3600.0, new VectorN(tmp2,6).x);
				VectorN ErrCov = new VectorN(tmp);
				stateOut = new VectorN(vecTime,ErrCov);
				new PrintStream(covariance[numSats]).println (stateOut.toString());
			}
	
	}

	/**
	 * Propagation method for an individual spacecraft.  Increments the 
	 * model held in the spacecraft flight computer by a time 'dt'.
	 * Updates the computer's models according to the progression of time.
	 * Updates the spacecraft state.
	 * 
	 */
	public void step(SpacecraftModel sm){
		double t = sm.get_sc_t();		
		if(verbose_timestep ){
			System.out.println("step: "+t+" / "+tf+"    stepsize: "+dt);
		}

		//rk8.setStepSize(sm.get_sc_dt());
		//rk8.setStepSize(dt);
		//* update models
		//double mjd_utc = spacetime.get_mjd_utc();
		double[] X = new double[6];
		double[] Xnew = new double[6];
		//double[] thrust = new double[3];
		//VectorN rnew;
		//VectorN vnew;
		//double[] tmp = new double[6];
		double num_sc = 1;
		for(int i=0; i<num_sc; i++){

			Spacecraft s = sm.get_spacecraft();
			X = s.toStateVector(false);
			Xnew = rk8.step(t, X, sm);
			//* store new values
			//rnew = new VectorN(Xnew[0],Xnew[1],Xnew[2]);
			//vnew = new VectorN(Xnew[3],Xnew[4],Xnew[5]);
			s.updateState(Xnew,false);
		}
		//* update simulation time
//		if(t > (tf-sm.get_sc_dt()) && t != tf){
//		sm.set_sc_dt(tf-t);
//		}
//		t=t+sm.get_sc_dt();
		if(t > (tf-dt) && t != tf){
			dt=(tf-t);
		}
		t=t+dt;
		//* add to trajectory
		//traj.add(sm.get_sc_mjd_utc(),sm.get_spacecraft().toStateVector());
		//* update the universe
		sm.update(t);
		//iteration++;
	}


	/**
	 * 
	 */
	protected void filter()
	{

		/*Provide the algorithm to run the Kalman filter one step at
		 a time.  This may need to be placed into another file later on*/

		/*Determine the num ber of measurements that we will be processing.
		 * Since this is a scalar update, we can simply loop over the 
		 * measurement update step that many times
		 */
		int numMeas = created_meas.getNumberMeasurements(); 

		/*Loop over the number of measurements being carefull to 
		 * omit measurements of 0
		 */

		VectorN newState = new VectorN(ref[0].get_spacecraft().toStateVector());
		int processedMeasurements = 0;
		for(int i = 0;i<numMeas;i++)
		{			
			if(simTime.get_sim_time()%(created_meas.frequency[i]) ==0 )
			{
				//String tmp = "MEAS."+i+".satellite";
				//int sat = initializer.parseInt(this.input,tmp);

				if(((OpticalMeasurementModel)created_meas.mm[i]).get_type()==OpticalMeasurementModel.TYPE_YANGLE_LOS){
					newState = filter.estimate(simTime.get_sim_time(),i,0,useMeas);
					processedMeasurements ++;
					newState = filter.estimate(simTime.get_sim_time(),i,1,this.useMeas);
					processedMeasurements ++;
					//System.out.println("Processing Measurement at time: " + simTime.get_sim_time());
				} else if(((OpticalMeasurementModel)created_meas.mm[i]).get_type()==OpticalMeasurementModel.TYPE_LANDMARK){
					OpticalMeasurementModel opt = (OpticalMeasurementModel)created_meas.mm[i];						
					int num_landmarks = opt.get_num_landmarks();
					for(int nmark=0; nmark < num_landmarks; nmark++){
						opt.currentLandmark(nmark);
						newState = filter.estimate(simTime.get_sim_time(),i,0,useMeas);
						processedMeasurements ++;
						newState = filter.estimate(simTime.get_sim_time(),i,1,useMeas);
						processedMeasurements ++;
					}						
				} else {
					newState = filter.estimate(simTime.get_sim_time(),i,0,this.useMeas);
					processedMeasurements ++;
				}

			}	

		}

		//catch the case where there are no measurements, set the measurement
		//flag to false to tell the filter to just propagate
		if(processedMeasurements  == 0)
			newState = filter.estimate(simTime.get_sim_time(), 0,0,false);


		//Update the current state with the output of the filter
		//Write the current state information to files

		double tmpState[] = new double[newState.length];
		for(int numSats = 0; numSats < numSpacecraft; numSats ++)
		{
			//Extract the state of the current satellite
			for(int i = 0;i < newState.length; i++)
			{
				tmpState[i]=newState.x[numSats*6 + i];
			}
			ref[numSats].get_spacecraft().updateMotion(tmpState);	
			ref[numSats].update(simTime.get_sim_time());
//			Write out the current True States
			double [] true_state = truth[numSats].get_spacecraft().toStateVector();

			//Print out simStep + 1 because we didn't output the initial state
			VectorN vecTime =  new VectorN(1,(simStep)*dt);
			VectorN trueState = new VectorN(true_state);
			VectorN truthOut = new VectorN(vecTime,trueState);
			if(!runMonteCarlo)
				new PrintStream(truths[numSats]).println (truthOut.toString());

			//Write out the current State estimates
			double [] ref_state  = tmpState;//ref[numSats].get_spacecraft().toStateVector();
			VectorN vecState = new VectorN(ref_state);
			VectorN stateOut = new VectorN(vecTime,vecState);
			if(!runMonteCarlo)
				new PrintStream(trajectories[numSats]).println (stateOut.toString());

			//Output the current ECI error
			VectorN error_out = new VectorN(6);
			for(int i = 0; i < 6; i++)
			{
				error_out.x[i] = true_state[i] - ref_state[i];
			}
			VectorN ErrState = new VectorN(error_out);
			stateOut = new VectorN(vecTime,ErrState);
			new PrintStream(ECIError[numSats]).println (stateOut.toString());

//			Output the current Covariances
//			Matrix Covariance = EKF.pold;
			Matrix Covariance = filter.get_pold();
			VectorN var = new VectorN(6);
			if(!runMonteCarlo || MONTE_Plot_override){
				RSW_Frame trans = new RSW_Frame(trueState.get(0, 3),trueState.get(3, 3));
				Matrix covtrans = trans.ECI2RSW();
				Matrix covr =   covtrans.times((Covariance.getMatrix(0, 2, 0, 2)).times(covtrans.transpose()));
				Matrix covv =   covtrans.times((Covariance.getMatrix(3, 5, 3, 5)).times(covtrans.transpose()));
				var = new VectorN(covr.diagonal(),covv.diagonal());
			}
			int numStates = filter.get_numStates();
			//double[] tmp = new double[numStates*numStates];
			double[] tmp = new double[numStates];
			double[] tmp2 = new double[numStates];
			int k = 0;
			for(int i = 0; i < numStates; i++)
			{
				//for(int j = 0; j < numStates; j++)
				//{
				tmp[k] = Covariance.get(i,i);
				try{
					if(!runMonteCarlo || MONTE_Plot_override)
						tmp2[k] = Math.sqrt(var.get(i));
				}catch(Exception e){
					//out of range
				}
				k++;
				//}
			}			
			if(!runMonteCarlo || MONTE_Plot_override)
				sim_cov.add(simTime.get_sim_time()/3600.0, new VectorN(tmp2,6).x);
			VectorN ErrCov = new VectorN(tmp);
			stateOut = new VectorN(vecTime,ErrCov);
			new PrintStream(covariance[numSats]).println (stateOut.toString());
		}

	}
	protected void openFiles()
	{
		/*The number and types of files that are created are based
		 * upon the number of spacecraft and the simulation mode*/
		if(!runMonteCarlo){
			trajectories = new FileOutputStream [numSpacecraft];
			if(JAT_runtruth){
				truths       = new FileOutputStream [numSpacecraft]; 
			}
		}
		ECIError     = new FileOutputStream [numSpacecraft];
		covariance   = new FileOutputStream [numSpacecraft];


		String fs, dir_in;
		fs = FileUtil.file_separator();
		try{
			dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
		}catch(Exception e){
			dir_in = "";
		}

		//String fileName5 = dir_in+"Visible.txt";
		/*	try {
		 visableSats = new FileOutputStream(fileName5);
		 } catch (FileNotFoundException e1) {
		  e1.printStackTrace();
		  }*/
		for(int numSats = 0; numSats < numSpacecraft; numSats++)
		{
			try
			{
				// Open an output stream
//				String fileName = dir_in+"Sat"+numSats+"ECI.txt";
//				String fileName2 = dir_in+"TRUE"+numSats+"ECI.txt";
//				String fileName3 = dir_in+"Error"+numSats+"ECI.txt";
//				String fileName4 = dir_in+"Covariance"+numSats+"ECI.txt";
				String fileName = dir_in+"JAT"+JAT_case+"LCI.eci";
				String fileName3 = dir_in+"JAT"+JAT_case+"LCI.err";
				String fileName4 = dir_in+"JAT"+JAT_case+"LCI.cov";
				if(!runMonteCarlo){
					trajectories[numSats] = new FileOutputStream (fileName);
					if(JAT_runtruth){
						String fileName2 = dir_in+"JAT_true_"+JAT_case+".txt";
						truths[numSats]       = new FileOutputStream (fileName2);
					}
				}
				ECIError[numSats]     = new FileOutputStream(fileName3);
				covariance[numSats]   = new FileOutputStream(fileName4);
			}
			catch (IOException e)
			{
				e.printStackTrace();
				System.err.println ("Unable to write to file");
				System.exit(-1);
			}
		}
	}
	protected void closeFiles()
	{
		/*The number and types of files that are created are based
		 * upon the number of spacecraft and the simulation mode*/

		for(int numFiles = 0; numFiles < numSpacecraft; numFiles++)
		{
			{
				// Close the files
				try 
				{
					if(!runMonteCarlo){
						trajectories[numFiles].close();					
						if(JAT_runtruth){
							truths[numFiles].close();
						}
					}
					ECIError[numFiles].close();
					covariance[numFiles].close();
				} 
				catch (IOException e) 
				{
					e.printStackTrace();
				}	
			}	
			//EKF.residuals.close();
			filter.closeFiles();
			/*try {
			 //visableSats.close();
			  } catch (IOException e) {
			   e.printStackTrace();
			   }*/
		}

		System.out.println("Closing files and exiting . . ");
	}
	protected void propagate(double simStep)
	{
		/*In the propagation step we want to move the orbits of all the
		 satellites forward the desired timestep.  Output the states of 
		 each vehicle*/

		for(int numSats = 0; numSats < numSpacecraft; numSats ++)
		{
			//Step the trajectory forward the indicated amount

			//step(truth[numSats],truth_traj[numSats]);
			step(truth[numSats]); 

			//Extract the State Vector
			//Only need to propagate the reference state if 
			//We aren't estimating it
			//* NOTE * originally == 0  ???
			//* TODO Add flag to avoid running the 'truth' trajectory
			if(initializer.parseDouble(this.input,"init.mode") == 0)
			{
				double [] true_state = truth[numSats].get_spacecraft().toStateVector();

				//Print out simStep + 1 because we didn't output the initial state
				VectorN vecTime =  new VectorN(1,(simStep)*dt);
				VectorN trueState = new VectorN(true_state);
				VectorN truthOut = new VectorN(vecTime,trueState);
				new PrintStream(truths[numSats]).println (truthOut.toString());


				//step(ref[numSats],ref_traj[numSats]);
				step(ref[numSats]);
				double [] ref_state  = ref[numSats].get_spacecraft().toStateVector();
				VectorN vecState = new VectorN(ref_state);
				VectorN stateOut = new VectorN(vecTime,vecState);
				new PrintStream(trajectories[numSats]).println (stateOut.toString());
			}

		}

	}

	public void runloop(){
		//Main loop for Filtering Routine.  This includes a Kalman filter 
		//and, trajectory generation, measurement generation as well as
		//guidance and control hooks.
		double start = System.currentTimeMillis();

		//Open files for saving off the generated data.  
		openFiles();


		//Initialize
		//This step should read in all the variables and set up all the required
		//models and functionality accordingly.  One should try to be careful
		//to flag any inconsistancies in the input data to save crashes later on.

		initialize();

		//for(int i=0; i<numSpacecraft; i++){
		//truth_traj[i].add(truth[i].get_sc_t(),truth[i].get_spacecraft().toStateVector());
		//ref_traj[i].add(ref[i].get_sc_t(),ref[i].get_spacecraft().toStateVector());
		LunaRef lref = new LunaRef();
		EarthRef eref = new EarthRef(simTime);
		LunaFixedRef lfr = new LunaFixedRef();
		//ReferenceFrameTranslater trans = new ReferenceFrameTranslater(lref,eref,simTime);
		ReferenceFrameTranslater trans = new ReferenceFrameTranslater(lref,lfr,simTime);
		VectorN xt = new VectorN(truth[0].get_spacecraft().toStateVector());
		VectorN xr = new VectorN(ref[0].get_spacecraft().toStateVector());
		VectorN vt = trans.translateVelocity(xt.get(3,3),xt.get(0,3));
		VectorN vr = trans.translateVelocity(xr.get(3,3),xr.get(0,3));			
		VectorN rt = trans.translatePoint(xt.get(0,3));
		VectorN rr = trans.translatePoint(xr.get(0,3));
		truth_traj[0].add(simTime.mjd_utc(),new VectorN(rt,vt).x);
		ref_traj[0].add(simTime.mjd_utc(),new VectorN(rr,vr).x);
		truth_traj[1].add(simTime.get_sim_time()/3600.0,xt.x);
		ref_traj[1].add(simTime.get_sim_time()/3600.0,xr.x);
		//}

		/*Cache off the simulation mode */
		int filterMode = initializer.parseInt(this.input,"init.mode");

		//Compute the length of the simulation in seconds
		double MJD0 =  initializer.parseDouble(this.input,"init.MJD0");
		double MJDF =  initializer.parseDouble(this.input,"init.MJDF");
		double T0   =  initializer.parseDouble(this.input,"init.T0");
		double TF   =  initializer.parseDouble(this.input,"init.TF");
		double simLength = Math.round((MJDF - MJD0)*86400 + TF - T0);

		MJD0 = MJD0 + T0/86400.0;
		MJDF = MJDF + TF/86400.0;
		//simTime = 0; //* this is done in call to "initialize()"
		simTime = new Time(MJD0);

		this.tf = simLength;
		set_verbose(this.verbose_estimation);
		//if(!Flag_GPS && !Flag_GPSState && !Flag_Cross ) this.useMeas = false;
		//double simLength = Math.round(TF - T0);
		//ObservationMeasurement obs = obs_list.getFirst();

		for( simStep = 1; simStep < simLength/dt; simStep ++)
		{
			//if(this.verbose_estimation) 
			//System.out.println("running..."+(dt*simStep)+" / "+simLength);
			//if(simStep%100 == 0)
			//	System.out.println(simStep*5);

			//simTime = simStep*dt;			
			propagate(simStep*dt);
			simTime.update(simStep*dt);

			filter();

//			if(Double.isNaN(ref[0].get_spacecraft().toStateVector()[0])){// || simTime.get_sim_time()>4620){
//			int donothing = 0;
//			donothing++;
//			}
			//System.out.println("SimTime: " + simTime.get_sim_time() + " SimStep: " + simStep);

			//for(int i=0; i<numSpacecraft; i++){
			//lref = new LunaRef();
			eref = new EarthRef(simTime);
			//trans = new ReferenceFrameTranslater(lref,eref,simTime);
			trans = new ReferenceFrameTranslater(lref,lfr,simTime);
			xt = new VectorN(truth[0].get_spacecraft().toStateVector());
			xr = new VectorN(ref[0].get_spacecraft().toStateVector());
			vt = trans.translateVelocity(xt.get(3,3),xt.get(0,3));
			vr = trans.translateVelocity(xr.get(3,3),xr.get(0,3));			
			rt = trans.translatePoint(xt.get(0,3));
			rr = trans.translatePoint(xr.get(0,3));
			truth_traj[0].add(simTime.mjd_utc(),new VectorN(rt,vt).x);
			ref_traj[0].add(simTime.mjd_utc(),new VectorN(rr,vr).x);			
			truth_traj[1].add(simTime.get_sim_time()/3600.0,xt.x);
			ref_traj[1].add(simTime.get_sim_time()/3600.0,xr.x);
			//}			
		}

		/*Close all output files*/
		closeFiles();
		System.gc();

		double elapsed = (System.currentTimeMillis()-start)*0.001/60;
		System.out.println("Elapsed time [min]: "+elapsed);

		/* Post Processing */
		if(!runMonteCarlo || MONTE_Plot_override){
			LinePrinter lp = new LinePrinter();
			RelativeTraj[] reltraj = new RelativeTraj[3];
			Trajectory[] boeing = new Trajectory[2];
			if(PlotBoeing){
				boeing = parseBoeing();
			}
//			LinePrinter traj_lp = new LinePrinter("C:/Code/Misc/truth_LCF.txt");
//			truth_traj[0].printAll(traj_lp);
//			traj_lp.close();
//			LinePrinter traj_lp2 = new LinePrinter("C:/Code/Misc/truth_LCI.txt");
//			truth_traj[1].printAll(traj_lp2);
//			traj_lp2.close();

			double mismatch_tol = 0.00001;
			//* TODO Plot marker
			if(PlotJAT){
				reltraj[0] = new RelativeTraj(ref_traj[1],truth_traj[1],lp,""+JAT_case+" Jat(Ref) v Jat(Truth)");
				reltraj[0].setVerbose(false);
				if(PlotCov)
					reltraj[0].process(mismatch_tol,sim_cov);
				else
					reltraj[0].process(mismatch_tol);
				reltraj[0].process_RSS(mismatch_tol);
				//reltraj[i].process_ECI(mismatch_tol);
//				reltraj[1] = new RelativeTraj(ref_traj[1],boeing[1],lp,""+JAT_case+" Jat (Ref) v Boeing");
//				reltraj[1].setVerbose(false);
//				reltraj[1].process(mismatch_tol);
			}
			if(PlotBoeing){
				reltraj[1] = new RelativeTraj(truth_traj[1],boeing[1],lp,""+JAT_case+" Jat (Truth) v Boeing");
				reltraj[1].setVerbose(false);
				reltraj[1].process(mismatch_tol);
			}

			try {
				Celestia cel = new Celestia("C:/Code/Celestia/");
				cel.set_trajectory_meters(ref_traj[0],MJD0);
				cel.write_trajectory("jat_ref_"+JAT_name+JAT_case,"jat_ref_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(MJD0));
				cel.set_trajectory_meters(truth_traj[0],MJD0);
				cel.write_trajectory("jat_truth_"+JAT_name+JAT_case,"jat_truth_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(MJD0));
				if(PlotBoeing){
					cel.set_trajectory_meters(boeing[0],MJD0);
					cel.write_trajectory("jat_boeing_"+JAT_name+JAT_case,"jat_boeing_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(MJD0));
				}
//				cel.set_trajectory_meters(truth_traj[1],MJD0);				
//				cel.write_trajectory("jat_lunar_"+JAT_name+JAT_case,"jat_lunar_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(MJD0));
			} catch (IOException e) {
				//e.printStackTrace();
				System.err.println("Couldn't write to Celestia.");
			}
		}
	}

	public void set_verbose(boolean b){
		this.verbose_estimation = b;
		if(filter != null)
			filter.set_verbose(b,this.tf);
	}

	public Trajectory[] parseBoeing(){
		Trajectory[] traj = new Trajectory[2];
		traj[0] = new Trajectory();
		traj[1] = new Trajectory();
		String file = "C:/Code/Jat/jat/sim/input/llo_ascii.txt";
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(file)));
			String line = "%init";
			for(int i=0; i<4; i++) line = in.readLine();
			StringTokenizer tok;
			double mjd = 58232;
			Time t = new Time(mjd);
			for(int i=0; i<4000; i++){
				line = in.readLine();
				tok = new StringTokenizer(line," ");
				double x[] = new double[6];
				double tmp = Double.parseDouble(tok.nextToken());
				t.update(tmp);
				for(int j=0; j<6; j++){
					x[j] = Double.parseDouble(tok.nextToken());
				}
				VectorN state = new VectorN(x);
				LunaFixedRef lfr = new LunaFixedRef();
				LunaRef lref = new LunaRef();
				ReferenceFrameTranslater trans = new ReferenceFrameTranslater(lfr,lref,t);
				VectorN r = trans.translatePoint(state.get(0,3).times(1000));
				VectorN v = trans.translateVelocity(state.get(3, 3).times(1000), state.get(0, 3).times(1000));
				traj[1].add(t.mjd_utc(), new VectorN(r,v).x);
				EarthRef eref = new EarthRef(t);
				trans = new ReferenceFrameTranslater(lfr,eref,t);
				r = trans.translatePoint(state.get(0,3).times(1000));
				v = trans.translateVelocity(state.get(3, 3).times(1000), state.get(0, 3).times(1000));
				traj[0].add(t.mjd_utc(), new VectorN(r,v).x);
			}
			return traj;
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e){

		}
		return traj;

	}
//	** Main **//

	public static void main(String[] args) {

		boolean useFilter = true;
		CEVLunarSim.MONTE_Plot_override  = true;
		CEVLunarSim.PlotJAT = true;
		CEVLunarSim.PlotBoeing = false;
		CEVLunarSim.PlotCov = true;
		double mc_start = System.currentTimeMillis();

		for(int mc = 3; mc< 4; mc++){
			//* TODO Flag marker			

			switch(mc){
			case 0:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_3KLM_001.txt";
				break;
			case 1:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_3KLM_01.txt";
				break;
			case 2:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_3KLM_1.txt";
				break;
			case 3:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_1ULMwD_001.txt";
				break;
			case 4:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_1ULMwD_01.txt";
				break;
			case 5:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_1ULMwD_1.txt";
				break;
			default:
				CEVLunarSim.InputFile = "initialConditions_cev_llo_3KLM_001.txt";
			break;
			}
			if (args.length > 0) {
				CEVLunarSim.InputFile = args[0];
			}

			String fs, dir_in;
			fs = FileUtil.file_separator();
			try{
				dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"input"+fs;
			}catch(Exception e){
				dir_in = "";
			}
			HashMap hm = initializer.parse_file(dir_in + CEVLunarSim.InputFile);
			boolean runningMonteCarlo = initializer.parseBool(hm,"init.runMonteCarlo");
			int num_runs = (runningMonteCarlo ? 
					initializer.parseInt(hm,"MONTE.num_runs") : 1);
			Integer startCase = initializer.parseInt(hm, "MONTE.start_case_num");
			if (runningMonteCarlo && (startCase != null)) {
				CEVLunarSim.JAT_case = startCase;
			}
			for(int c=0; c<num_runs; c++){

				CEVLunarSim.JAT_name = "lowlunar";

				//CEVLunarSim.InputFile = "initialConditions_cev_llo.txt";

				CEVLunarSim Sim = new CEVLunarSim(useFilter);
				Sim.set_verbose(true);
				Sim.runloop();

				try{
					OpticalMeasurementModel.fobs.close();
					OpticalMeasurementModel.fpred.close();
				}catch(NullPointerException ne){

				}
				JAT_case++;
			}
		}

		double mc_elapsed = (System.currentTimeMillis()-mc_start)*0.001/60;
		System.out.println("MonteCarlo time [min]: "+mc_elapsed);
		
		System.out.println("Finished.");
	}
}
