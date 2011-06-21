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
import jat.alg.integrators.LinePrinter;
import jat.cm.Constants;
import jat.forces.GravitationalBody;
import jat.forces.gravity.*;
import jat.forces.density.earth.*;
import jat.forces.Moon;
import jat.forces.SolarRadiationPressure;
import jat.forces.Sun;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.measurements.OpticalMeasurementModel;
import jat.measurements.createMeasurements;
import jat.spacecraft.Spacecraft;
import jat.spacecraft.SpacecraftModel;
import jat.spacetime.LunaFixedRef;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.spacetime.UniverseModel;
import jat.traj.RelativeTraj;
import jat.traj.Trajectory;
import jat.util.Celestia;
import jat.util.FileUtil;
import jat.forces.gravity.earth.*;

import java.io.IOException;
import java.io.PrintStream;

public class CEVSim extends EstimatorSimModel {

//** Static Variables **//
	
	public static String JAT_name;
	//** Object Variables **//
	private int numStates;
	
	public CEVSim(){
		super();
	}
	
	public CEVSim(boolean useFilter) {
		super(useFilter);
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
		
		double MJD0 =  initializer.parseDouble(input,"init.MJD0");
		double T0 = initializer.parseDouble(input, "init.T0");
		double MJDF =  initializer.parseDouble(input, "init.MJDF");
		double TF = initializer.parseDouble(input, "init.TF");
		this.mjd_utc_start = MJD0+T0/86400;
		simTime = new Time(MJD0+T0/86400);
		
		//geons_truth = parseGEONSTruth(simTime.get_epoch_mjd_utc(),MJDF+TF/86400);
			
		numSpacecraft = initializer.parseInt(input,"prop.NumSpacecraft");
		numStates = initializer.parseInt(input,"FILTER.states");
		dt = initializer.parseInt(input,"init.dt");
		if(JAT_runtruth){
			truth = new SpacecraftModel[numSpacecraft];
			truth_traj = new Trajectory[numSpacecraft];
		}
		ref   = new SpacecraftModel[numSpacecraft];
		ref_traj = new Trajectory[numSpacecraft];
		for(int i=0; i<numSpacecraft; i++){
			if(JAT_runtruth) truth_traj[i] = new Trajectory();
			ref_traj[i] = new Trajectory();
		}
		created_meas = new createMeasurements(input);
		filter = new EKF(input);
		
	}

	protected void initialize()
	{
		double[] r = new double[3];
		double[] tr = new double[3];  //variable for true trajectory
		double[] v = new double[3];
		double[] tv = new double[3];  //variable for true trajectory
		double cr,cd,area,mass,dt;
		
		double MJD0 =  initializer.parseDouble(input,"init.MJD0")+initializer.parseDouble(input, "init.T0")/86400.0;
		
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
						
			/*Read in the appropriate model flags*/
			boolean[] force_flag = createForceFlag(i); 
			VectorN rr = new VectorN(r);
			VectorN vv = new VectorN(v);
			Spacecraft s = new Spacecraft(rr,vv,cr,cd,area,mass);
			s.set_use_params_in_state(false);
			UniverseModel spacetime = createUniverseModel(MJD0,s,force_flag, gravityModel, "HP");
			spacetime.set_use_iers(false);
			ref[i] = new SpacecraftModel(s,spacetime);
			
			rr = new VectorN(tr);
			vv = new VectorN(tv);
			s = new Spacecraft(rr,vv,cr,cd,area,mass);
			s.set_use_params_in_state(false);
			truth[i] = new SpacecraftModel(s,spacetime);
			
			
			/*Set the step size for the trajectory generation*/
			/*Set the integrator Step size*/
			dt = initializer.parseInt(this.input,"init.dt");
			//ref[i].set_sc_dt(dt);
			//truth[i].set_sc_dt(dt);
			
			/* Read in the value for the spacecraft receiver noise */
			double range_noise=0;
			double state_noise[] = new double[3];
			str = "MEAS.types";
			String str2;
			int numMeas = initializer.parseInt(this.input,str);
			for(int m=0; m<numMeas; m++){
				str = "MEAS."+m+".desc";
				str2 = initializer.parseString(this.input,str);
				if(str2.equalsIgnoreCase("GPS")){
					str = "MEAS."+m+".R";
					range_noise = initializer.parseDouble(input,str);					
				}else if(str2.equalsIgnoreCase("pseudoGPS")){
					int end = initializer.parseInt(input,"MEAS."+m+".size");
					state_noise = new double[end];
					for(int snum=0; snum<end; snum++){
						str = "MEAS."+m+".R."+snum;
						state_noise[snum]=initializer.parseDouble(input,str);
					}						
				}
			}
			ref[i].set_GPS_noise(state_noise,range_noise);
		}		
	}
	public UniverseModel createUniverseModel(double mjd_utc,Spacecraft sc, boolean[] force_flag, boolean use_JGM2, String drag_model){
		
        boolean use_LP165P = false;
		UniverseModel umodel = new UniverseModel(mjd_utc);
		//ForceModelList forces = new ForceModelList();
		VectorN zero = new VectorN(0,0,0);
		if(force_flag[0]){
			System.out.println("Earth");
			GravitationalBody earth =
				new GravitationalBody(398600.4415e+9);
			umodel.addForce(earth);
		} else {
			if(use_JGM2){
				System.out.println("JGM2");
				GravityModel earth_grav = new GravityModel(2,2,EarthGravityType.JGM2);
				umodel.addForce(earth_grav);
			}else{
				System.out.println("JGM3");
				GravityModel earth_grav = new GravityModel(20,20,EarthGravityType.JGM3);
				umodel.addForce(earth_grav);
			}
			
		}
		if(force_flag[1]){
			System.out.println("Sun");
			umodel.set_compute_sun(true);
			Sun sun =
				new Sun(Constants.GM_Sun,zero,zero);
			umodel.addForce(sun);
		}
		if(force_flag[2]){
            umodel.set_compute_moon(true);
            if (use_LP165P) 
            {
                System.out.println("Moon - LP165P");
                GravityModel moon_grav = 
                  new GravityModel(2, 2, new LunaFixedRef(), 
                      "jat/forces/moonGravity/LP165P.grv");
                umodel.addForce(moon_grav);              
            }
            else 
            {
              System.out.println("Moon");
              Moon moon =
                new Moon(Constants.GM_Moon,zero,zero);
              umodel.addForce(moon);
            }
		}
		if(force_flag[3]){
			double ap_opt = 14.918648166;
			double f107_opt = 150;
			double n_param_opt = 6;
			umodel.set_compute_sun(true);
			if(drag_model.endsWith("NRL") || drag_model.endsWith("A") || drag_model.endsWith("C")){
				System.out.println("NRLMSISE");
				NRLMSISE_Drag drag = new NRLMSISE_Drag(sc);
				drag.setAP(ap_opt);
				drag.setF107Daily(f107_opt);
				drag.setF107Average(f107_opt);
				umodel.addForce(drag);
			}else{
				umodel.set_compute_sun(true);
				System.out.println("HarrisPriester");
				HarrisPriester atmos = new HarrisPriester(sc,150);//145.8480085177176);
				//atmos.setF107(145.8480085177176);//148.715);//99.5);
				atmos.setParameter(n_param_opt);
//				if(drag_model.equalsIgnoreCase("Sun-Sync"))
//				atmos.setParameter(6);
//				else if(drag_model.equalsIgnoreCase("ISS"))
//				atmos.setParameter(4);
				umodel.addForce(atmos);
				
			}
		}
		if(force_flag[4]){
			umodel.set_compute_sun(true);
			System.out.println("SolarRadiationPressure");
			SolarRadiationPressure srp = new SolarRadiationPressure(sc);
			umodel.addForce(srp);
		}
		return umodel;
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
		if(verbose_timestep){
			System.out.println("step: "+t+" / "+tf+"    stepsize: "+dt);
		}
		
		//rk8.setStepSize(sm.get_sc_dt());
		rk8.setStepSize(dt);
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
		iteration++;
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
		
		double tmpState[] = new double[6];
		for(int numSats = 0; numSats < numSpacecraft; numSats ++)
		{
			//Extract the state of the current satellite
			for(int i = 0;i < 6; i++)
			{
				tmpState[i]=newState.x[numSats*6 + i];
			}
			ref[numSats].get_spacecraft().updateMotion(tmpState);			
//			Write out the current True States
			double [] true_state = truth[numSats].get_spacecraft().toStateVector();
			
			//Print out simStep + 1 because we didn't output the initial state
			VectorN vecTime =  new VectorN(1,(simStep)*dt);
			VectorN trueState = new VectorN(true_state);
			VectorN truthOut = new VectorN(vecTime,trueState);
			new PrintStream(truths[numSats]).println (truthOut.toString());
			
			//Write out the current State estimates
			double [] ref_state  = ref[numSats].get_spacecraft().toStateVector();
			VectorN vecState = new VectorN(ref_state);
			VectorN stateOut = new VectorN(vecTime,vecState);
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
			int numStates = filter.get_numStates();
			double[] tmp = new double[numStates*numStates];
			int k = 0;
			for(int i = 0; i < numStates; i++)
			{
				for(int j = 0; j < numStates; j++)
				{
					tmp[k] = Covariance.get(i,j);
					k++;
				}
			}
			VectorN ErrCov = new VectorN(tmp);
			stateOut = new VectorN(vecTime,ErrCov);
			new PrintStream(covariance[numSats]).println (stateOut.toString());
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
		
		for(int i=0; i<numSpacecraft; i++){
			truth_traj[i].add(truth[i].get_sc_t(),truth[i].get_spacecraft().toStateVector());
			ref_traj[i].add(ref[i].get_sc_t(),ref[i].get_spacecraft().toStateVector());
		}
		
		/*Cache off the simulation mode */
		int filterMode = initializer.parseInt(this.input,"init.mode");
		
		//Compute the length of the simulation in seconds
		double MJD0 =  initializer.parseDouble(this.input,"init.MJD0");		
		double MJDF =  initializer.parseDouble(this.input,"init.MJDF");
		double T0   =  initializer.parseDouble(this.input,"init.T0");
		double TF   =  initializer.parseDouble(this.input,"init.TF");
		//simTime = 0; //* this is done in call to "initialize()"
		this.mjd_utc_start = MJD0+T0/86400.0;
		simTime = new Time(MJD0+T0/86400.0);
		double simLength = Math.round((MJDF - MJD0)*86400 + TF - T0);
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
//				int donothing = 0;
//				donothing++;
//			}
			//System.out.println("SimTime: " + simTime.get_sim_time() + " SimStep: " + simStep);
			
			for(int i=0; i<numSpacecraft; i++){
				truth_traj[i].add(truth[i].get_sc_t(),truth[i].get_spacecraft().toStateVector());
				ref_traj[i].add(ref[i].get_sc_t(),ref[i].get_spacecraft().toStateVector());
			}
			
		}
		
		/*Close all output files*/
		closeFiles();
		System.gc();
		
		double elapsed = (System.currentTimeMillis()-start)*0.001/60;
		System.out.println("Elapsed time [min]: "+elapsed);
		
		/* Post Processing */
		
		LinePrinter lp = new LinePrinter();
		RelativeTraj[] reltraj = new RelativeTraj[numSpacecraft];
		double mismatch_tol = 0.00001;
		//* TODO Plot marker
		for(int i=0; i<numSpacecraft; i++){
			if(PlotJAT){
				reltraj[i] = new RelativeTraj(ref_traj[i],truth_traj[i],lp,"Jat(Ref) v Jat(Truth)");
				reltraj[i].setVerbose(false);
				reltraj[i].process(mismatch_tol);
				reltraj[i].process_RSS(mismatch_tol);
				//reltraj[i].process_ECI(mismatch_tol);
			}

			try {
				Celestia cel = new Celestia("C:/Code/Celestia/");
				cel.set_trajectory_meters(ref_traj[i],MJD0);
				cel.write_trajectory("jat_ref_"+JAT_name+JAT_case,"jat_ref_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(this.mjd_utc_start));
				cel.set_trajectory_meters(truth_traj[i],MJD0);
				cel.write_trajectory("jat_truth_"+JAT_name+JAT_case,"jat_truth_"+JAT_name+JAT_case,TimeUtils.MJDtoJD(this.mjd_utc_start));
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
	
//	** Main **//
	
	public static void main(String[] args) {
		
		boolean useFilter = true;
		
		//* TODO Flag marker
		EstimatorSimModel.JAT_case = 44;
		
		CEVSim.JAT_name = "moon2earth_";
		//CEVSim.InputFile = "initialConditions_cev_m2e.txt";
		CEVSim.InputFile = "initialConditions_cev_m2e_HIGH.txt";
		
//		CEVSim.JAT_name = "earth2moon";
//		CEVSim.InputFile = (args.length > 0 ? args[0] : "initialConditions_cev_e2m.txt");
		
//		CEVSim.JAT_name = "earth2moon";
//		CEVSim.InputFile = "initialConditions_cev_e2m_bias.txt";
		
		CEVSim.PlotJAT = true;
		
		CEVSim Sim = new CEVSim(useFilter);
		Sim.set_verbose(true);
		Sim.runloop();
        if (OpticalMeasurementModel.fobs != null) {
          OpticalMeasurementModel.fobs.close();
        }
        if (OpticalMeasurementModel.fpred != null) {
          OpticalMeasurementModel.fpred.close();
        }
        System.out.println("Finished.");
	}
}
