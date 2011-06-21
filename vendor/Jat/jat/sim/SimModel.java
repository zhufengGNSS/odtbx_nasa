/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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
 */
package jat.sim;

import jat.spacecraft.*;
import jat.spacetime.UniverseModel;
import jat.traj.*;
import jat.math.MathUtils;
//import jat.matlabInterface.MatlabControl;
import jat.matlabInterface.MatlabFunc;
import jat.matvec.data.VectorN;
import jat.alg.integrators.*;
import jat.cm.Constants;
//import jat.forces.ForceModelList;
import jat.forces.ForceModel;
import jat.forces.GravitationalBody;
import jat.forces.gravity.*;
import jat.forces.gravity.earth.*;
import jat.forces.density.earth.*;
import jat.forces.Moon;
import jat.forces.SolarRadiationPressure;
import jat.forces.Sun;
import jat.util.*;

/**
 * This is the primary helper class for Simulation.java.  It encapsulates the
 * various spacecraft to be propagated along with the reference frame and time
 * information needed in the simulation.
 *
 * @author Richard C. Page III
 *
 */
public class SimModel implements Derivatives {

	/**
	 * Flag indicating whether to print to file
	 */
	protected boolean doPrint = true;
    /**
     * Eigth Order Runge Kutta integration algorithm
     */
    protected RungeKutta8 rk8;
    /**
     * Spacetime model containing: time, reference frames, and forces
     */
    public UniverseModel spacetime;
    /**
     * Spacecraft model used when propagating a single spacecraft
     */
    public SpacecraftModel sc;
    /**
     * Flag indicating whether multiple spacecraft are present
     */
    protected boolean use_formation;
    /**
     * A container class carrying the various spacecraft and their relationships
     */
    protected Formation sc_formation;
    /**
     * Printer to write data to a file or the command line
     */
    protected LinePrinter lp;
    /**
     * Counter to keep track of the current spacecraft being propagated
     */
    protected int current_sc=0;
    /**
     * Initial time offset of the simulation [sec] (default is zero seconds)
     */
    public double t0 = 0;
    /**
     * Final time of the simulation [sec] (default is 86400 seconds)
     */
    public double tf = 86400;
    /**
     * Simulation time [sec]
     */
    public double t = t0;
    /**
     * The start date of the simulation in Modified Julian Date Universal Coordinated Time (UTC)
     */
    public double mjd_utc_start;
    /**
     * Integrator stepsize
     */
    public double stepsize=60;
    /**
     * Counter for the number of iterations
     */
    public int iteration = 0;
    /**
     * Reference trajectories for doing analysis on single spacecraft
     */
    private Trajectory jat,truth;
    /**
     * Reference trajectories for doing analysis on spacecraft formations
     */
    private Trajectory[] traj_formation;
    /**
     * Flag to tell whether to print the progress of the simulation.
     */
    protected boolean verbose_timestep=false;
    /**
     * Properties member which loads and stores ephemeris data.
     */
    protected SimProperties properties;
    /**
     * Default ephemeris filename.
     */
    protected String default_file;
    
    /**
     * Default Constructor initializes the universe model, integrator, and prints
     * to the command line.
     */
    public SimModel(){
        rk8 = new RungeKutta8();
        lp = new LinePrinter();
        spacetime = new UniverseModel();
    }
    /**
     * Constructor initializes a single spacecraft simulation given relevant parameters
     * @param r Position vector [m]
     * @param v Velocity vector [m/s]
     * @param cr Coefficient of Reflectivity used for Solar Radiation Pressure
     * @param cd Coefficient of Drag used for Atmospheric Drag calculations
     * @param area Cross sectional area used both for drag and Solar Radiation Pressure
     * @param mass Mass of the spacecraft
     */
    public SimModel(double[] r, double[] v, double cr, double cd, double area, double mass){
        VectorN rr = new VectorN(r);
        VectorN vv = new VectorN(v);
        Spacecraft s = new Spacecraft(rr,vv,cr,cd,area,mass);
        s.set_use_params_in_state(false);
        sc = new SpacecraftModel(s);
        rk8 = new RungeKutta8();
        use_formation = false;
        lp = new LinePrinter();
        spacetime = new UniverseModel();
    }

    
    /**
     * Constructor initializes a single spacecraft simulation given relevant parameters
     * @param r Position vector [m]
     * @param v Velocity vector [m/s]
     * @param cr Coefficient of Reflectivity used for Solar Radiation Pressure
     * @param cd Coefficient of Drag used for Atmospheric Drag calculations
     * @param area Cross sectional area used both for drag and Solar Radiation Pressure
     * @param mass Mass of the spacecraft
     */
    public SimModel(double[] r, double[] v, double cr, double cd, double area, double mass, double utc){
        VectorN rr = new VectorN(r);
        VectorN vv = new VectorN(v);
        Spacecraft s = new Spacecraft(rr,vv,cr,cd,area,mass);
        s.set_use_params_in_state(false);
        sc = new SpacecraftModel(s);
        rk8 = new RungeKutta8();
        use_formation = false;
        lp = new LinePrinter();
        spacetime = new UniverseModel(utc);
    }

    
    
    /**
     * Constructor initializes a formation of spacecraft given relevant parameters
     * @param r Array of position vectors r[number_of_sc][3] [m]
     * @param v Array of velocity vectors v[number_of_sc][3] [m/s]
     * @param cr Array of Coefficients of Reflectivity used for Solar Radiation Pressure
     * @param cd Array of Coefficients of Drag used for Atmospheric Drag calculations
     * @param area Array of Cross-sectional area for both Drag and Solar Radiation
     * @param mass Array of spacecraft masses
     */
    public SimModel(double[][] r, double[][] v, double[] cr, double[] cd,
            double[] area, double[] mass){
        rk8 = new RungeKutta8();
        VectorN rr,vv;
        Spacecraft s;
        sc_formation = new Formation();
        for(int i=0; i<mass.length; i++){
            rr = new VectorN(r[i]);
            vv = new VectorN(v[i]);
            s = new Spacecraft(rr,vv,cr[i],cd[i],area[i],mass[i]);
            s.set_use_params_in_state(false);
            sc_formation.add_spacecraft(s);
        }
        use_formation = true;
        lp = new LinePrinter();
        spacetime = new UniverseModel();
    }
    


    /**
     * This is the generic initialization method for the spacecraft and simulation
     * properties.  It allows the user to create their own input parameters or
     * parsers and feed them into this method.
     *
     * @param sm An array of spacecraft models.  The first is considered the primary
     * spacecraft of the formation.
     * @param start The start offset for the simulation time [sec]
     * @param finish The final time of the simulation [sec]
     * @param mjd_utc Simulation epoch in Modified Julian Date of Universal Coordinated Time
     * @param step Integrator timestep [sec]
     * @param thin Thinning parameter.  Indicates how many timesteps between printing data.
     * (thin = 1 means print all timesteps, thin = 2 means print every other timestep)
     * @param output Filename for the tab-delimited output data
     */
    public void initialize(SpacecraftModel[] sm, double start, double finish,
            double mjd_utc, double step, int thin, String output){
        this.use_formation = true;
        traj_formation = new Trajectory[sm.length];
        sc_formation = new Formation();
        for(int i=0; i<sm.length; i++){
            sc_formation.add_spacecraft(sm[i]);
            traj_formation[i] = new Trajectory();
        }
        this.t0 = start;
        this.tf = finish;
        this.mjd_utc_start = mjd_utc;
        this.stepsize = step;
        lp = new LinePrinter(output);
        lp.setThinning(thin);
        //* Load Simulation
        spacetime.set_time(mjd_utc_start);
        jat = new Trajectory();
        truth = new Trajectory();
    }

    /**
     * This is the generic initialization method for the spacecraft and simulation
     * properties.  It allows the user to create their own input parameters or
     * parsers and feed them into this method.
     *
     * @param sm An array of spacecraft models.  The first is considered the primary
     * spacecraft of the formation.
     * @param start The start offset for the simulation time [sec]
     * @param finish The final time of the simulation [sec]
     * @param mjd_utc Simulation epoch in Modified Julian Date of Universal Coordinated Time
     * @param step Integrator timestep [sec]
     * @param thin Thinning parameter.  Indicates how many timesteps between printing data.
     * (thin = 1 means print all timesteps, thin = 2 means print every other timestep)
     * @param output Filename for the tab-delimited output data
     */
    public void initialize(SpacecraftModel sm, double start, double finish,
            double mjd_utc, double step, int thin, String output){
        this.use_formation = false;
        sc = sm;
        this.t0 = start;
        this.tf = finish;
        this.mjd_utc_start = mjd_utc;
        this.stepsize = step;
        lp = new LinePrinter(output);
        lp.setThinning(thin);
        //* Load Simulation
        spacetime.set_time(mjd_utc_start);
        jat = new Trajectory();
        truth = new Trajectory();
    }

    /**
     * This is the initialize method without the print method or line thining.
     *
     * @param sm An array of spacecraft models.  The first is considered the primary
     * spacecraft of the formation.
     * @param start The start offset for the simulation time [sec]
     * @param finish The final time of the simulation [sec]
     * @param mjd_utc Simulation epoch in Modified Julian Date of Universal Coordinated Time
     * @param step Integrator timestep [sec]
     * This initialize method does not use print method or line thining.
     */
    public void initialize(SpacecraftModel sm, double start, double finish,
            double mjd_utc, double step){
        this.use_formation = false;
        sc = sm;
        this.t0 = start;
        this.tf = finish;
        this.mjd_utc_start = mjd_utc;
        this.stepsize = step;
        //* Load Simulation
        spacetime.set_time(mjd_utc_start);
        jat = new Trajectory();
        truth = new Trajectory();
    }

    /**
     * Initialize the Simulation from a Matlab session.  Note this method does not work
     * when calling it from a standalone Java session.
     * @param init_sc Name of the Matlab function which returns the spacecraft state and parameters
     * @param init_integ Name of the Matlab function which returns the integration parameters
     * @param output Name of the output file of the simulation with extension
     * @throws InterruptedException
     */
    public void initializeMatlab(String init_sc, String init_integ, String output) throws InterruptedException{
//        MatlabControl input = new MatlabControl();
//        input.eval("disp('print')");
        MatlabFunc input_sc = new MatlabFunc(init_sc);
        MatlabFunc input_integ = new MatlabFunc(init_integ);
        Object[] args = new Object[1];
        args[0] = new Object();
        double[] get = (double[])input_sc.call(args);
        double[] param = (double[])input_integ.call(args);
        System.out.println("getlength: "+get.length);
        int numsc = get.length/10;
        Spacecraft s;
        double[] X = new double[10];
        if(numsc>1){
            sc_formation = new Formation();
            this.use_formation = true;
            traj_formation = new Trajectory[numsc];
            for(int i=0; i<numsc; i++){
                traj_formation[i] = new Trajectory();
                for(int j=0; j<10; j++){
                    X[j] = get[i*10+j];
                    //System.out.println("get["+(i*10+j)+"]: "+get[i*10+j]);
                }
                s = new Spacecraft(X);
                sc_formation.add_spacecraft(s);
            }
        } else {
            s = new Spacecraft(get);
            sc = new SpacecraftModel(s);
        }
        //System.out.println("param: "+param[0]+" "+param[1]+" "+param[2]+" "+param[3]);
        mjd_utc_start = param[0];
        t0 = param[1];
        tf = param[2];
        stepsize = param[3];
        lp = new LinePrinter(output);
        lp.setThinning((int) param[4]);
        //* Load Simulation
        spacetime.set_time(mjd_utc_start);
    }

    /**
     * Initialize the Simulation from a Matlab session.  Note this method does not work
     * when calling it from a standalone Java session.
     * @param param Integration parameters [mjd_utc_start t0 tf stepsize thinning]
     * where t0 and tf are in seconds and thinning is the number of steps between
     * output (ie thinning = 1 means print every step, thinning = 2 every other step)
     * @param output Name of the output file of the simulation with extension
     * @throws InterruptedException
     */
    public void initializeMatlab(double[] param, String output){
        this.use_formation = true;
        	int numsc = sc_formation.get_num_sc();
            traj_formation = new Trajectory[numsc];
            for(int i=0; i<numsc; i++){
                traj_formation[i] = new Trajectory();
            }
        mjd_utc_start = param[0];
        t0 = param[1];
        tf = param[2];
        stepsize = param[3];
        lp = new LinePrinter(output);
        lp.setThinning((int) param[4]);
        //* Load Simulation
        spacetime.set_time(mjd_utc_start);
    }

    /**
     * Initialize the forces present in the universe model during the Simulation.
     * @param force_flag An array of boolean values indicating the forces in order:
     *  [0] true: two-body gravity false: nonspherical gravity
     *  [1] true: Solar gravity [2] true: Lunar gravity [3] true: Atmospheric drag
     *  [4] true: Solar Radiation Pressure
     * @param use_JGM2 If using nonspherical gravity true selects JGM2 instead of JGM3
     * @param drag_model "NRL" for NRLMSISE2000 or "HP" for Harris Priester
     */
    public void initializeForces(boolean[] force_flag, boolean use_JGM2, String drag_model){

        //ForceModelList forces = new ForceModelList();
        VectorN zero = new VectorN(0,0,0);
	    if(force_flag[0]){
	        GravitationalBody earth =
	            new GravitationalBody(398600.4415e+9);
	        spacetime.addForce(earth);
	    } else {
	        if(use_JGM2){
	            GravityModel earth_grav = new GravityModel(2,2,EarthGravityType.JGM2);
	            spacetime.addForce(earth_grav);
	            spacetime.set_use_iers(true);
	        }else{
	            GravityModel earth_grav = new GravityModel(20,20,EarthGravityType.JGM3);
	            spacetime.addForce(earth_grav);
	            spacetime.set_use_iers(true);
	        }

	    }
	    if(force_flag[1]){
	        spacetime.set_compute_sun(true);
	        Sun sun =
	            new Sun(Constants.GM_Sun,zero,zero);
	        spacetime.addForce(sun);
	    }
	    if(force_flag[2]){
	        spacetime.set_compute_moon(true);
	        Moon moon =
	            new Moon(Constants.GM_Moon,zero,zero);
	        spacetime.addForce(moon);
	    }
	    if(force_flag[3]){
	        double ap_opt = 14.918648166;
            double f107_opt = 150;
            double n_param_opt = 6;
            spacetime.set_compute_sun(true);
	        if(drag_model.endsWith("NRL") || drag_model.endsWith("A") || drag_model.endsWith("C")){
	            NRLMSISE_Drag drag = new NRLMSISE_Drag(sc.get_spacecraft());
				drag.setAP(ap_opt);
				drag.setF107Daily(f107_opt);
				drag.setF107Average(f107_opt);
	            spacetime.addForce(drag);
	            spacetime.set_use_iers(true);
	        }else{
	            spacetime.set_compute_sun(true);
	            HarrisPriester atmos = new HarrisPriester(sc.get_spacecraft(),150);//145.8480085177176);
	            //atmos.setF107(145.8480085177176);//148.715);//99.5);
	            atmos.setParameter(n_param_opt);
//	            if(drag_model.equalsIgnoreCase("Sun-Sync"))
//	                atmos.setParameter(6);
//	            else if(drag_model.equalsIgnoreCase("ISS"))
//	                atmos.setParameter(4);
	            spacetime.addForce(atmos);
	            spacetime.set_use_iers(true);

	        }
	    }
	    if(force_flag[4]){
	        spacetime.set_compute_sun(true);
	        SolarRadiationPressure srp = new SolarRadiationPressure(sc.get_spacecraft());
	        spacetime.addForce(srp);
	    }
    }

    public void initializeForcesMatlab(int[] flags,String drag){
        boolean[] case_flags = new boolean[flags.length-1];
        for(int i=0; i<flags.length-1; i++){
            if(flags[i]==0) case_flags[i] = false;
            else case_flags[i] = true;
        }
        boolean use_JGM2;
        if(flags[flags.length-1]==0) use_JGM2 = false;
        else use_JGM2 = true;
        initializeForces(case_flags, use_JGM2, drag);
    }

    /**
     * Initialize the forces present in the universe model during the Simulation.
     * @param force_flag An array of boolean values indicating the forces in order:
     *  [0] true: two-body gravity false: nonspherical gravity
     *  [1] true: Solar gravity [2] true: Lunar gravity [3] true: Atmospheric drag
     *  [4] true: Solar Radiation Pressure
     * @param use_JGM2 If using nonspherical gravity true selects JGM2 instead of JGM3
     * @param drag_model "NRL" for NRLMSISE2000 or "HP" for Harris Priester
     * @param order DMS 5/18/07
     * @param degree DMS 5/18/07
     */
    public void initializeForces(boolean[] force_flag, boolean use_JGM2, String drag_model, int order, int degree){

        //ForceModelList forces = new ForceModelList();
        VectorN zero = new VectorN(0,0,0);
	    if(force_flag[0]){
	        GravitationalBody earth =
	            new GravitationalBody(398600.4415e+9);
	        spacetime.addForce(earth);
	    } else {
	        if(use_JGM2){
	            GravityModel earth_grav = new GravityModel(degree,order,EarthGravityType.JGM2);
	            spacetime.addForce(earth_grav);
	            spacetime.set_use_iers(true);
	        }else{
	            GravityModel earth_grav = new GravityModel(degree,order,EarthGravityType.JGM3);
	            spacetime.addForce(earth_grav);
	            spacetime.set_use_iers(true);
	        }

	    }
	    if(force_flag[1]){
	        spacetime.set_compute_sun(true);
	        Sun sun =
	            new Sun(Constants.GM_Sun,zero,zero);
	        spacetime.addForce(sun);
	    }
	    if(force_flag[2]){
	        spacetime.set_compute_moon(true);
	        Moon moon =
	            new Moon(Constants.GM_Moon,zero,zero);
	        spacetime.addForce(moon);
	    }
	    if(force_flag[3]){
	        double ap_opt = 14.918648166;
            double f107_opt = 150;
            double n_param_opt = 6;
            spacetime.set_compute_sun(true);
	        if(drag_model.endsWith("NRL") || drag_model.endsWith("A") || drag_model.endsWith("C")){
	            NRLMSISE_Drag drag = new NRLMSISE_Drag(sc.get_spacecraft());
				drag.setAP(ap_opt);
				drag.setF107Daily(f107_opt);
				drag.setF107Average(f107_opt);
	            spacetime.addForce(drag);
	            spacetime.set_use_iers(true);
	        }else{
	            spacetime.set_compute_sun(true);
	            HarrisPriester atmos = new HarrisPriester(sc.get_spacecraft(),150);//145.8480085177176);
	            //atmos.setF107(145.8480085177176);//148.715);//99.5);
	            atmos.setParameter(n_param_opt);
//	            if(drag_model.equalsIgnoreCase("Sun-Sync"))
//	                atmos.setParameter(6);
//	            else if(drag_model.equalsIgnoreCase("ISS"))
//	                atmos.setParameter(4);
	            spacetime.addForce(atmos);
	            spacetime.set_use_iers(true);

	        }
	    }
	    if(force_flag[4]){
	        spacetime.set_compute_sun(true);
	        SolarRadiationPressure srp = new SolarRadiationPressure(sc.get_spacecraft());
	        spacetime.addForce(srp);
	    }
    }

    public void initializeForcesMatlab(int[] flags,String drag, int degree, int order){
        boolean[] case_flags = new boolean[flags.length-1];
        for(int i=0; i<flags.length-1; i++){
            if(flags[i]==0) case_flags[i] = false;
            else case_flags[i] = true;
        }
        boolean use_JGM2;
        if(flags[flags.length-1]==0) use_JGM2 = false;
        else use_JGM2 = true;
        initializeForces(case_flags, use_JGM2, drag, order, degree);
    }

    /**
     * Add a force model to the simulation.
     * @param f Force Model
     */
    public void add_ForceModel(ForceModel f){
        spacetime.addForce(f);
    }

    /**
     * Add a spacecraft model to the simulation.
     * @param sm Spacecraft Model
     */
    public void add_Spacecraft(SpacecraftModel sm){
        this.sc_formation.add_spacecraft(sm);
    }

    /**
     * Add a spacecraft from a state vector.
     * @param X State Vector.
     */
    public void add_Spacecraft(double[] X){
        if(sc_formation == null) sc_formation = new Formation();
        this.sc_formation.add_spacecraft(new Spacecraft(X));
    }

    /**
     * Add a spacecraft to the simulation.
     * @param s Spacecraft.
     */
    public void add_Spacecraft(Spacecraft s){
        this.sc_formation.add_spacecraft(s);
    }
    /**
     * Get the spacecraft model with the given index.
     * @param i Index
     * @return Spacecraft Model
     */
    public SpacecraftModel get_SpacecraftModel(int i){
        return sc_formation.get_spacecraftmodel(i);
    }

    /**
     * Get the spacecraft with the given index.
     * @param i Index
     * @return Spacecraft
     */
    public Spacecraft get_Spacecraft(int i){
        return sc_formation.get_spacecraft(i);
    }
    /**
     * Get the Spacecraft model with the given ID or return null
     * @param id ID
     * @return Spacecraft Model
     */
    public SpacecraftModel get_SpacecraftModel(String id){
        return sc_formation.get_spacecraftmodel(id);
    }
    /**
     * Get the Spacecraft with the given ID or return null
     * @param id ID
     * @return Spacecraft Model
     */
    public Spacecraft get_Spacecraft(String id){
        return sc_formation.get_spacecraft(id);
    }

    /**
     * Get the number of spacecraft in the current simulation
     * @return The number of spacecraft
     */
    public int get_number_of_spacecraft(){
    	return sc_formation.get_num_sc();
    }

    /**
     * Set the integrator stepsize
     * @param s step size [sec]
     */
    public void set_stepsize(double s){
        rk8.setStepSize(s);
    }

    /**
     * Choose a controller for a particular spacecraft, numbered as they were added
     * to the current simulation
     * @param sc_num The spacecraft's number
     * @param c The control law to be stored to spacecraft 'sc_num'
     */
    public void set_controller(int sc_num, ControlLaw c){
        if(!use_formation){
            sc.set_controller(c);
        } else {
            SpacecraftModel s = sc_formation.get_spacecraftmodel(sc_num);
            s.set_controller(c);
        }
    }

    /**
     * If propagating only one spacecraft apply the control law.  If propagating
     * multiple spacecraft, apply the control law to the primary spacecraft or
     * to the first spacecraft added.
     * @param c The control law to be applied.
     */
    public void set_controller(ControlLaw c){
        if(!use_formation){
            sc.set_controller(c);
        } else {
            SpacecraftModel s = sc_formation.get_primarymodel();
            s.set_controller(c);
        }
    }

    /**
     * For use with Matlab - sets a controller from an m-file.
     * @param num Spacecraft number to apply controller
     * @param cmd m-file name
     */
    public void set_Matlab_controller(int num, String cmd){
        MatlabControlLaw mcontroller = new MatlabControlLaw(cmd);
        sc_formation.get_spacecraftmodel(num).set_controller(mcontroller);
    }

    /**
     * Choose an estimator for a particular spacecraft, numbered as they were added
     * to the current simulation
     * @param sc_num The spacecraft's number
     * @param e Estimator to be set
     */
    public void set_estimator(int sc_num, StateEstimation e){
        if(!use_formation){
            sc.set_estimator(e);
        } else {
            SpacecraftModel s = sc_formation.get_spacecraftmodel(sc_num);
            s.set_estimator(e);
        }
    }

    /**
     * If propagating only one spacecraft apply the estimator.  If propagating
     * multiple spacecraft, apply the estimator to the primary spacecraft or
     * to the first spacecraft added.
     * @param e Estimator to be applied
     */
    public void set_estimator(StateEstimation e){
        if(!use_formation){
            sc.set_estimator(e);
        } else {
            SpacecraftModel s = sc_formation.get_primarymodel();
            s.set_estimator(e);
        }
    }

    /**
     * Set the 'truth' trajectory to analyze the simulation output error
     * @param filename Tab delimited file containing [MJD_UTC x y z xdot ydot zdot]
     */
    public void set_truth_traj(String filename){ truth.readFromFile(filename);}

    /**
     * Set whether to show the time progression during the simulation loop
     * @param b = true to print the time
     */
    public void set_showtimestep(boolean b){
    	this.verbose_timestep = b;
    }
    
    /**
     * Sets whether to print to file or not.
     * @param b = true to print to file
     */
    public void set_doPrint(boolean b){
    	this.doPrint = b;
    }

    /**
     * Get the Trajectory obtained from propagating a single spacecraft
     * @return The trajectory: [MJD_UTC x y z xdot ydot zdot]
     */
    public Trajectory get_traj(){ return jat;}

    /**
     * Get the Trajectory of the requestes spacecraft
     * @param sc_num The spacecraft's number in the sequence as added to the simulation
     * @return The trajectory: [MJD_UTC x y z xdot ydot zdot]
     */
    public Trajectory get_traj(int sc_num){ return traj_formation[sc_num];}

    /**
     * Get the Relative Trajectory comparison between the 'truth' and the simulation.
     * The 'truth' trajectory must be set.
     * @see jat.sim.SimModel#set_truth_traj(String)
     * @param printer Allows writing to a file or the command line
     * @return The relative trajectory comparison
     */
    public RelativeTraj get_rel_traj(LinePrinter printer){
        if(this.use_formation)
            return new RelativeTraj(traj_formation[0],truth,printer);
        else
            return new RelativeTraj(jat,truth,printer);
    }

    /**
     * Get the Relative Trajectory information of the requested spacecraft in the
     * formation.  This method requires more than one spacecraft.
     * @param l Printer used to write to file or command line
     * @param sc The number of the spacecraft to compare against the primary (spacecraft zero)
     * @return The relative trajectory
     */
    public RelativeTraj get_rel_traj(LinePrinter l, int sc){
        return new RelativeTraj(traj_formation[0],traj_formation[sc],l);
    }

    /**
     * Get the Relative Trajectory information of the requested spacecraft in the
     * formation.  This method requires more than one spacecraft.
     * @param l Printer used to write to file or command line
     * @param scOne The number of the first spacecraft to compare
     * @param scTwo The number of the second spacecraft to compare
     * @return The relative trajectory
     */
    public RelativeTraj get_rel_traj(LinePrinter l, int scOne, int scTwo){
        return new RelativeTraj(traj_formation[scOne],traj_formation[scTwo],l);
    }

    public void plot_rel_traj(int scOne, int scTwo){
        LinePrinter p = new LinePrinter();
        RelativeTraj rel = new RelativeTraj(traj_formation[scOne],traj_formation[scTwo],p);
        rel.setVerbose(false);
        rel.process();
    }

    /**
     * Primary propagation method.  Increments the simulation by a time 'dt'.
     * Updates the universe model according to the progression of time.
     * Updates the spacecraft state.
     * @param dt Timestep in seconds
     */
    public void step(double dt){

    	if(verbose_timestep){
    		System.out.println("step: "+t+" / "+tf+"    stepsize: "+dt);
    	}

        rk8.setStepSize(dt);
        //* update models
        double mjd_utc = spacetime.get_mjd_utc();
        double[] X = new double[6];
        double[] Xnew = new double[6];
        double[] thrust = new double[3];
        VectorN rnew;
        VectorN vnew;
        double[] tmp = new double[6];
        double num_sc = 1;
        if(use_formation){
            num_sc = sc_formation.get_num_sc();
            sc_formation.compute_control(t);
            for(int i=0; i<num_sc; i++){
                current_sc = i;
                //* call integrator
                SpacecraftModel sm = sc_formation.get_spacecraftmodel(current_sc);
                Spacecraft s = sm.get_spacecraft();
                X = s.toStateVector(false);
                Xnew = rk8.step(t, X, this);
                //* store new values
                rnew = new VectorN(Xnew[0],Xnew[1],Xnew[2]);
                vnew = new VectorN(Xnew[3],Xnew[4],Xnew[5]);
                s.updateMotion(rnew,vnew);

            }
        } else {
            //* call integrator
            Spacecraft s = sc.get_spacecraft();
            X = s.toStateVector(false);
            Xnew = rk8.step(t, X, this);
            //* store new values
            rnew = new VectorN(Xnew[0],Xnew[1],Xnew[2]);
            vnew = new VectorN(Xnew[3],Xnew[4],Xnew[5]);
            s.updateMotion(rnew,vnew);
        }
        //* update simulation time
        if(t > (tf-dt) && t != tf){
            dt = tf-t;
        }
        t = t+dt;
        //* update the universe
        spacetime.update(t);
        iteration++;
    }

    /**
     * Primary propagation method.  Increments the simulation by a time 'dt'.
     * Updates the universe model according to the progression of time.
     * Updates the spacecraft state.
     * @param dt Timestep in seconds
     */
    public double[] stepMatlabAdaptor(double dt, int neqns){

    	if(verbose_timestep){
    		System.out.println("step: "+t+" / "+tf+"    stepsize: "+dt);
    	}

        rk8.setStepSize(dt);
        //* update models
        //double mjd_utc = spacetime.get_mjd_utc();
        double[] X = new double[neqns];
        double[] Xnew = new double[neqns];
        //double[] thrust = new double[3];
        VectorN rnew;
        VectorN vnew;
        //double[] tmp = new double[neqns];
        double num_sc = 1;
        if(use_formation){
            num_sc = sc_formation.get_num_sc();
            sc_formation.compute_control(t);
            for(int i=0; i<num_sc; i++){
                current_sc = i;
                //* call integrator
                SpacecraftModel sm = sc_formation.get_spacecraftmodel(current_sc);
                Spacecraft s = sm.get_spacecraft();
                X = s.toStateVector(false);
                Xnew = rk8.step(t, X, this);
                //* store new values
                rnew = new VectorN(Xnew[0],Xnew[1],Xnew[2]);
                vnew = new VectorN(Xnew[3],Xnew[4],Xnew[5]);
                s.updateMotion(rnew,vnew);

            }
        } else {
            //* call integrator
            Spacecraft s = sc.get_spacecraft();
            X = s.toStateVector(false);
            Xnew = rk8.step(t, X, this);
            //* store new values
            rnew = new VectorN(Xnew[0],Xnew[1],Xnew[2]);
            vnew = new VectorN(Xnew[3],Xnew[4],Xnew[5]);
            s.updateMotion(rnew,vnew);
        }
        //* update simulation time
        if(t > (tf-dt) && t != tf){
            dt = tf-t;
        }
        t = t+dt;
        //* update the universe
        spacetime.update(t);
        iteration++;
        return Xnew;
    }

    /**
     * Used to increment the simulation and output to a file or the command line.
     * @see jat.sim.SimModel#step(double)
     * @param dt Timestep in seconds
     * @param lp Printer to file or command line
     */
    public void step(double dt, LinePrinter lp){
        step(dt);
        if(MathUtils.mod(iteration,lp.getThinning())== 0)
            print(lp);
    }

    /**
     * Method used to propagate the simulation from initial to final states.
     * Calls SimModel.stepDefault(dt)
     */
    public double[][] runloopMatlabAdaptor(double[] x0, double[] time){

        double[][] newX0 = new double[time.length][x0.length];
        double[] timeArray = new double[time.length];
        double[][] output = new double[time.length][x0.length+1];
        int neqns = x0.length;
        int counter = 1;
        t = t0;
        timeArray[0]=time[0];
        newX0[0]= x0;
//      * update the universe
        spacetime.update(t);
        //* loop over time equal to the sim duration

        // check both the time and the counter index because floating-point
        // errors can cause the counter to accidentally overrun
        while((t < tf) && (counter < timeArray.length)){
            newX0[counter]=stepMatlabAdaptor(stepsize, neqns);
            timeArray[counter]=t;
            counter++;
        }

        for (int i=0; i<time.length; i++)
        {
        	for (int j=0; j<x0.length+1; j++)
        	{
        		if (j==0)
   			 {
   				 output[i][j]= timeArray[i];
   			 }
   			 else
   			 {
   				 output[i][j] = newX0[i][j-1];
   			 }
        	}
        }
        //System.out.println("Loop Finished");
        return output;
    }

    /**
     * Method used to propagate the simulation from initial to final states.
     * Calls SimModel.step(double,lp)
     */
    public void runloop(){
        //* print initial state
    	print(lp);
        t = t0;
//      * update the universe
        spacetime.update(t);
        //* loop over time equal to the sim duration
        while(t < tf){
            step(stepsize,lp);
        }
        System.out.println("Loop Finished");
        lp.close();
    }

    /**
     * Print the current state(s) to file or command line.
     * Store state(s) to trajectory containers.
     * @param lp Printer to file or command line
     */
    protected void print(LinePrinter lp){
        double[] X = new double[6];
        Spacecraft s;
        double mjd_utc;
        if(this.use_formation){
            int num_sc = sc_formation.get_num_sc();
            for(int i=0;i<num_sc;i++){
                s = sc_formation.get_spacecraft(i);
                mjd_utc = spacetime.get_mjd_utc();
                X[0] = s.r().get(0)/1000.0;
                X[1] = s.r().get(1)/1000.0;
                X[2] = s.r().get(2)/1000.0;
                X[3] = s.v().get(0)/1000.0;
                X[4] = s.v().get(1)/1000.0;
                X[5] = s.v().get(2)/1000.0;
                if(this.doPrint){
                	lp.println(mjd_utc+"\t"+X[0]+"\t"+X[1]+"\t"+X[2]+"\t"
                			+X[3]+"\t"+X[4]+"\t"+X[5]);
                }
                traj_formation[i].add(mjd_utc,X);
            }
        } else {
            s = sc.get_spacecraft();
            mjd_utc = spacetime.get_mjd_utc();
            X[0] = s.r().get(0)/1000.0;
            X[1] = s.r().get(1)/1000.0;
            X[2] = s.r().get(2)/1000.0;
            X[3] = s.v().get(0)/1000.0;
            X[4] = s.v().get(1)/1000.0;
            X[5] = s.v().get(2)/1000.0;
            if(this.doPrint){
            	lp.println(mjd_utc+"\t"+X[0]+"\t"+X[1]+"\t"+X[2]+"\t"
            			+X[3]+"\t"+X[4]+"\t"+X[5]);
            }
            jat.add(mjd_utc,X);
        }
    }

    /**
     * Implements the Derivatives interface for use with the RungeKutta integrator.
     * Given a time in seconds since epoch and a state vector, return the derivatives
     * of the state vector.
     *
     * @param t Time since epoch [sec]
     * @param X State vector. [x y z xdot ydot zdot ... other parameters optional]
     */
    public double[] derivs(double t, double[] x) {
        //* Update spacecraft
        SpacecraftModel sm;
        if(use_formation){
            sm = sc_formation.get_spacecraftmodel(current_sc);
            sm.update(t,x);
        }else{
            sm = this.sc;
            sm.update(t,x);
        }
        spacetime.update(t);
        //* Get non-control derivatives
        double[] xdot = spacetime.derivs(t,sm.get_spacecraft());
        //* Get control derivatives
        double[] xdot2 = sm.control_thrust_derivs(t,xdot);
       
        //NOTE:  as of last revision, this was not correct.  
        //A more thorough was of adding controls may be required
        //* Compile derivatives
        //double[] out = MathUtils.plus(xdot,xdot2);
        double[] out = xdot;
        return out;
    }

    /**
     * Main method to test the class.  For a generic way to run the simulation, see
     * jat.demo.simulation.Simulation.java
     *
     * @see jat.demo.simulation.Simulation
     * @param args
     */
    public static void main(String[] args) {
    	//* Initialize the simulation model
        SimModel sim = new SimModel();
        //* Keep track of the runtime
        double start = System.currentTimeMillis();
        //* Find the proper directory for input/output
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        //* Name of the various spacecraft
        String[] tests = {"demo"};
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag =
        	{{false,true,true,false,false}};
        //* Optional string which appends test numbers to the spacecraft
        String[][] test_nums = 	{{""}};
        //* If plot_traj == true , then plot the output to Celestia
        boolean plot_traj = true;
        int i=0,j=0;
        //* Initial radius [km]
        VectorN r = new VectorN(-4453.783586,-5038.203756,-426.384456);
        //* convert to SI units [m]
        r = r.times(1000);
        //* Initial velocity [km/s]
        VectorN v = new VectorN(3.831888,-2.887221,-6.018232);
        //* convert to SI units [m/s]
        v = v.times(1000);
        //* Simulation parameters t0: initial time [sec]  tf: final time [sec]
        double t0 = 0, tf = 86400;
        //* Modified Julian Date (in UTC) at the initial time
        double mjd_utc = 53157.5;  //* June 1, 2004 12:00 UTC
        //* Simulation stepsize [sec]
        double stepsize = 60;
        //* Output directory (for the textual output data)
        String out = dir+"output"+fs+tests[i]+test_nums[j][i]+".txt";
        //* Coeff. of reflectivity, Drag Coeff., Cross-section [m*m], Mass [kg]
        double cr=1.2, cd=2.2, area=20, mass=1000;
        //* Add the spacecraft parameters to the spacecraft model
        SpacecraftModel sm = new SpacecraftModel(r,v,cr,cd,area,mass);
        //* Loop over the number of spacecraft and test-runs (here only one)
        for(j=0; j<1; j++){
            for(i=0; i<1; i++){
            	//* Initialize the simulation (input simulation parameters)
                sim.initialize(sm,t0,tf,mjd_utc, stepsize, 1, out);
                String test = tests[i]+test_nums[j][i];
                //* Add the force model to the simulation
                sim.initializeForces(force_flag[j], false, test);
                //* Here, add any other custom properties, controllers, or other models
                //* Run the simulation loop
                sim.runloop();
            }
        }
        //* Stop the runtime counter
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        //* Format and output the Celestia files for visualization
        if(plot_traj){
        	jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_Dev/my_celestia/");
        	try{
        		i--;
        		j--;
        		celestia.set_trajectory(sim.get_traj());
        		//* append '_jat' to identify objects originating from JAT
        		String name = tests[i]+test_nums[j][i]+"_jat";
        		celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
        		System.out.println("Wrote to Celestia");
        	}catch(java.io.IOException ioe){}
        }
        System.out.println("Finished");


    }

    /**
     * Initialize the DE405 Ephemerides for use with the Moon.
     *
     */
    public void initializeMoonEphem(){
    	this.spacetime.initializeMoonEphem();
    }

    /**
     * Initialize the DE405 Ephemerides for use with the Sun.
     *
     */
    public void initializeSunEphem(){
    	this.spacetime.initializeSunEphem();
    }
}
