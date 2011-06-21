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
 */
package jat.demo.simulation;

import jat.matvec.data.VectorN;
import jat.sim.SimModel;
import jat.spacecraft.SpacecraftModel;
import jat.util.FileUtil;
import jat.cm.*;
import jat.forces.*;

/**
 * TODO Javadoc
 * @author Richard C. Page III
 */
public class EarthMoonSim {

	
    public EarthMoonSim(){ }
    
    /**
     * Simulation of an Earth-Moon Free return maneuver from LEO
     *
     */
    public static void runSimulation(){
    	SimModel sim = new SimModel();
    	sim.set_showtimestep(true);
    	//* Simulation parameters t0: initial time [sec]  tf: final time [sec]
        double t0 = 0, t1 = 0.0505*86400, t2 = 6*86400; 
        //* Modified Julian Date (in UTC) at the initial time
        double mjd_utc = 53157.5;  //* June 1, 2004 12:00 UTC
        //* Simulation stepsize [sec]
        double stepsize = 1*60; 
    	//* Keep track of the runtime
        double start = System.currentTimeMillis();
        //* Find the proper directory for input/output
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        //* Name of the various spacecraft
        String[] tests = {"demo"};
        //* Optional string which appends test numbers to the spacecraft
        String[][] test_nums = 	{{""}};  									
        //* If plot_traj == true , then plot the output to Celestia
        boolean plot_traj = true;
        int i=0,j=0;
        GravitationalBody earth = new GravitationalBody(Constants.GM_Earth);
        Moon moon = new Moon(mjd_utc);
        Sun sun = new Sun();
        sim.add_ForceModel(earth);
        sim.add_ForceModel(moon);
        sim.initializeMoonEphem();
        sim.add_ForceModel(sun);
        sim.initializeSunEphem();
        //* Initial radius [km]
        VectorN rold = new VectorN(-4453.783586,-5038.203756,-426.384456);
        VectorN r = moon.getPosition().unitVector().times(rold.mag());
        //* convert to SI units [m]
        r = r.times(1000);
        //* Initial velocity [km/s]
        VectorN vold = new VectorN(3.831888,-2.887221,-6.018232);
        VectorN v = moon.getVelocity().unitVector().times(vold.mag());
        //* convert to SI units [m/s]
        v = v.times(1000);
        
        //* Output directory (for the textual output data)
        String out = dir+"output"+fs+tests[i]+test_nums[j][i]+".txt";
        //* Coeff. of reflectivity, Drag Coeff., Cross-section [m*m], Mass [kg]
        double cr=1.2, cd=2.2, area=20, mass=1000;
        //* Add the spacecraft parameters to the spacecraft model
        SpacecraftModel sm = new SpacecraftModel(r,v,cr,cd,area,mass);
        	//* Initialize the simulation (input simulation parameters)
                sim.initialize(sm,t0,t1,mjd_utc, stepsize, 1, out);
                String test = tests[i]+test_nums[j][i];
                //* Add the force model to the simulation
                //* Here, add any other custom properties, controllers, or other models
                //* Run the simulation loop
                sim.runloop();
//                if(plot_traj){
//                	jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_Dev/my_celestia/");
//                	try{
//                		//i--;
//                		//j--;
//                		celestia.set_trajectory(sim.get_traj());
//                		//* append '_jat' to identify objects originating from JAT
//                		String name = tests[i]+"_1_jat";
//                		celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
//                		System.out.println("Wrote to Celestia");
//                	}catch(java.io.IOException ioe){}
//                }
                sim.t0=t1;
                sim.tf=t2;
                double dvmag = 3100;
                //VectorN dv = new VectorN(sim.sc.get_abs_vel().times(0.419));
                VectorN dv = new VectorN(sim.sc.get_abs_vel().times(0.41905));
                //dv.x[2] = dv.x[2]-0;
                //dv = (dv.unitVector().times(dvmag));
                sm.applyDeltaV(dv);
                sim.runloop();
                System.out.println("r0_mag: "+r.mag()+" m");
                System.out.println("v0_mag: "+v.mag()+" m/s");
                System.out.println("dv_mag: "+dv.mag()+" m/s/s");
                
        //* Stop the runtime counter
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        //* Format and output the Celestia files for visualization
        if(plot_traj){
        	jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_Dev/my_celestia/");
        	try{
        		//i--;
        		//j--;
        		celestia.set_trajectory(sim.get_traj());
        		//* append '_jat' to identify objects originating from JAT
        		String name = tests[i]+"_3_jat";
        		celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
        		System.out.println("Wrote to Celestia");
        	}catch(java.io.IOException ioe){}
        }
        System.out.println("Finished");    	
 
    }
    
	public static void main(String[] args) {
		runSimulation();
	}
}
