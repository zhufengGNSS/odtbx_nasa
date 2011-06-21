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
package jat.demo.simulation;

import jat.sim.*;
import jat.alg.integrators.LinePrinter;
import jat.matlabInterface.MatlabControl;
import jat.matlabInterface.MatlabFunc;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacecraft.SpacecraftModel;
import jat.spacetime.EarthRef;
import jat.traj.RelativeTraj;
import jat.util.FileUtil;

/**
 * This is a demo class for a generic simulation scenario.  It creates the necessary
 * objects, propagates the trajectories, and prints the trajectories into the 'Celestia'
 * format for visualization.  'Celestia' can be obtained at: 
 * @link http://celestia.sourceforge.net
 * 
 * @author Richard C. Page III
 *
 */
public class Simulation {
    
	/**
	 * Simulation Model created as the worker object to run the simulation
	 */
    public SimModel sim;
    
    /**
     * Simple constructor which instantiates the simulation model.
     */
    public Simulation(){
        sim = new SimModel();
    }

    /**
     * This method is an example of using a hardcoded ephemeris to intialize the
     * simulation.  Values for initial radius, velocity, and relevant spacecraft
     * characteristics are added to the model.  The SimModel class then does the 
     * work propagating the orbit.  Finally, the output is written to a Celestia
     * spacecraft (*.ssc) and trajectory file (*.xyz). 
     *
     */
    public void runSimManual(){
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
        		System.out.println("Wrote " + name + " to Celestia");
        	}catch(java.io.IOException ioe){}
        }
        System.out.println("Finished");    	
 
    }
    
    /**
     * This method shows how to use Matlab to initialize a simulation run by JAT
     * @throws InterruptedException
     */
    public void runSimMatlab() throws InterruptedException{
        MatlabControl input = new MatlabControl();
        input.eval("disp('runSimMatlab init')");
        SimModel sim = new SimModel();
        //double start = System.currentTimeMillis();
        //String fs = FileUtil.file_separator();
        //String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag = 
        {{false,false,false,false,false},						//JGM3		0
                {true,  true,  false,     false,          false},		//Sun		1
                {true,  false,  true,     false,          false},		//Moon		2
                {true,  false, false,     true,           false},		//HP		3
                {true,  false, false,     true,           false},		//NRL		4
                {true,  false, false,     false,          true},		//SRP		5
                {false, true, true, true, true},						//ALL HP	6
                {false, true, true, true, true},						//ALL NRL	7
                {true, false, false, false, false}};					//two body  8
        
        
        int force_case = 8;
        boolean plot_traj = true;
        String name = "Matlab";
        input.eval("disp('runSimMatlab initMatlab')");
        sim.initializeMatlab("initJAT_sc2","initJAT_integ",
                "output"+"/"+name+".txt");
        input.eval("disp('runSimMatlab initMatlab finished')");
        boolean use_JGM2 = false;
        String test = name;
        sim.initializeForces(force_flag[force_case], use_JGM2, test);
        input.eval("disp('runSimMatlab initMatlab finished')");
        sim.runloop();
        //double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        //System.out.println("Elapsed time [min]: "+elapsed);
        VectorN K  = new VectorN(0,0,1);
        double[] X0 = sim.get_traj(1).get(0);
        VectorN r = new VectorN(X0[1],X0[2],X0[3]);
        r.times(1000);
        r.print("position: ");
        VectorN v = new VectorN(X0[4],X0[5],X0[6]);
        v.times(1000);
        v.print("velocity: ");
        VectorN hv = r.crossProduct(v);
        VectorN nv = K.crossProduct(hv);
        double n  = Math.sqrt(nv.mag()*nv.mag());
        double h2 = (hv.mag()*hv.mag());
        double v2 = v.mag()*v.mag();
        double rmag  = r.mag();
        System.out.println("rmag: "+rmag);
        //ev = 1/EarthRef.GM_Earth *( (v2-EarthRef.GM_Earth/r)*rv - (rv'*vv)*vv );
        VectorN ev = (r.times(v2-EarthRef.GM_Earth/rmag).minus(v.times(r.dotProduct(v)))).times(1.0/EarthRef.GM_Earth);
        double p  = h2/EarthRef.GM_Earth;
        double e = ev.mag();
        double a = p/(1-e*e);
        double period = 2*Math.PI*Math.sqrt(a*a*a/EarthRef.GM_Earth);
        System.out.println("Period [hr]: "+period+"  a: "+a+"  e: "+e);
        
        if(plot_traj){
            jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_dev/celestia/");
            try{
                for(int k=0; k<sim.get_number_of_spacecraft(); k++){
                    celestia.set_trajectory(sim.get_traj(k));
                    String celname = "formation_"+k;
                    celestia.write_trajectory(celname,celname,sim.mjd_utc_start);
                    System.out.println("Wrote to Celestia "+k);
                }
            }catch(java.io.IOException ioe){}
            LinePrinter lp2 = new LinePrinter();
            RelativeTraj rel = sim.get_rel_traj(lp2,0,1);
            rel.setVerbose(false);
            rel.process();
            rel = sim.get_rel_traj(lp2,0,2);
            rel.setVerbose(false);
            rel.process();
        }
        System.out.println("Finished");
    }
    
    /**
     * Main method.
     * @param args
     * @throws InterruptedException
     */
    public static void main(String[] args) throws InterruptedException {
        Simulation sim = new Simulation();
        sim.runSimManual();
        //sim.runSimXML();
        //sim.runSimMatlab();
    }
}
