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

import java.util.Random;

import jat.alg.integrators.LinePrinter;
import jat.matlabInterface.MatlabControl;
import jat.matlabInterface.MatlabFunc;
import jat.matvec.data.RandomNumber;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacecraft.SpacecraftModel;
import jat.spacetime.CalDate;
import jat.spacetime.EarthRef;
import jat.spacetime.FitIERS;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.test.propagator.PlotTrajectory;
import jat.traj.RelativeTraj;
import jat.traj.Trajectory;
import jat.util.FileUtil;

/**
 * This is the primary class for a simulation scenario.  It obtains input from files
 * adhering to the specified XML Schema jat.sim.xml.input.sim_input.xsd, creates the necessary
 * objects, and propagates the trajectories.
 * 
 * @see jat.demo.simulation.Simulation
 * @author Richard C. Page III
 *
 */
public class Simulation {
    
    public SimModel sim;
    Random rand;
    
    public Simulation(){
        sim = new SimModel();
        rand = new Random();
    }
        
    public void runSimTwo(){
        SimModel sim = new SimModel();
        double start = System.currentTimeMillis();
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        
        String[] tests = {"ISS","GEO","ISS","GEO","GEO","ISS","Molniya"};
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag = {{false,true,true,true,true},
        		{false,true,true,true,true},
        		{true,false,false,true,false},
        		{true,false,false,false,true},
        		{false,false,false,false,false},
        		{true,false,false,false,false},
        		{false,true,true,true,true}};					//JGM3		0
        String[] test_nums = 
        		{"7_HP","35_HP","4_HP","33","34","1","28_NRL"};  							//JGM3		
        
        boolean plot_traj = true;
        int i=6;
        //*ISS
        //VectorN r = new VectorN(-4453.783586,-5038.203756,-426.384456);
        //VectorN v = new VectorN(3.831888,-2.887221,-6.018232);
        //*Molniya 
        VectorN r = new VectorN(-1529.894287,-2672.877357,-6150.115340);
        VectorN v = new VectorN(8.717518,-4.989709,0);
        //*GEO
        //VectorN r = new VectorN(36607.358256,-20921.723703,0.000000);
        //VectorN v = new VectorN(1.525636,2.669451,0);
        //*GPS VectorN r = new VectorN(5525.33668,-15871.18494,-20998.992446);
        //*GPS VectorN v = new VectorN(2.750341,2.434198,-1.068884);
        //*SunSync VectorN r = new VectorN(-2290.301063,-6379.471940,0);
        //*SunSync VectorN v = new VectorN(-0.883923,0.317338,7.610832);
        
        r = r.times(1000);
        //v = v.times(1755);
        v = v.times(1000);
        double t0 = 0, tf = 3*86400; //604800;
        double mjd_utc = 53157.5;
        //double mjd_utc = 53683;
        double stepsize = 60;
        //String out = dir+"output"+fs+tests[i]+test_nums[j][i]+".txt";
        String out = "C:/Code/Jat/jat/test/propagator/output/"+tests[i]+test_nums[i]+"_jat.txt";
        SpacecraftModel sm = new SpacecraftModel(r,v,1.2,2.2,20,1000);
        //for(j=0; j<1; j++){
          //  for(i=0; i<2; i++){
                sim.initialize(sm,t0,tf,mjd_utc, stepsize, 1, out);
                sim.set_showtimestep(true);
                sim.set_doPrint(true);
                boolean use_JGM2 = false;
                String test = tests[i]+test_nums[i];
                sim.initializeForces(force_flag[i], use_JGM2, test);
                sim.runloop();
            //}	        
        //}
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        plot_traj=false;
        if(plot_traj){
            	        jat.util.Celestia celestia = new jat.util.Celestia("C:/Code/Celestia/");
            	        try{
            	            //i--;
            	            celestia.set_trajectory(sim.get_traj());
            	            String name = tests[i]+test_nums[i];
            	            celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
            	            System.out.println("Wrote to Celestia");
            	        }catch(java.io.IOException ioe){}
            //LinePrinter lp2 = new LinePrinter();
//            RelativeTraj rel = sim.get_rel_traj(lp2);
//            rel.setVerbose(false);
//            double err = rel.get_max_error()*1000;
//            System.out.println("error:  "+err);
//            rel.process();
        }
        String sjat = tests[i]+test_nums[i]+"_jat.txt";
        String sstk = tests[i]+test_nums[i]+"_stk.txt";
        PlotTrajectory.plot(sjat,sstk);
        System.out.println("Finished");    	
 
    }
    
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
                for(int k=0; k<sim.sc_formation.get_num_sc(); k++){
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
    
    public void runSimFormation(){
        SimModel sim = new SimModel();
        double start = System.currentTimeMillis();
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        
        String[] tests = {"ISS","Sun-Sync","GPS","Molniya","GEO"};
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag = 
        {{false,false,false,false,false},								//JGM3		0
                {true,  true,  false,     false,          false},		//Sun		1
                {true,  false,  true,     false,          false},		//Moon		2
                {true,  false, false,     true,           false},		//HP		3
                {true,  false, false,     true,           false},		//NRL		4
                {true,  false, false,     false,          true},		//SRP		5
                {false, true, true, true, true},						//ALL HP	6
                {false, true, true, true, true},						//ALL NRL	7
                {true, false, false, false, false}};					//two body  8
        String[][] test_nums = 
        {{"6","13","20","27","34"},  									//JGM3		
                {"3","10","17","24","31"},  							//Sun		
                {"2","9","16","23","30"},   							//Moon		
                {"4_HP","11_HP","18_HP","25_HP","32_HP"}, 				//HP
                {"4_NRL","11_NRL","18_NRL","25_NRL","32_NRL"},  		//NRL
                {"5","12","19","26","33"},								//SRP 
                {"7_HP","14_HP","21_HP","28_HP","35_HP"},				//All HP
                {"7_NRL","14_NRL","21_NRL","28_NRL","35_NRL"},			//All NRL
                {"1","8","15","22","29"}};								//two body
        
        boolean plot_traj = true;
        int i=0,j=8;
        VectorN r0 = new VectorN(8948.2,  0.0,   0.0);
        VectorN r1 = new VectorN(9198.2,    -4.2328e-012,       -1.075e-012);
        VectorN r2 = new VectorN(9198.2,    -4.2328e-012,       -1.075e-012);
        VectorN v0 = new VectorN(0,           5.8667,	            3.1854);
        VectorN v1 = new VectorN(3.0093e-015,   5.7027,             3.0963);
        VectorN v2 = new VectorN(3.0093e-015,   5.7027,             3.0964);
        double acr = 1.2;
        double acd = 2.2;
        double aa = 20;
        double am = 1000;
        Spacecraft s0 = new Spacecraft(r0.times(1000), v0.times(1000), acr, acd, aa, am);
        Spacecraft s1 = new Spacecraft(r1.times(1000), v1.times(1000), acr, acd, aa, am);
        Spacecraft s2 = new Spacecraft(r2.times(1000), v2.times(1000), acr, acd, aa, am);
        SpacecraftModel[] sc_arg = new SpacecraftModel[3];
        sc_arg[0] = new SpacecraftModel(s0);
        sc_arg[1] = new SpacecraftModel(s1);
        sc_arg[2] = new SpacecraftModel(s2);
        String out = dir+"formation"+fs+"formation.txt";
        
        sim.initialize(sc_arg,0,60000,53157.5,60, 1, out);
        String test = tests[i]+test_nums[j][i];
        boolean use_JGM2 = false;
        sim.initializeForces(force_flag[j], use_JGM2, test);
//        if(plot_traj){
//            String stkfile = "C:/STK_Test_Files/delim/"+tests[i]+test_nums[j][i]+".txt";
//            sim.set_truth_traj(stkfile);
//        }	            
        sim.runloop();	        
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        if(plot_traj){
            jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_dev/celestia/");
            try{
                for(int k=0; k<3; k++){
                    celestia.set_trajectory(sim.get_traj(k));
                    String name = "formation_"+k;
                    celestia.write_trajectory(name,name,sim.mjd_utc_start);
                    System.out.println("Wrote to Celestia "+k);
                }
            }catch(java.io.IOException ioe){}
            LinePrinter lp2 = new LinePrinter();
            RelativeTraj rel = sim.get_rel_traj(lp2,0,1);
            rel.setVerbose(false);
            double err = rel.get_max_error()*1000;
            System.out.println("error:  "+err);
            rel.process();
        }
        System.out.println("Finished");            
    }
    public void runSimBatch(){
        SimModel sim = new SimModel();
        double start = System.currentTimeMillis();
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        
        String[] tests = {"ISS","Sun-Sync","GPS","Molniya","GEO"};
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag = 
        {{false,false,false,false,false},								//JGM3		0
                {true,  true,  false,     false,          false},		//Sun		1
                {true,  false,  true,     false,          false},		//Moon		2
                {true,  false, false,     true,           false},		//HP		3
                {true,  false, false,     true,           false},		//NRL		4
                {true,  false, false,     false,          true},		//SRP		5
                {false, true, true, true, true},						//ALL HP	6
                {false, true, true, true, true},						//ALL NRL	7
                {true, false, false, false, false}};					//two body  8
        String[][] test_nums = 
        {{"6","13","20","27","34"},  									//JGM3		
                {"3","10","17","24","31"},  							//Sun		
                {"2","9","16","23","30"},   							//Moon		
                {"4_HP","11_HP","18_HP","25_HP","32_HP"}, 				//HP
                {"4_NRL","11_NRL","18_NRL","25_NRL","32_NRL"},  		//NRL
                {"5","12","19","26","33"},								//SRP 
                {"7_HP","14_HP","21_HP","28_HP","35_HP"},				//All HP
                {"7_NRL","14_NRL","21_NRL","28_NRL","35_NRL"},			//All NRL
                {"1","8","15","22","29"}};								//two body
        
        boolean plot_traj = true;
        int i=0,j=8;
        VectorN r[] = new VectorN[5];
        VectorN v[] = new VectorN[5];
//      *ISS
        r[0] = new VectorN(-4453.783586,-5038.203756,-426.384456).times(1000);
        v[0] = new VectorN(3.831888,-2.887221,-6.018232).times(1000);
        //* SunSync
        r[1] = new VectorN(-2290.301063,-6379.471940,0).times(1000);
        v[1] = new VectorN(-0.883923,0.317338,7.610832).times(1000);
        //* GPS
        r[2] = new VectorN(5525.33668,-15871.18494,-20998.992446).times(1000);
        v[2] = new VectorN(2.750341,2.434198,-1.068884).times(1000);
        //* Molniya
        r[3] = new VectorN(-1529.894287,-2672.877357,-6150.115340).times(1000);
        v[3] = new VectorN(8.717518,-4.989709,0).times(1000);
        //* GEO
        r[4] = new VectorN(36607.358256,-20921.723703,0.000000).times(1000);
        v[4] = new VectorN(1.525636,2.669451,0).times(1000);
        
        
        double acr = 1.2;
        double acd = 2.2;
        double aa = 20;
        double am = 1000;
        
        double t0 = 0;
        double[] tf = {.01*86400,.01*86400,0.02*86400,.03*86400,.07*86400}; //604800;
        double[] stepsize = {5,5,60,5,60};
        double mjd_utc = 53157.5;
        //double mjd_utc = 53683;
        //String out = dir+"output"+fs+tests[i]+test_nums[j][i]+".txt";        
        for(j=6; j<8; j++){
            for(i=0; i<5; i++){
            	String out = "C:/Code/Jat/jat/test/propagator/output/"+tests[i]+test_nums[j][i]+"_jat.txt";
        		SpacecraftModel sm = new SpacecraftModel(new VectorN(r[i]),new VectorN(v[i]),acr,acd,aa,am);
                sim.initialize(sm,t0,tf[i],mjd_utc, stepsize[i], 1, out);
                sim.set_showtimestep(true);
                boolean use_JGM2 = false;
                String test = tests[i]+test_nums[j][i];
                sim.initializeForces(force_flag[j], use_JGM2, test);
                sim.runloop();
              String sjat = tests[i]+test_nums[j][i]+"_jat.txt";
              String sstk = tests[i]+test_nums[j][i]+".txt";
              PlotTrajectory.plot(sjat,sstk);
            }	        
        }
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        
        System.out.println("Finished");    	
    }

    public void test(){
        jat.matlabInterface.MatlabControl test = new MatlabControl();
        test.eval("disp('testing')");
    }
    
    public void runSim4Datsim(){
        SimModel sim = new SimModel();
        double start = System.currentTimeMillis();
        String fs = FileUtil.file_separator();
        String dir = FileUtil.getClassFilePath("jat.sim", "SimModel");
        
        //String[] tests = {"jat_ecef"};
        String[] tests = {"jat_j2k"};
        //* force_flag = {2-Body, Sun,   Moon, Harris Priester, Solar Radiation}
        boolean[][] force_flag = {{false,true,true,false,true},
        						  {false,true,true,false,true},
        						  {false,false,false,true,false}};					//JGM3		0
        String[][] test_nums = 
        		{{"1_8","1_1","1_6"}};  							//JGM3		
        double ptsol_noise = 3.0;
        int i=0,j=0;
        //*ISS
        VectorN r[] = new VectorN[3];
        VectorN v[] = new VectorN[3];
        r[0] = new VectorN(-1.86631428726119e+06,  4.21156716602933e+07, -8.23659306996206e+03);
        v[0] = new VectorN(-3.07220560283457e+03, -1.35729347644219e+02,  1.35175509960231e+01);
        r[1] = new VectorN(-4.02622553013133e+007,-1.24965505638999e+007,-4.27951891301160e+003);
        v[1] = new VectorN(9.11185854997334e+002,-2.93710765139044e+003,1.39092766633154e+001);
        r[2] = new VectorN(2.72701600600060e+06,-4.56452737925970e+06,-4.35529550517323e+06);
        v[2] = new VectorN(4.45912406220447e+03,-3.63297132340426e+03,6.59775392945695e+03);
        
        double cd = 2.2, cr = 1.2, m = 2100, area = 54;
        double t0 = 0, tf = 4*86400;
        double mjd_utc[] = {50985,52187,51013.0};  
        double stepsize = 60;
        for(i=0; i<1; i++){
        	String out = dir+"output"+fs+tests[j]+test_nums[j][i]+".txt";
        	SpacecraftModel sm = new SpacecraftModel(r[i],v[i],cr,cd,area,m);
        	sim.initialize(sm,t0,tf,mjd_utc[i], stepsize, 1, out);
        	sim.set_showtimestep(true);
        	boolean use_JGM2 = false;
        	String test = tests[j]+test_nums[j][i];
        	sim.initializeForces(force_flag[i], use_JGM2, test);
        	//sim.set_doPrint(false);
        	if(i>0){
        		sim.runloop();
        		Trajectory traj = sim.get_traj();
//        		FitIERS iers = new FitIERS();
//                iers.process();
//                Time time = new Time(mjd_utc[i]);
//        		EarthRef earth = new EarthRef(time);
//        		double[] param = iers.search(time.mjd_tt());
//    		    //earth.setIERS(param[0],param[1]);
//    		    //time.set_UT1_UTC(param[2]);
//        		Trajectory traj_ecf = new Trajectory();
//        		for(int k=0; k<traj.size(); k++){
//        			double t = traj.getTimeAt(k);
//        			VectorN x = new VectorN(traj.getState(k));
//        			time.update(t);
//        			//param = iers.search(time.mjd_tt());
//        		    //earth.setIERS(param[0],param[1]);
//        		    //time.set_UT1_UTC(param[2]);
//        			VectorN xecf = earth.eci2ecf(x.get(0,3),x.get(3,3),new Time(t));
//        			traj_ecf.add(t,xecf.x);
//        		}
        		//parse4datsim(traj_ecf,out);
        		parse4datsim(traj,out);
        	}else{
        		Trajectory meas = new Trajectory();        
        		double t = t0;
        		double mjd = mjd_utc[i] + t/86400.0;
        		VectorN x = new VectorN(sim.sc.get_spacecraft().toStateVector());
        		Time time = new Time(mjd_utc[i]);
        		EarthRef earth = new EarthRef(time);
        		FitIERS iers = new FitIERS();
                //iers.process(); //* don't need this               
        		double[] param = iers.search(time.mjd_tt());
    		    earth.setIERS(param[0],param[1]);
    		    time.set_UT1_UTC(param[2]);
        		//earth.update(time);
        		RotationMatrix eci2ecf = new RotationMatrix(earth.eci2ecef(time));
        		VectorN xecef = eci2ecf.transform(x.get(0,3));
        		//meas.add(mjd,make_ptsol(xecef.x,ptsol_noise));
        		meas.add(mjd,make_ptsol(x.x,ptsol_noise));
        		while(t<tf){
        			sim.step(stepsize);
        			t=t+stepsize;
        			mjd = mjd_utc[i] + t/86400.0;
        			x = new VectorN(sim.sc.get_spacecraft().toStateVector());
        			time.update(t);
        			//param = iers.search(time.mjd_tt());
        		    //earth.setIERS(param[0],param[1]);
        		    //time.set_UT1_UTC(param[2]);        			        		    
        		    //earth.update(time);
            		eci2ecf = new RotationMatrix(earth.eci2ecef(time));
            		xecef = eci2ecf.transform(x.get(0,3));
            		//meas.add(mjd,make_ptsol(xecef.x,ptsol_noise));
            		meas.add(mjd,make_ptsol(x.x,ptsol_noise));
        		}
        		parsePtSol(meas,dir+"output"+fs+"test1_8_jat_eci.rnx");
        	}
        }	        
        double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
        boolean plot_traj = false;
        if(plot_traj){
            	        jat.util.Celestia celestia = new jat.util.Celestia("C:/games/Celestia_Dev/my_celestia/");
            	        try{
            	            i--;
            	            j--;
            	            celestia.set_trajectory(sim.get_traj());
            	            String name = tests[i]+test_nums[j][i];
            	            celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
            	            System.out.println("Wrote to Celestia");
            	        }catch(java.io.IOException ioe){}
            //LinePrinter lp2 = new LinePrinter();
//            RelativeTraj rel = sim.get_rel_traj(lp2);
//            rel.setVerbose(false);
//            double err = rel.get_max_error()*1000;
//            System.out.println("error:  "+err);
//            rel.process();
        }
        System.out.println("Finished");    	
 
    }
    
    private double[] make_ptsol(double[] state,double noise) {
    	double[] out = new double[3];
    	double tmp;
		for(int i=0; i<3; i++){
			tmp = rand.nextGaussian()*noise*noise;
			out[i] = state[i] + tmp;
		}
		return out;
	}

	private void parse4datsim(Trajectory traj, String filename){
    	LinePrinter lp = new LinePrinter(filename);
    	lp.println(""+filename+" Ephemeris file");
    	lp.println("1998 6  21  0  0      0.00000     60.00000 YYYY MM DD HH MM SS Interval");
    	lp.println("1998 172      0.00000 YYYY DOY SEC");
    	lp.println(" 50985.00000 Modified Julian date");
    	lp.println(" MJD         Seconds    Position (meters)                                                 Velocity (meter/second)");
    	for(int i=0; i<traj.size(); i++){
    		double time = traj.getTimeAt(i);
    		double day = Math.floor(time);
    		String sday = String.valueOf((int)day);
    		double sec = (time - day)*86400.0;
    		String ssec = String.valueOf(Math.round(sec));//new jat.util.PrintfFormat("%12.6f").sprintf(sec);
    		String[] state = new String[6];
    		for(int j=0; j<6; j++)
    			state[j] = String.valueOf(traj.getState(i)[j]*1000);//new jat.util.PrintfFormat("%-23.14e").sprintf(traj.getState(i)[j]*1000.0);
    		lp.println(" "+sday+" "+ssec+"  "+state[0]+"  "+state[1]+"  "+state[2]+"  "+state[3]+"  "+state[4]+"  "+state[5]);
    	}
    	lp.println("");
    	lp.close();
    }
	private void parsePtSol(Trajectory traj, String filename){
    	LinePrinter lp = new LinePrinter(filename);
    	for(int i=0; i<traj.size(); i++){
    		double t = traj.getTimeAt(i);
    		CalDate time = new CalDate(t);
    		int year = time.year();
    		String syear = String.valueOf(year);//new jat.util.PrintfFormat("%4.0d").sprintf(year);
    		int doy = time.doy();
    		String sdoy = String.valueOf(doy);//new jat.util.PrintfFormat("%4.0d").sprintf(doy);
    		double sec = Math.round(time.sec_of_day());
    		String ssec = String.valueOf(sec);//new jat.util.PrintfFormat("%13.6f").sprintf(sec);
    		String[] state = new String[6];
    		double clock = 0.0;
    		String misc = String.valueOf(clock);//new jat.util.PrintfFormat("%13.10f").sprintf(clock);
    		for(int j=0; j<3; j++)
    			state[j] = String.valueOf(traj.getState(i)[j]);//new jat.util.PrintfFormat("%-13f").sprintf(traj.getState(i)[j]);
    		lp.println(""+syear+"  "+sdoy+"  "+ssec+"         "+state[0]+"         "+state[1]+"         "+state[2]+"         "+misc);
    	}
    	lp.println("");
    	lp.close();
    }
    public static void main(String[] args) throws InterruptedException {
        Simulation sim = new Simulation();
        //sim.runSim4Datsim();
        sim.runSimTwo();
        //sim.runSimBatch();
        //sim.runSimMatlab();
        //sim.runSimFormation();
    }
}
