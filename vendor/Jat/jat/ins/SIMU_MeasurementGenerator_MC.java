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
 * File Created on Aug 29, 2003
 */
 
package jat.ins;
//import jat.gps_ins.*;
import jat.cm.*;
//import jat.plot.*;
import jat.traj.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.math.*;
//import jat.forces.*;
import jat.timeRef.*;

/** 
 * SIMU_MeasurementGenerator generates SIMU measurements and truth trajectories for
 * the chaser and ISS for a Monte Carlo study.
 * 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */ 
public class SIMU_MeasurementGenerator_MC {
	public Trajectory chaser = new Trajectory();
	
	public Trajectory iss = new Trajectory();
	
	private ChaserEOM chaserEOM = new ChaserEOM(chaser, 51969.0, 104328.0, 454.4, 2.0);
	private ChaserBurnEOM burnEOM = new ChaserBurnEOM(chaser, 51969.0, 104328.0, 454.4, 2.0);
	
	private JGM3DragEOM issEOM = new JGM3DragEOM(iss, 51969.0, 128990.0, 640.7, 2.35);
	
	private String title;
	
	private String chaserFile;
	
	private String issFile;
	
	private String insMeasFile;
	
	private double dt = 0.5;
		
	private RungeKutta8 rk8 = new RungeKutta8(this.dt);
		
	public FiniteBurnList burnlist;
	private int dv_index = 0;
	
	private INS_MeasurementList list = new INS_MeasurementList();
	private double tau = 3.0 * 3600.0;
	private VectorN gyro_sigma;
	private VectorN accel_sigma;
	private GaussianVector gyro_noise;
	private GaussianVector accel_noise;
	private Matrix accel_sfma;  // later
	private Matrix gyro_sfma;   // later
	
	private LinePrinter biaslp;
	
	private long seed = -1;
	
	
	private Matrix phi;

	/**
	 * Constructor
	 * @param t String containing the title
	 * @param file1 String containing the chaser trajectory file to be created
	 * @param file2 String containing the ISS trajectory file to be created
	 * @param file3 String containing the finite burn file name to be input
	 * @param file4 String containing the INS measurement file to be created
	 * @param file5 String containing the true bias file to be created
	 */			
	public SIMU_MeasurementGenerator_MC(String t, String file1, String file2, String file3, String file4, String file5, long sd) {
		this.seed = sd;
		this.chaserFile = file1;
		this.issFile = file2;
		this.title = t;
		this.burnlist = FiniteBurnList.recover(file3);
		this.insMeasFile = file4;
		this.biaslp = new LinePrinter(file5);
		
		// set up stochastic stuff
		this.phi = this.phi(this.dt, this.tau);
		this.accel_sigma = this.qAccel(this.dt).ebeSqrt();
		this.gyro_sigma = this.qGyro(this.dt).ebeSqrt();
		VectorN zero = new VectorN(6);
		this.accel_noise = new GaussianVector(zero, this.accel_sigma, seed);
		this.gyro_noise = new GaussianVector(zero, this.gyro_sigma, (seed-1));
		
		// set scalefactor-misalignment matrix
		this.accel_sfma = this.sfmaAccel();
		this.gyro_sfma = this.sfmaGyro();
		
	}
					
	private Matrix phi(double dt, double tau) {
		Matrix out = new Matrix(6);
		double exponent = -1.0 * dt / tau;
		double term1 = Math.exp(exponent);
		double term2 = (term1 - 1.0) * tau;
		out.set(0, 0, term1);
		out.set(1, 1, term1);
		out.set(2, 2, term1);
		out.set(3, 0, term2);
		out.set(4, 1, term2);
		out.set(5, 2, term2);
		return out;
	}
	
	private VectorN qGyro(double dt){
		
		double biasSigma = 5.0E-05 * MathUtils.DEG2RAD / 3600.0;
		double exponent = -2.0*dt/tau;
		double biasQ = biasSigma*biasSigma*(1.0 - Math.exp(exponent));
		
		double noiseSigma = 7.9E-05 * MathUtils.DEG2RAD;
		double noiseQ = noiseSigma * noiseSigma*dt / 3600.0;
		VectorN out = new VectorN(6);
		for (int i = 0; i < 3; i++) {
			out.set(i, biasQ);
			out.set(i+3, noiseQ);
		}
		return out;
	}
	
	private Matrix sfmaGyro(){
		long sd = this.seed - 2;
		RandomNumber rn = new RandomNumber(sd);
		double sigma_ma = 1.0 * MathUtils.ARCSEC2RAD;  // 1 arcsec
		double sigma_sf = 1.0E-06;                // 1 parts per million
		Matrix temp = new Matrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double x = rn.normal(0.0, sigma_ma);
				if (i == j) x = rn.normal(0.0, sigma_sf);
				temp.set(i, j, x);
			}
		}

		Matrix eye = new Matrix(3);
		Matrix out = eye.minus(temp);
		return out;
				
	}

	private Matrix sfmaAccel(){
		long sd = this.seed - 3;
		RandomNumber rn = new RandomNumber(sd);
		double sigma_ma = 2.5E-05;  // 100 microrad
		double sigma_sf = 5.0E-05;                // 310 parts per million
		Matrix temp = new Matrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double x = rn.normal(0.0, sigma_ma);
				if (i == j) x = rn.normal(0.0, sigma_sf);
				temp.set(i, j, x);
			}
		}
		Matrix eye = new Matrix(3);
		Matrix out = eye.minus(temp);
		return out;
				
	}

	
	private VectorN qAccel(double dt){
		double biasSigma = 1.0E-09;
		double exponent = -2.0*dt/tau;
		double biasQ = biasSigma*biasSigma*(1.0 - Math.exp(exponent));

		double noiseSigma = 1.0E-10;
		double noiseQ = noiseSigma * noiseSigma * dt;
		VectorN out = new VectorN(6);
		for (int i = 0; i < 3; i++) {
			out.set(i, biasQ);
			out.set(i+3, noiseQ);
		}
		return out;
	}
	
	private VectorN initBiases(){
		long sd = this.seed - 4;
		RandomNumber rn = new RandomNumber(sd);
		VectorN out = new VectorN(12);
		double gyro_sigma = 0.0003 * MathUtils.DEG2RAD / 3600.0;
		double accel_sigma = 1.0E-06 * 9.81;
		for (int i = 0; i < 3; i++) {
			double xa = rn.normal(0.0, accel_sigma);
			double xg = rn.normal(0.0, gyro_sigma);
			out.set(i, xa);
			out.set((i+6), xg);
		}
		return out;
	}
		
	/**
	 * Generate the data
	 * @param t0 initial time (sim time in seconds)
	 * @param tf final time (sim time in seconds)
	 * @param x0 VectorN containing the initial chaser state
	 * @param xref0 VectorN containing the initial ISS state
	 */					
	public void generate (double t0, double tf, VectorN x0, VectorN xref0) {
		
		// propagate the orbit
		int neqns = x0.length;

        double t = t0;
        double[] newchaser = new double[neqns];
        double[] oldchaser = new double[neqns];
        double[] newiss = new double[neqns];
        double[] oldiss = new double[neqns];

        // put initial conditions into the previous state array

        for (int i = 0; i < neqns; i++) {
            oldchaser[i] = x0.x[i];
            oldiss[i] = xref0.x[i];            
        }

        if ((t + dt) > tf) {
            dt = tf - t;
        }

        chaserEOM.print(t, oldchaser);
        issEOM.print(t, oldiss);
        
        boolean burnstarted = false;
        
        // initialize INS states, add an init function to start with non-zero biases
        VectorN insStates = this.initBiases();
		this.biaslp.println(t+"\t"+insStates.toString());
                
        // main integration loop
        while (t < tf) {
        	
        	// reset accelerometer and gyro accumulators
        	for (int i = 0; i < 6; i++) {
        		oldchaser[i+10] = 0.0;
        	}
        	        	
        	// check next burn time
        	FiniteBurn burn = burnlist.get(dv_index);
        	double tstart = burn.tstart;
        	double tstop = burn.tstop;        	        	
        	if ((t >= tstart) && (t < tstop)) {
        		burnEOM.setBurn(burn);
        		newchaser = rk8.step(t, oldchaser, burnEOM);
        		burnstarted = true;
        	}
        	else {
        		newchaser = rk8.step(t, oldchaser, chaserEOM);
        		if ((t > tstop) && (burnstarted)) {
        			System.out.println(t+" burn "+dv_index+" completed, incrementing");
        			if (burnlist.hasNext(dv_index)) {
        				 dv_index = dv_index + 1;
        			}
        			burnstarted = false;
        		}
        	}
                        
            newiss = rk8.step(t, oldiss, issEOM);
            
            // normalize the quaternion
            VectorN r = new VectorN(newchaser[0], newchaser[1], newchaser[2]);
            VectorN v = new VectorN(newchaser[3], newchaser[4], newchaser[5]);
            RSW_Frame rsw = new RSW_Frame(r, v);
            Matrix x = rsw.ECI2RSW().transpose();
            Quaternion q = new Quaternion(x);
            // check for sign flip
            for (int i = 0; i < 4; i++) {
            	double mult = newchaser[i+6] * q.x[i];
            	newchaser[i+6] = q.x[i];
            	if ((mult < 0.0)&&(Math.abs(newchaser[i+6]) > 1.0E-12)) newchaser[i+6] = -1.0*q.x[i];
            }
                                    
            for (int i = 0; i < neqns; i++) {
	            oldchaser[i] = newchaser[i];
	            oldiss[i] = newiss[i];
            }
                       
            t = t + dt;
            
            insStates = this.INSmeas(t, oldchaser, insStates);
            
            if (MathUtils.Frac(t) == 0.0){
		        chaserEOM.print(t, oldchaser);
		        issEOM.print(t, oldiss);
            }

            if ((t + dt) > tf) {
                dt = tf - t;
            }

        }


		// set the attributes
		chaser.setTitle(this.title);
		chaser.setCentralBody(CentralBody.EARTH);
		chaser.setCoordinateSystem(CoordinateSystem.INERTIAL);
		chaser.setDistanceUnits(DistanceUnits.METERS);
		chaser.setEpoch(2001, 3, 1, 0, 0, 0.0);
		chaser.setTimeUnits(TimeUnits.SECONDS);
		String[] labels = {"t","x","y","z","xdot","ydot","zdot","q1","q2","q3","q4","fx","fy","fz"};
		chaser.setLabels(labels);

		iss.setTitle(this.title);
		iss.setCentralBody(CentralBody.EARTH);
		iss.setCoordinateSystem(CoordinateSystem.INERTIAL);
		iss.setDistanceUnits(DistanceUnits.METERS);
		iss.setEpoch(2001, 3, 1, 0, 0, 0.0);
		iss.setTimeUnits(TimeUnits.SECONDS);
		iss.setLabels(labels);


		// Serialize the trajectory
		chaser.serialize(this.chaserFile);
		iss.serialize(this.issFile);
		list.serialize(this.insMeasFile);
		System.out.println("trajectory serialized");
		this.biaslp.close();

	}
	
	private VectorN INSmeas(double t, double[] state, VectorN in) {
        // get the true accel and omega
        VectorN accel_true = new VectorN(state[10], state[11], state[12]);
        VectorN omega_true = new VectorN(state[13], state[14], state[15]);
        
        // strip out the bias states
        VectorN accel_bias = in.get(0, 3); 
        VectorN gyro_bias = in.get(6, 3);
        
        // form the state vectors, zeroing out the accel and gyro accumulators
        VectorN zero = new VectorN(3);
        VectorN accel_state = accel_bias.append(zero);
        VectorN gyro_state = gyro_bias.append(zero);
       
        // get the next noise values 
        this.accel_noise.nextSet();
        this.gyro_noise.nextSet(); 
        
        // propagate the states forward in time
        accel_state = this.phi.times(accel_state).plus(this.accel_noise);
        gyro_state = this.phi.times(gyro_state).plus(this.gyro_noise);
        
        // get the noise accumulators
        VectorN accel_noise_acc = accel_state.get(3, 3);
        VectorN gyro_noise_acc = gyro_state.get(3, 3);
        
        // add in scalefactor-misalignment
        VectorN accel_temp = this.accel_sfma.times(accel_true);
        VectorN gyro_temp = this.gyro_sfma.times(omega_true);
        
        
        // form the measurements        
        VectorN accel_meas = accel_temp.plus(accel_noise_acc);
        VectorN gyro_meas = gyro_temp.plus(gyro_noise_acc);
        
        // add the measurements to the list
//        INS_Measurement meas = new INS_Measurement(t, accel_meas, gyro_meas);
        INS_Measurement meas = new INS_Measurement(t, accel_true, omega_true);
        list.add(meas);
//    	System.out.println(t+"\t"+accel_meas+"\t"+gyro_meas);

		// send the biases to the lineprinter
		if (MathUtils.Frac(t) == 0.0){
			this.biaslp.println(t+"\t"+accel_state.toString()+gyro_state.toString());
		}

		// form the output state vector
		VectorN out = accel_state.append(gyro_state);
		return out;
		
	}
	

	/** Runs the example.
	 * @param args Arguments.
	 */
	public VectorN initChaser () {

		double d = 15.0/6765.5000;
		double theta = 360.0 - d * Constants.rad2deg;
		
		TwoBody sat =
			new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, theta);
		
		VectorN r = sat.getR();
		VectorN v = sat.getV();
		VectorN rv = new VectorN(r, v);
		
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
						
		Quaternion q0 = new Quaternion(Cb2i);
		VectorN x0 = new VectorN(rv, q0);
		VectorN zero = new VectorN(6);
		x0 = x0.append(zero);
		
		return x0;
	}

	public VectorN initISS () {

		// create a TwoBody orbit using orbit elements
		TwoBody sat =
			new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, 0.0);
				
		VectorN r = sat.getR();
		VectorN v = sat.getV();
		VectorN rv = new VectorN(r, v);
				
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
		Quaternion q0 = new Quaternion(Cb2i);
		VectorN x0 = new VectorN(rv, q0);
		VectorN zero = new VectorN(6);
		x0 = x0.append(zero);
		
		return x0;
	}

	
	public static void main(String[] args){
		long counter = 1;
		int j = 1;
		for (int i = 0; i < 30; i++) {
			
			String num = Integer.toString(j);
			
			String dotJat = ".jat";
			String dotTxt = ".txt";
		
			// set up the generator
			String title = "STS/ISS Rendezvous Trajectory";
			String file1 = "C:\\Jat\\jat\\output\\monte\\rvtraj_simu_rbar.jat";
			String file2 = "C:\\Jat\\jat\\output\\monte\\isstraj_simu_rbar.jat";
			String file3 = "C:\\Jat\\jat\\input\\burns\\rbar_burns.jat";
			String file4 = "C:\\Jat\\jat\\input\\monte\\simu_rbar_" + num + dotJat;
			String file5 = "C:\\Jat\\jat\\output\\monte\\simu_rbar_bias_"+ num + dotTxt;
	//		String file4 = "C:\\Jat\\jat\\input\\sigi_ins.jat";
	//		String file5 = "C:\\Jat\\jat\\output\\sigi_bias.txt";
	
			System.out.println("Generating: "+file4);

			long seed = -1 * counter;
			
			SIMU_MeasurementGenerator_MC x = new SIMU_MeasurementGenerator_MC(title, file1, file2, file3, file4, file5, seed);
			VectorN x0 = x.initChaser();
			VectorN xref0 = x.initISS();
			
			
			double t0 = 0.0;
	//		double tf = 5579.0;
			double tf = 9400.0;
			x.generate(t0, tf, x0, xref0);
			System.out.println("Done Generating");
			
//			LinePrinter lp1 = new LinePrinter("C:\\Jat\\jat\\traj\\rvtraj_simu_rbar.txt");
//			x.chaser.sendToLinePrinter(lp1);
//			LinePrinter lp2 = new LinePrinter("C:\\Jat\\jat\\traj\\isstraj_simu_rbar.txt");
//			x.iss.sendToLinePrinter(lp2);
			
			counter = counter + 1;
			j = j + 1;
			
		}
				
	}
	

}
