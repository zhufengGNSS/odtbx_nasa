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
 */

package jat.cm.rendezvous;
import jat.cm.*;
import jat.traj.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.timeRef.*;
/**
 * <P>
 * The DeltaV_Generator Class provides an example of how to generate a rendezvous trajectory 
 * using CW_Guidance.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class DeltaV_Generator {
	
	/** Chaser trajectory */
	public Trajectory chaser = new Trajectory();
	
	/** Target trajectory */
	public Trajectory iss = new Trajectory();
	
	private JGM3DragEOM chaserEOM = new JGM3DragEOM(chaser, 51969.0, 104328.0, 454.4, 2.0);
	
	private JGM3DragEOM issEOM = new JGM3DragEOM(iss, 51969.0, 128990.0, 640.7, 2.35);
	
	private String title;
	
	private String chaserFile;
	
	private String issFile;
		
	private RungeKutta8 rk8 = new RungeKutta8(1.0);
	
	/** Guidance object */
	public CW_Guidance guid;
	
	/** Delta V List */
	public DeltaV_List dvlist = new DeltaV_List();
	
	private String dvfile;
	
			
//    private HarrisPriester hp = new HarrisPriester(); // earth atmosphere model
//	private static final double t_mjd0 = 51969.0;
//	private double mass = 104328.0;
//	private double area = 454.4;
//	private double cd = 2.0;
	
	private int dv_index = 0;
	
	/**
	 * Contructor
	 * @param t title
	 * @param file1 Chaser trajectory file
	 * @param file2 Target trajectory file
	 * @param file3 output file
	 */
	public DeltaV_Generator(String t, String file1, String file2, String file3) {
		this.chaserFile = file1;
		this.issFile = file2;
		this.title = t;
		this.dvfile = file3;
	}
	
				
//	private double currentMJD (double t){
//		double out = this.t_mjd0 + t/86400.0;
//		return out;
//	}	
		
	
	/**
	 * Generate the rendezvous trajectory
	 * @param t0 initial time
	 * @param tf final time
	 * @param x0 initial chaser state vector
	 * @param xref0 initial target state vector
	 */
	public void generate (double t0, double tf, VectorN x0, VectorN xref0) {
		
		// propagate the orbit
		int neqns = x0.length;

        double dt = 1.0;
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
                
        // main integration loop

        while (t < tf) {

        	VectorN dv = new VectorN(3);
        	
        	// check next burn time
        	if (this.dv_index < guid.burntimes.length) {
	        	double tburn = guid.burntimes[this.dv_index];
	        	if (t == tburn) {
	        		VectorN r = new VectorN(oldchaser[0], oldchaser[1], oldchaser[2]);
	        		VectorN v = new VectorN(oldchaser[3], oldchaser[4], oldchaser[5]);        		
	        		Quaternion q = new Quaternion(oldchaser[6], oldchaser[7], oldchaser[8], oldchaser[9]);
	        		VectorN rref = new VectorN(oldiss[0], oldiss[1], oldiss[2]);
	        		VectorN vref = new VectorN(oldiss[3], oldiss[4], oldiss[5]);

	        		dv = guid.computeBurn(t, r, v, rref, vref, this.dv_index);
//	        		Matrix rsw2eci = q.quat2DCM();
	        		RSW_Frame rsw = new RSW_Frame(r, v);
	        		Matrix rsw2eci = rsw.ECI2RSW().transpose();
	        		VectorN dv_eci = rsw2eci.times(dv);
	        		v = v.plus(dv_eci);
	        		
	        		DeltaV deltav = new DeltaV(t, dv);
	        		this.dvlist.add(deltav);
	        		
	        		System.out.println(t+" implementing burn: "+dv);
	        		for (int i = 0; i < 3; i++){
	        			oldchaser[i+3] = v.x[i];
	        		}
	        		this.dv_index = this.dv_index + 1;
	        	}
        	}
        	        	
            newchaser = rk8.step(t, oldchaser, chaserEOM);
            newiss = rk8.step(t, oldiss, issEOM);
            
                        
            for (int i = 0; i < neqns; i++) {
	            oldchaser[i] = newchaser[i];
	            oldiss[i] = newiss[i];
            }
                       
            t = t + dt;

//            if (MathUtils.Frac(t) == 0.0){
		        chaserEOM.print(t, oldchaser);
		        issEOM.print(t, oldiss);
//            	System.out.println(t+"\t"+oldstate[0]+"\t"+oldstate[1]+"\t"+oldstate[2]);
//            }
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
		System.out.println("trajectory serialized");
		dvlist.serialize(this.dvfile);

	}
	

	/** Initialize the chaser.
	 */
	public VectorN initChaser () {

		// create a TwoBody orbit using orbit elements
//		TwoBody sat =
//			new TwoBody(Constants.GM_Earth, 6765500.0, 0.001, 51.8, 0.0, 0.0, 359.8729677);
		double d = 15.0/6765.5000;
		double theta = 360.0 - d * Constants.rad2deg;
		
		TwoBody sat =
			new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, theta);
//		Matrix Cb2i = sat.RSW2ECI();
		
		
		VectorN r = sat.getR();
		VectorN v = sat.getV();
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
		VectorN omega = rsw.omega();
		
		// place it 15 km behind
//		VectorN dr = new VectorN(0.0, -15000.0, 0.0);
//		VectorN dr_eci = Cb2i.times(dr);
//		VectorN dv_eci = omega.crossProduct(dr);
//		r = r.plus(dr_eci);
//		v = v.plus(dv_eci);
		
		VectorN rv = new VectorN(r, v);
				
		Quaternion q0 = new Quaternion(Cb2i);
		VectorN x0 = new VectorN(rv, q0);
		VectorN zero = new VectorN(6);
		x0 = x0.append(zero);
		
		return x0;
	}
	
	/**
	 * Initialize the target
	 */
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
		
		this.guid = new CW_Guidance(sat, 1000.0);
		
		return x0;
	}

	
	public static void main(String[] args){
		
		// set up the generator
		String title = "STS/ISS Rendezvous Trajectory";
		String file1 = "C:\\Jat\\jat\\traj\\rvtraj_vbar.jat";
		String file2 = "C:\\Jat\\jat\\traj\\isstraj_vbar.jat";
		String file3 = "C:\\Jat\\jat\\output\\vbar_delta_v.jat";
		
		DeltaV_Generator x = new DeltaV_Generator(title, file1, file2, file3);
		VectorN x0 = x.initChaser();
		VectorN xref0 = x.initISS();
		
		// set up guidance
		VectorN rtgt = new VectorN(3);
		
		// vbar case
		VectorN gs0 = new VectorN(0.0, -183.0, 0.0);
		VectorN mc4 = new VectorN(3);
		mc4.set(0, 0.0);
		mc4.set(1, -200.0);

		//rbar case
//		VectorN gs0 = new VectorN(-183.0, 0.0, 0.0);
//		VectorN mc4 = new VectorN(-549.0, -274.0, 0.0);
		
		
		double tof = 4620.0;
		double tof1 = 780.0;
		double tof2 = 3000.0;
		double rhoDot0 = -0.2;
		double rhoDotT = 0.0;
		x.guid.computeTargets(mc4, tof, gs0, tof1, rtgt, rhoDot0, rhoDotT, tof2, 4);
		x.guid.tgtList.sendToFile("C:\\Jat\\jat\\output\\vbar_targets.txt");		
		
		double t0 = 0.0;
//		double tf = 5579.0;
		double tf = 9400.0;
		x.generate(t0, tf, x0, xref0);
		System.out.println("Done Generating");
		
        DeltaV_List dvl = DeltaV_List.recover("C:\\Jat\\jat\\output\\vbar_delta_v.jat");
        System.out.println("Printing Recovered DeltaV's");
        int index = 0;
        for(int i = 0; i < dvl.size(); i++) {
	        double t = dvl.time(index);
	        VectorN delv = dvl.dv(index);
	        System.out.println("deltav: "+index+" "+t+" "+delv);
	        index = index + 1;
        }
		
	}
}
