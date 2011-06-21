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
 * 
 */

package jat.alg.estimators;

import jat.alg.integrators.LinePrinter;
import jat.alg.integrators.RungeKutta8;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.sim.*;
import jat.timeRef.RSW_Frame;
import jat.traj.CentralBody;
import jat.traj.CoordinateSystem;
import jat.traj.DistanceUnits;
import jat.traj.TimeUnits;
import jat.traj.Trajectory;
import jat.util.FileUtil;

import java.util.HashMap;
import jat.alg.estimators.*;


/**
* The J2DragProcess.java Class provides the dynamics model for an orbit 
* including J2 and effects of exponential drag.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class JGM4x4SRPProcess15state implements ProcessModel {


	
	//Construct the required classes
	private VectorN xref;
	private Matrix phi;
	private RungeKutta8 rk8 = new RungeKutta8(1.0);
	private JGM4x4SRPEOM15state eom;
	private LinePrinter lp1;
	private LinePrinter lp2;
	public Trajectory traj = new Trajectory();
	HashMap hm;
	public int n;
	private Matrix Q;
	
	/*This file is not generic and has to be modified for any
	 * change in the state.
	 */
	public JGM4x4SRPProcess15state(LinePrinter lp, LinePrinter lp_2, HashMap hm){
		this.hm = hm;
		this.eom = new JGM4x4SRPEOM15state(hm);
		this.lp1 = lp;
		this.lp2 = lp_2;
	    traj.setTitle("Test Trajectory 1");
	    traj.setCentralBody(CentralBody.EARTH);
	    traj.setCoordinateSystem(CoordinateSystem.INERTIAL);
	    traj.setDistanceUnits(DistanceUnits.METERS);
	    traj.setEpoch(1998, 6, 21, 00, 00, 0.0);
	    traj.setTimeUnits(TimeUnits.SECONDS);
	    String[] labels = {"t","x","y","z","xdot","ydot","zdot"};
	    traj.setLabels(labels);
	    
		String fs, dir_in;
        fs = FileUtil.file_separator();
        try{
            dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"input"+fs;
        }catch(Exception e){
            dir_in = "";
        }
		//hm = initializer.parse_file(dir_in+"initialConditions.txt");
	    n = initializer.parseInt(hm,"FILTER.states");
	    Q = parse_Q();
		xref = new VectorN(n);
		phi = new Matrix(n);
	    
	}
  
	/** 
	 * Returns the initial covariance matrix.
	 * @return the initial covariance matrix.
	 */
	public Matrix P0() {

		// initialize sigmas
		double[] sigmas = new double[n];
		
		for(int i = 0; i < initializer.parseInt(hm,"prop.NumSpacecraft"); i++)
		{
			String tmp = "P0."+i+".X";
			sigmas[6*i + 0] = initializer.parseDouble(hm,tmp);
			tmp = "P0."+i+".Y";
			sigmas[6*i + 1] = initializer.parseDouble(hm,tmp);
			tmp = "P0."+i+".Z";
			sigmas[6*i + 2] = initializer.parseDouble(hm,tmp);
			tmp = "P0."+i+".VX";
			sigmas[6*i + 3] = initializer.parseDouble(hm,tmp);
			tmp = "P0."+i+".VY";
			sigmas[6*i + 4] = initializer.parseDouble(hm,tmp);
			tmp = "P0."+i+".VZ";
			sigmas[6*i + 5] = initializer.parseDouble(hm,tmp);
		}
				
		sigmas[12] = initializer.parseDouble(hm,"P0.0.Cr");
		sigmas[13] = initializer.parseDouble(hm,"P0.1.Cr");
		sigmas[14] = initializer.parseDouble(hm,"P0.0.clockBias");

		
		// square the sigmas
		for (int j = 0; j < n; j++) {
			sigmas[j] = sigmas[j] * sigmas[j];
		}

		// initialize covariance
		VectorN sig = new VectorN(sigmas);
		Matrix p = new Matrix(sig);
		return p;
	}
	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	private Matrix parse_Q() {
		
		Matrix Q = new Matrix(n,n);

		for(int i = 0; i < initializer.parseInt(hm,"prop.NumSpacecraft"); i++)
		{
			String tmp = "Q."+i+".X";
			Q.set((6*i + 0),(6*i + 0), initializer.parseDouble(hm,tmp));
			tmp = "Q."+i+".Y";
			Q.set((6*i + 1),(6*i + 1), initializer.parseDouble(hm,tmp));
			tmp = "Q."+i+".Z";
			Q.set((6*i + 2),(6*i + 2), initializer.parseDouble(hm,tmp));
			tmp = "Q."+i+".VX";
			Q.set((6*i + 3),(6*i + 3), initializer.parseDouble(hm,tmp));
			tmp = "Q."+i+".VY";
			Q.set((6*i + 4),(6*i + 4), initializer.parseDouble(hm,tmp));
			tmp = "Q."+i+".VZ";
			Q.set((6*i + 5),(6*i + 5), initializer.parseDouble(hm,tmp));
		}
				
		Q.set(12,12,initializer.parseDouble(hm,"Q.0.Cr"));
		Q.set(13,13,initializer.parseDouble(hm,"Q.1.Cr"));
		Q.set(14,14, initializer.parseDouble(hm,"Q.0.clockBias"));
		return Q;
	}
	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	public Matrix Q(double t, double dt, EstSTM x) {
		
		return this.Q;
	}

	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	public Matrix QRIC(VectorN rECI, VectorN vECI, double t) {
		Matrix Q = new Matrix(9,9);

		Q.set(0,0,1e-10);
		Q.set(1,1,1e-10);
		Q.set(2,2,1e-10);
		Q.set(3,3,1e-12);
		Q.set(4,4,1e-12);
		Q.set(5,5,1e-12);
		Q.set(6,6,1e2);
		Q.set(7,7,1e-8);
		Q.set(8,8,1e-8);
		
		
		
			RSW_Frame rsw = new RSW_Frame(rECI, vECI);
			Matrix M = rsw.ECI2RIC(rECI, vECI);
			Matrix MT = M.transpose();
			
			Matrix Qtmp = new Matrix(3,3);
			Qtmp.set(0,0,Q.get(0,0));
			Qtmp.set(1,1,Q.get(1,1));
			Qtmp.set(2,2,Q.get(2,2));
			
			Matrix Qout = new Matrix(3,3);
			Matrix Qtmp2 = new Matrix(3,3);
			Qtmp2 = MT.times(Qtmp);
			Qout = Qtmp2.times(M);
			
			Q.setMatrix(0,2,0,2,Qout);		
			
			Qtmp.set(0,0,Q.get(3,3));
			Qtmp.set(1,1,Q.get(4,4));
			Qtmp.set(2,2,Q.get(5,5));
			
			Qtmp2 = MT.times(Qtmp);
			Qout = Qtmp2.times(M);

			Q.setMatrix(3,5,3,5,Qout);
			
		return Q;
	}
	
	
	/** 
	 * Returns the initial reference state.
	 * @return the initial reference state.
	 */
	public VectorN xref0() {
		VectorN out = new VectorN(n);
		
		//Note:  Need to decide where the best place to "randomly" insert
		//initial errors should be.  This could be a good place
		//double tmp = Math.random();
		
		//first set in the Satellite States
		for (int i = 0; i < initializer.parseInt(hm, "prop.NumSpacecraft"); i++) {
			String ref = "REF_STATE.";
			String tmp = ref+i +".X";
			out.x[6 * i + 0] = initializer.parseDouble(hm, tmp);
			tmp = ref+ i + ".Y";
			out.x[6 * i + 1] = initializer.parseDouble(hm, tmp);
			tmp = ref+ i + ".Z";
			out.x[6 * i + 2] = initializer.parseDouble(hm, tmp);
			tmp = ref+ i + ".VX";
			out.x[6 * i + 3] = initializer.parseDouble(hm, tmp);
			tmp = ref + i + ".VY";
			out.x[6 * i + 4] = initializer.parseDouble(hm, tmp);
			tmp = ref+ i + ".VZ";
			out.x[6 * i + 5] = initializer.parseDouble(hm, tmp);
			
			// Other clock states can be initialized here . .
			ref = "jat.";
			tmp = ref+"0.Cr";
			out.x[12] = initializer.parseDouble(hm, tmp);
			tmp = ref+"1.Cr";
			out.x[13] = initializer.parseDouble(hm, tmp);
			tmp = ref+"0.clockBias";
			out.x[14] = initializer.parseDouble(hm, tmp);
		}
		
		return out;
	}

	/**
	 * Returns the number of states.
	 * @return the number of states.
	 */
	public int numberOfStates() {
		return this.n;
	}

	/** 
	 * Propagate the state and state transition matrix to the next measurement time.
	 * @param t0 previous time
	 * @param xin array containing state and state transition matrix at previous time.
	 * @param tf next time
	 */
	public double[] propagate(double t0, double[] x, double tf) {
		rk8.setStepSize(tf-t0);
		double[] out = rk8.step(t0, x, eom);
		return out;
	}

	/**
	 * Returns the state transition matrix.
	 * @return the state transition matrix (after propagation).
	 */
	public Matrix phi() {
		return this.phi;
	}

	/**
	 * Returns the current reference state (after propagation).
	 * @return the current reference state.
	 */
	public VectorN xref() {
		return this.xref;
	}

	public void print(double t, VectorN state, Matrix cov) {
		VectorN sigmas = cov.diagonal();
		sigmas = sigmas.ebeSqrt();
		VectorN printvector = new VectorN(state, sigmas);
		VectorN temp = state.get(0,6);
		lp1.println(printvector.toString());
		traj.add(t,temp.x);
		
	}
	
	public void printResiduals(double t, double r1, double r2) {
		double[] y = new double[3];
		y[0] = t;
		y[1] = r1;
		y[2] = r2;
		lp2.print(y);
		System.out.println(t);
	} 
	
	
	public void closeLinePrinter(){
		lp1.close();
        // Serialize the Ouptput
        System.out.println("Searializing the Output . . .");
        traj.serialize("C:\\GOESR\\output\\esttraj.jat");
        System.out.println("trajectory serialized");
	}

}
