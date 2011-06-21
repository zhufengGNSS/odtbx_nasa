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
package jat.alg.estimators;

import jat.alg.integrators.RungeKutta8;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.sim.initializer;
import jat.spacetime.RSW_Frame;

import java.util.HashMap;

public class JGM4x4DragProcess9state implements ProcessModel {

//	Construct the required classes
	private VectorN xref;
	private Matrix phi;
	private RungeKutta8 rk8 = new RungeKutta8(1.0);
	private JGM4x4SRPEOM9state eom;// = new JGM4x4SRPEOM9state();
	//public Trajectory traj = new Trajectory();
	HashMap hm;
	public int n;
	private Matrix Q;
	private Matrix QXYZ;
	
	/*This file is not generic and has to be modified for any
	 * change in the state.
	 */
	public JGM4x4DragProcess9state(HashMap hm){
		this.hm = hm;
		this.eom = new JGM4x4SRPEOM9state(hm);
		//traj.setTitle("Test Trajectory 1");
	    //traj.setCentralBody(CentralBody.EARTH);
	    //traj.setCoordinateSystem(CoordinateSystem.INERTIAL);
	    //traj.setDistanceUnits(DistanceUnits.METERS);
	    //traj.setEpoch(1998, 6, 21, 00, 00, 0.0);
	    //traj.setTimeUnits(TimeUnits.SECONDS);
	    //String[] labels = {"t","x","y","z","xdot","ydot","zdot"};
	    //traj.setLabels(labels);
	    
		//hm = initializer.parse_file(dir_in+"initialConditions.txt");
	    n = initializer.parseInt(hm,"FILTER.states");
	    //num_sc = initializer.parseInt(hm,"prop.NumSpacecraft");
	    Q = parse_Q();
	    QXYZ = parse_QRIC();
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
				
		sigmas[6] = initializer.parseDouble(hm,"P0.0.clockBias");
		sigmas[7] = initializer.parseDouble(hm,"P0.0.clockDrift");
		//sigmas[8] = initializer.parseDouble(hm,"P0.0.Cr");

		
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
				
		Q.set(6,6,initializer.parseDouble(hm,"Q.0.clockBias"));
		Q.set(7,7,initializer.parseDouble(hm,"Q.0.clockDrift"));
		//Q.set(8,8, initializer.parseDouble(hm,"Q.0.Cr"));
		
		return Q;
	}
	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	public Matrix Q(double t, double dt, EstSTM x) {
		return Q;
	}
	private Matrix parse_QRIC(){
		Matrix QXYZ = new Matrix(n,n);

		QXYZ.set(0,0,10e-11);
		QXYZ.set(1,1,10e-11);
		QXYZ.set(2,2, 10e-11);
		QXYZ.set(3,3,1e-13);
		QXYZ.set(4,4,1e-13);
		QXYZ.set(5,5,1e-13);
				
		QXYZ.set(6,6,initializer.parseDouble(hm,"Q.0.clockBias"));
		QXYZ.set(7,7,initializer.parseDouble(hm,"Q.0.clockDrift"));
		//QXYZ.set(8,8, initializer.parseDouble(hm,"Q.0.Cr"));
		return QXYZ;
	}
	/**
	 * Returns the process noise matrix.
	 * @param t time
	 * @param dt dt = current time - previous time.
	 * @return the process noise matrix.
	 */
	public Matrix QRIC(VectorN rECI, VectorN vECI, double t) {
		Matrix QRIC = new Matrix(n,n);
//		Matrix QXYZ = new Matrix(n,n);
//
//		QXYZ.set(0,0,10e-11);
//		QXYZ.set(1,1,10e-11);
//		QXYZ.set(2,2, 10e-11);
//		QXYZ.set(3,3,1e-13);
//		QXYZ.set(4,4,1e-13);
//		QXYZ.set(5,5,1e-13);
//				
//		QXYZ.set(6,6,initializer.parseDouble(hm,"Q.0.clockBias"));
//		QXYZ.set(7,7,initializer.parseDouble(hm,"Q.0.clockDrift"));
//		QXYZ.set(8,8, initializer.parseDouble(hm,"Q.0.Cr"));
		
		
		
			//RSW_Frame rsw = new RSW_Frame(rECI, vECI);
			Matrix M = RSW_Frame.ECI2RIC(rECI, vECI);
			Matrix MT = M.transpose();
			
			Matrix Qtmp = new Matrix(3,3);
			Qtmp.set(0,0,QXYZ.get(0,0));
			Qtmp.set(1,1,QXYZ.get(1,1));
			Qtmp.set(2,2,QXYZ.get(2,2));
			
			Matrix Qout = new Matrix(3,3);
			Matrix Qtmp2 = new Matrix(3,3);
			Qtmp2 = MT.times(Qtmp);
			Qout = Qtmp2.times(M);
			
			QRIC.setMatrix(0,2,0,2,Qout);		
			
			Qtmp.set(0,0,QXYZ.get(3,3));
			Qtmp.set(1,1,QXYZ.get(4,4));
			Qtmp.set(2,2,QXYZ.get(5,5));
			
			Qtmp2 = MT.times(Qtmp);
			Qout = Qtmp2.times(M);

			QRIC.setMatrix(3,5,3,5,Qout);
			QRIC.set(6,6,QXYZ.get(6,6));
			QRIC.set(7,7,QXYZ.get(7,7));
			QRIC.set(8,8,QXYZ.get(8,8));
			
		return QRIC;
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
			tmp = ref+"0.clockBias";
			out.x[6] = initializer.parseDouble(hm, tmp);
			tmp = ref+"0.clockDrift";
			out.x[7] = initializer.parseDouble(hm, tmp);
			//tmp = ref+"0.Cr";
			//out.x[8] = initializer.parseDouble(hm, tmp);
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
//		VectorN sigmas = cov.diagonal();
//		sigmas = sigmas.ebeSqrt();
//		VectorN printvector = new VectorN(state, sigmas);
//		VectorN temp = state.get(0,6);
//		lp1.println(printvector.toString());
//		traj.add(t,temp.x);
//		
	}
	
	public void printResiduals(double t, double r1, double r2) {
//		double[] y = new double[3];
//		y[0] = t;
//		y[1] = r1;
//		y[2] = r2;
//		lp2.print(y);
//		System.out.println(t);
	} 
	
	
	public void closeLinePrinter(){
//		lp1.close();
//        // Serialize the Ouptput
//        System.out.println("Searializing the Output . . .");
//        traj.serialize("C:\\GOESR\\output\\esttraj.jat");
//        System.out.println("trajectory serialized");
	}

}
