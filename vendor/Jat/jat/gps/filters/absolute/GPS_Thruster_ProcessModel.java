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
 * File Created on May 19, 2003
 */

package jat.gps.filters.absolute;

import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.gps.*;
import jat.gps.filters.*;
import jat.cm.*;
//import jat.math.*;
//import jat.gps_ins.*;
//import jat.timeRef.*;
//import jat.forces.*;

/**
* The GPS_ProcessModel.java Class provides the process model
* for the GPS-only EKF.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class GPS_Thruster_ProcessModel implements ProcessModel {
	
	private IonoModel iono = new IonoModel();
	
	private ReceiverFilterModel rcvr = new ReceiverFilterModel();
	
			
	private RungeKutta8 rk8;
	private double dt = 1.0;
	
	private GPS_Constellation gpscon;
	
	private URE_Model ure;
		
	private int nsv;
	
	private LinePrinter lp1;
	private LinePrinter lp2;

	private GPS_EOM coastEOM;
	private GPS_Thruster_EOM burnEOM;
	
	public FiniteBurnList burnlist;
	private int dv_index = 0;
    private boolean burnstarted = false;
	
	/**
	 * Constructor
	 * @param gps GPS_Constellation
	 * @param l1 LinePrinter for state, covariance output
	 * @param l2 LinePrinter for residual output
	 * @param burnfile String containing the finite burns filename
	 */
	public GPS_Thruster_ProcessModel(GPS_Constellation gps, LinePrinter l1, LinePrinter l2, String burnfile) {
		this.gpscon = gps;
		this.nsv = gps.size();
		
		this.burnlist = FiniteBurnList.recover(burnfile);
		
		this.ure = new URE_Model(this.nsv);
		this.lp1 = l1;
		this.lp2 = l2;
		this.rk8 = new RungeKutta8(this.dt);
		this.coastEOM = new GPS_EOM(this.nsv, this.numberOfStates(), this.iono, this.rcvr, this.ure);
		this.burnEOM = new GPS_Thruster_EOM(this.nsv, this.numberOfStates(), this.iono, this.rcvr, this.ure);
		
	}
	
	/**
	 * return the initial state vector
	 */
	public VectorN xref0() {
		int n = this.numberOfStates();
		
		double d = 15.0/6765.5000;
		double theta = 360.0 - d * Constants.rad2deg;
		
		// set initial position, velocity, attitude
		TwoBody orbit1 = new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, theta);
		VectorN r = orbit1.getR();
		VectorN v = orbit1.getV();
		
		// perturb ICs
//		double [] pert = new double[6];
//		pert[0] = 0.0;
//		pert[1] = 0.0;
//		pert[2] = 0.0;
//		pert[3] = 0.0;
//		pert[4] = 0.0;
//		pert[5] = 0.0;
//		VectorN pt = new VectorN(pert);
//		rv = rv.plus(pt);
//		RSW_Frame rsw = new RSW_Frame(r, v);
//		Matrix Cb2i = rsw.ECI2RSW().transpose();

		//perturb ICs
		VectorN mean = new VectorN(3);
		VectorN sigp = new VectorN(3);
		VectorN sigv = new VectorN(3);
		sigp.set(10.0);
		sigv.set(0.1);
		GaussianVector ppert = new GaussianVector(mean, sigp);
		r = r.plus(ppert);
		GaussianVector vpert = new GaussianVector(mean, sigv);
		v = v.plus(vpert);

		
		VectorN out = new VectorN(n);
		out.set(0, r);
		out.set(3, v);
//        Quaternion q0 = new Quaternion(Cb2i);
//        out.set(6, q0);
                
        // set initial clock bias, assume perfect
        out.set(6, 0.0);
        out.set(7, 0.0);
        
        // set initial drag state
        out.set(8, 0.0);

		// set initial iono state
		out.set(9, 0.0);
		
		// set initial ure states to zero
		
		
		return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#P0()
	 */
	public Matrix P0() {
		Matrix out = new Matrix(this.numberOfStates());
		
		double sigma_r = 10.0;
		double sigma_v = 0.1;
		double sigma_bc = 1.0E-06;
		double sigma_dc = 1.0E-06;
		double sigma2_r = sigma_r * sigma_r;
		double sigma2_v = sigma_v * sigma_v;
		double sigma2_bc = sigma_bc * sigma_bc;
		double sigma2_dc = sigma_dc * sigma_dc;
		double sigma2_ure = ure.sigma()*ure.sigma();
		
		double sigma2_iono = iono.sigma * iono.sigma;
//		double sigma2_ia = 1.0E+12*GPS_Utils.lambda*GPS_Utils.lambda;

		double sigma2_drag = 0.8 * 0.8;
		
		out.set(0, 0, sigma2_r);
		out.set(1, 1, sigma2_r);
		out.set(2, 2, sigma2_r);
		out.set(3, 3, sigma2_v);
		out.set(4, 4, sigma2_v);
		out.set(5, 5, sigma2_v);
		out.set(6, 6, sigma2_bc);
		out.set(7, 7, sigma2_dc);
		out.set(8, 8, sigma2_drag);
		out.set(9, 9, sigma2_iono);

		for (int i = 0; i < this.nsv; i++){
			out.set((i+10),(i+10), sigma2_ure);
		}		
		
		return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#numberOfStates()
	 */
	public int numberOfStates() {
		int n = 10 + this.nsv;
		return n;
	}
	
	
	private Matrix abar(VectorN f) {
		f.checkVectorDimensions(3);
		double [] x = f.getArray();
        Matrix out = new Matrix(4, 3);
        out.set(0, 1,  -x[2]);
        out.set(0, 2,  x[1]);
        out.set(1, 0,  x[2]);
        out.set(1, 2,  -x[0]);
        out.set(2, 0,  -x[1]);
        out.set(2, 1,  x[0]);
        out.set(3, 0,  -x[0]);
        out.set(3, 1,  -x[1]);
        out.set(3, 2,  -x[2]);
        return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#Q(double, double)
	 */
	public Matrix Q(double t, double dt, EstSTM stm) {
		int n = this.numberOfStates();

		VectorN x = stm.state();
		Matrix phi = stm.phi();
		Matrix phiT = phi.transpose();


		Matrix q = new Matrix(n, n);
		
		double sp = 1.0E-06;
		q.set(3, 3, sp);
		q.set(4, 4, sp);
		q.set(5, 5, sp);

		q.setMatrix(6, 6, rcvr.Q());

		q.set(8, 8, DragProcessModel.Q());

		q.set(9, 9, iono.Q());
		
		for (int i = 0; i < this.nsv; i++) {
			q.set((i+10),(i+10), ure.Q());
		}

		Matrix out = phi.times(q.times(phiT));
		out = out.times(dt);			
		
				
		return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#propagate(double, double[], double)
	 */
	public double[] propagate(double t0, double[] xin, double tf) {
		
		if ((tf - t0) != this.dt) {
			System.out.println("propagate step size messed up");
		}
		
		double [] xnew = new double[xin.length];

        	// check next burn time
        	FiniteBurn burn = burnlist.get(dv_index);
        	double tstart = burn.tstart;
        	double tstop = burn.tstop;        	        	
        	if ((t0 >= tstart) && (t0 < tstop)) {
        		burnEOM.setBurn(burn);
        		xnew = rk8.step(t0, xin, burnEOM);
//        		xnew = rk8.step(t0, xin, coastEOM);
        		this.burnstarted = true;
        	}
        	else {
        		xnew = rk8.step(t0, xin, coastEOM);
        		if ((t0 > tstop) && (burnstarted)) {
        			System.out.println(t0+" burn "+dv_index+" completed, incrementing");
        			if (burnlist.hasNext(dv_index)) {
        				 dv_index = dv_index + 1;
        			}
        			this.burnstarted = false;
        		}
        	}
				
		
		return xnew;
	}
		
	/**
	 * Print output
	 * @param t sim time in seconds
	 * @param state VectorN containing the state vector
	 * @param cov Matrix containing the covariance matrix
	 */
	public void print(double t, VectorN state, Matrix cov) {
		
		// absolute state processing
		VectorN sigmas = cov.diagonal();
		sigmas = sigmas.ebeSqrt();
		VectorN printvector = new VectorN(state, sigmas);
		lp1.print(t, printvector.x);
		
	}
	
	/**
	 * Print residuals
	 * @param t sim time in seconds
	 * @param r1 residual before measurement update
	 * @param r2 residual after measurement update
	 */	
	public void printResiduals(double t, double r1, double r2) {
		double[] y = new double[3];
		y[0] = t;
		y[1] = r1;
		y[2] = r2;
		lp2.print(y);
	} 
	

	/**
	 * close the LinePrinters
	 */
	public void closeLinePrinter(){
		lp1.close();
	}
				

}
