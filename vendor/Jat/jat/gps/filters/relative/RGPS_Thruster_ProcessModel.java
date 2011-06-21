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

package jat.gps.filters.relative;

import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.gps.*;
import jat.gps.filters.*;
import jat.cm.*;
//import jat.math.*;
//import jat.gps_ins.*;
import jat.timeRef.*;
//import jat.forces.*;

/**
* The RGPS_ProcessModel.java Class provides the process model
* for the RGPS-only EKF.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RGPS_Thruster_ProcessModel implements ProcessModel {
	
	private IonoModel iono = new IonoModel();
	
	private ReceiverFilterModel rcvr = new ReceiverFilterModel();
	
//	private JGM3 jgm3 = new JGM3(12,12);
			
	private RungeKutta8 rk8;
	private double dt = 1.0;
	
	private GPS_Constellation gpscon;
	
	private URE_Model ure;
	
	private int issIndex;
	private int iaIndex;
	
	private int nsv;
	
	private LinePrinter lp1;
	private LinePrinter lp2;
	private LinePrinter lp3;

	private RGPS_EOM coastEOM;
	private RGPS_Thruster_EOM burnEOM;
	
	public FiniteBurnList burnlist;
	private int dv_index = 0;
    private boolean burnstarted = false;
	
	
	/**
	 * Constructor
	 * @param gps GPS_Constellation
	 * @param l1 LinePrinter for state, covariance output
	 * @param l2 LinePrinter for relative state, covariance output
	 * @param l3 LinePrinter for residual output
	 * @param burnfile String containing the finite burns filename
	 */
	public RGPS_Thruster_ProcessModel(GPS_Constellation gps, LinePrinter l1, LinePrinter l2, LinePrinter l3, String burnfile) {
		this.gpscon = gps;
		this.nsv = gps.size();
		this.issIndex = 10 + this.nsv;
		this.iaIndex = 19 + this.nsv;
		
		this.burnlist = FiniteBurnList.recover(burnfile);
		

		this.ure = new URE_Model(this.nsv);
		this.lp1 = l1;
		this.lp2 = l2;
		this.lp3 = l3;
		this.rk8 = new RungeKutta8(this.dt);
		this.coastEOM = new RGPS_EOM(this.nsv, this.issIndex, this.numberOfStates(), this.iono, this.rcvr, this.ure);
		this.burnEOM = new RGPS_Thruster_EOM(this.nsv, this.issIndex, this.numberOfStates(), this.iono, this.rcvr, this.ure);
		
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
		TwoBody orbit2 = new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, 0.0);
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
		
		// set ISS initial states
		VectorN rv = orbit2.rv;
		
		// perturb ICs
//		rv = rv.plus(pt);
		out.set(this.issIndex, rv);
		
		// set initial ISS clock states, integer amb states to zero
		
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
			out.set((this.iaIndex+i),(this.iaIndex+i), 8.0);
		}
		
		int k = this.issIndex;
		out.set(k, k, sigma2_r);
		out.set(k+1, k+1, sigma2_r);
		out.set(k+2, k+2, sigma2_r);
		out.set(k+3, k+3, sigma2_v);
		out.set(k+4, k+4, sigma2_v);
		out.set(k+5, k+5, sigma2_v);
		out.set(k+6, k+6, sigma2_bc);
		out.set(k+7, k+7, sigma2_dc);
		out.set(k+8, k+8, sigma2_drag);
		
		
		return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#numberOfStates()
	 */
	public int numberOfStates() {
		int n = 19 + 2*this.nsv;
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
//	public Matrix Q(double t, double dt, VectorN x) {
//		int n = this.numberOfStates();
//		Matrix q = new Matrix(n, n);
//		
//		
//		// adaptive tuning
//		double tsw = 7800.0;
//		double sp = 0.0;		
//		if (t < tsw) {
//			sp = 1.0E-06;
//		}
//		else {
//			sp = 1.0E-12;
//		}
//
//		// nominal tuning
////		double sp = 1.0E-06;
//
//		// blockage tuning
////		double sp = 1.0E-12;
//
//		// common SV tuning
////		double sp = 1.0E-11;
//
//		double sp3 = sp/3.0;
//		double sp2 = sp/2.0;
//		q.set(0, 0, sp3);
//		q.set(1, 1, sp3);
//		q.set(2, 2, sp3);
//		q.set(3, 3, sp);
//		q.set(4, 4, sp);
//		q.set(5, 5, sp);
//		q.set(0, 3, sp2);
//		q.set(1, 4, sp2);
//		q.set(2, 5, sp2);
//		q.set(3, 0, sp2);
//		q.set(4, 1, sp2);
//		q.set(5, 2, sp2);
//
////		q.setMatrix(6, 6, rcvr.biasQ(dt).times(10.0));
////		q.setMatrix(6, 6, rcvr.biasQ(dt));
//		
//
//		q.set(8, 8, DragProcessModel.dragQ(dt));
//
////		q.set(9, 9, 10.0*iono.ionoQ(dt));
////		q.set(9, 9, iono.ionoQ(dt));
//		
//		// adaptive
//		if (t < tsw) {
//			q.setMatrix(6, 6, rcvr.biasQ(dt));
//			q.set(9, 9, iono.ionoQ(dt));
//		} else {
//			q.setMatrix(6, 6, rcvr.biasQ(dt).times(10.0));
//			q.set(9, 9, 10.0*iono.ionoQ(dt));
//		}
//		
//		for (int i = 0; i < this.nsv; i++) {
//			
//			// adaptive
//			if (t < tsw) {
//				q.set((i+10),(i+10), ure.biasQ());
//			} else {
//				q.set((i+10),(i+10), 10.0*ure.biasQ());
//			}
//
//			q.set((i+this.iaIndex),(i+this.iaIndex), 1.0E-6*dt);
//		}
//		
//		int kk = this.issIndex;
//
//		// adaptive tuning		
//		if (t < 7200.0) {		
//			sp = 1.0E-05;
//			sp3 = sp/3.0;
//			sp2 = sp/2.0;
//		}
//		
//		q.set(kk, kk, sp3);
//		q.set(kk+1, kk+1, sp3);
//		q.set(kk+2, kk+2, sp3);
//		q.set(kk+3, kk+3, sp);
//		q.set(kk+4, kk+4, sp);
//		q.set(kk+5, kk+5, sp);
//		q.set(kk, kk+3, sp2);
//		q.set(kk+1, kk+4, sp2);
//		q.set(kk+2, kk+5, sp2);
//		q.set(kk+3, kk, sp2);
//		q.set(kk+4, kk+1, sp2);
//		q.set(kk+5, kk+2, sp2);
//		q.setMatrix(kk+6, kk+6, rcvr.biasQ(dt));
//		q.set(kk+8, kk+8, DragProcessModel.dragQ(dt));
//				
//		return q;
//	}

	public Matrix Q(double t, double dt, EstSTM stm) {
		int n = this.numberOfStates();
		Matrix q = new Matrix(n, n);

		VectorN x = stm.state();
		Matrix phi = stm.phi();
		Matrix phiT = phi.transpose();
		
		
		// adaptive tuning
//		double tsw = 7200.0;
//		double sp = 0.0;		
//		if (t < tsw) {
//			sp = 1.0E-06;
//		}
//		else {
//			sp = 1.0E-09;
//		}

		double sp = 1.0E-09;


		q.set(3, 3, sp);
		q.set(4, 4, sp);
		q.set(5, 5, sp);
		
		q.set(8, 8, DragProcessModel.Q());
		
		// adaptive
//		if (t < tsw) {
//			q.setMatrix(6, 6, rcvr.Q());
//			q.set(9, 9, iono.Q());
//		} else {
			q.setMatrix(6, 6, rcvr.Q().times(10.0));
			q.set(9, 9, 10.0*iono.Q());
//		}
		
		for (int i = 0; i < this.nsv; i++) {
			
			// adaptive
//			if (t < tsw) {
//				q.set((i+10),(i+10), ure.Q());
//			} else {
				q.set((i+10),(i+10), 10.0*ure.Q());
//			}

			q.set((i+this.iaIndex),(i+this.iaIndex), 1.0E-6);
		}
		
		int kk = this.issIndex;

		q.set(kk+3, kk+3, sp);
		q.set(kk+4, kk+4, sp);
		q.set(kk+5, kk+5, sp);
		q.setMatrix(kk+6, kk+6, rcvr.Q());
		q.set(kk+8, kk+8, DragProcessModel.Q());
		
		// Relative Nav correlation
		double rho = 0.99;
		Matrix qab = new Matrix(6, 6);
		for (int i = 0; i < 6; i++) {
			double qa = q.get(i, i);
			double qb = q.get(kk+i, kk+i);
			double val = rho * Math.sqrt(qa * qb);
			q.set(i, (kk+i), val);
			q.set((kk+i), i, val);
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
	
	
	public void print(double t, VectorN state, Matrix cov) {
		
		// absolute state processing
		VectorN sigmas = cov.diagonal();
		sigmas = sigmas.ebeSqrt();
		VectorN printvector = new VectorN(state, sigmas);
		lp1.print(t, printvector.x);
		
		// relative state processing
		VectorN r_chaser = state.get(0,3);
		VectorN v_chaser = state.get(3,3);
		VectorN r_iss = state.get(this.issIndex,3);
		VectorN v_iss = state.get(this.issIndex+3,3);
		RSW_Frame rsw = new RSW_Frame(r_iss, v_iss);
		VectorN r_rel = r_chaser.minus(r_iss);
		VectorN v_rel = v_chaser.minus(v_iss);
		VectorN rv_rel = rsw.transform(r_rel, v_rel);
		
		Matrix p_chaser = cov.getMatrix(0, 5, 0, 5);
		Matrix p_iss = cov.getMatrix(this.issIndex, this.issIndex+5,this.issIndex, this.issIndex+5);
		Matrix p_a = cov.getMatrix(0, 5, this.issIndex, this.issIndex+5); 
		Matrix p_b = cov.getMatrix(this.issIndex, this.issIndex+5, 0, 5);
		Matrix temp1 = p_chaser.plus(p_iss);
		Matrix temp2 = p_a.plus(p_b);
		Matrix p_rel = temp1.minus(temp2);
		
		Matrix eci2rsw = rsw.ECI2RSW();
		
		Matrix T = new Matrix(6,6);
		T.setMatrix(0, 0, eci2rsw);
		T.setMatrix(3, 3, eci2rsw);
		Matrix Ttrans = T.transpose();
		
		Matrix temp3 = p_rel.times(Ttrans);
		Matrix p_rsw = T.times(temp3);
		sigmas = p_rsw.diagonal();
		sigmas = sigmas.ebeSqrt();
		printvector = new VectorN(rv_rel, sigmas);
		lp2.print(t, printvector.x);

	}
	
	public void printResiduals(double t, double r1, double r2) {
		double[] y = new double[3];
		y[0] = t;
		y[1] = r1;
		y[2] = r2;
		lp3.print(y);
	} 
	

	public void closeLinePrinter(){
		lp1.close();
		lp2.close();
		lp3.close();
	}
				

}
