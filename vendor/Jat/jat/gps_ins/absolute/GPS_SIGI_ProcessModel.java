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

package jat.gps_ins.absolute;

import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.gps.*;
import jat.gps.filters.*;
import jat.cm.*;
import jat.math.*;
import jat.gps_ins.*;
import jat.ins.*;
import jat.timeRef.*;
import jat.forces.*;
import jat.forces.gravity.earth.J2Gravity;

/**
* The GPS_SIGI_ProcessModel.java Class provides the process model for GPS/SIGI
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class GPS_SIGI_ProcessModel implements Derivs, ProcessModel {
	
	private INS_MeasurementList insData;
	
	private SIGIGyroFilterModel gyro = new SIGIGyroFilterModel();
	
	private SIGIAccelFilterModel accel = new SIGIAccelFilterModel();
	
	private IonoModel iono = new IonoModel();
	
	private DragProcessModel dragModel = new DragProcessModel();
	
	private ReceiverFilterModel rcvr = new ReceiverFilterModel();
	
	private RK4_INS rk4;
	private double dt = 1.0;
		
	private GPS_Constellation gpscon;
	
	private URE_Model ure;
	
	private double t_mjd0 = 51969.0;
	
	private int nsv;
	
	private LinePrinter lp1;
	private LinePrinter lp2;
	
	/**
	 * Constructor
	 * @param ins INS_MeasurementList
	 * @param gps GPS_Constellation
	 * @param l1 LinePrinter for state, covariance output
	 * @param l2 LinePrinter for residual output
	 */
	public GPS_SIGI_ProcessModel(INS_MeasurementList ins, GPS_Constellation gps, LinePrinter l1, LinePrinter l2) {
		this.insData = ins;
		this.gpscon = gps;
		this.nsv = gps.size();

		this.ure = new URE_Model(this.nsv);
		this.lp1 = l1;
		this.lp2 = l2;
		this.rk4 = new RK4_INS(this.insData);
	}
	
	/**
	 * return the initial state vector
	 * @return the initial state vector
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


		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
		
		VectorN out = new VectorN(n);
		out.set(0, r);
		out.set(3, v);
        Quaternion q0 = new Quaternion(Cb2i);
        out.set(6, q0);
                
        // set initial gyro bias, assume perfect
        out.set(10, 0.0);
        out.set(11, 0.0);
        out.set(12, 0.0);
        
        // set initial accel bias, assume perfect
        out.set(13, 0.0);
        out.set(14, 0.0);
        out.set(15, 0.0);
        
        // set initial clock bias, assume perfect
        out.set(16, 0.0);
        out.set(17, 0.0);
//        out.set(16, 1.0E-02*con.c);
//        out.set(17, 6.7E-07*con.c);

		// set initial iono state
		out.set(18, 0.0);
		
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
		double sigma_q = 0.001;
		double sigma_bg = 0.001 * MathUtils.DEG2RAD / 3600.0;
		double sigma_ba = 1.6E-04 * 9.81;
		double sigma_bc = 1.0E-06;
		double sigma_dc = 1.0E-06;
		double sigma2_r = sigma_r * sigma_r;
		double sigma2_v = sigma_v * sigma_v;
		double sigma2_q = sigma_q * sigma_q;
		double sigma2_bg = sigma_bg * sigma_bg;
		double sigma2_ba = sigma_ba * sigma_ba;
		double sigma2_bc = sigma_bc * sigma_bc;
		double sigma2_dc = sigma_dc * sigma_dc;
		double sigma2_ure = ure.sigma()*ure.sigma();
		
		double sigma2_iono = iono.sigma * iono.sigma;
		
		out.set(0, 0, sigma2_r);
		out.set(1, 1, sigma2_r);
		out.set(2, 2, sigma2_r);
		out.set(3, 3, sigma2_v);
		out.set(4, 4, sigma2_v);
		out.set(5, 5, sigma2_v);
		out.set(6, 6, sigma2_q);
		out.set(7, 7, sigma2_q);
		out.set(8, 8, sigma2_q);
		out.set(9, 9, sigma2_q);
		out.set(10, 10, sigma2_bg);
		out.set(11, 11, sigma2_bg);
		out.set(12, 12, sigma2_bg);
		out.set(13, 13, sigma2_ba);
		out.set(14, 14, sigma2_ba);
		out.set(15, 15, sigma2_ba);
		out.set(16, 16, sigma2_bc);
		out.set(17, 17, sigma2_dc);
		out.set(18, 18, sigma2_iono);

		for (int i = 0; i < this.nsv; i++){
			out.set((i+19),(i+19), sigma2_ure);
		}
				
		
		return out;
	}

	/**
	 * @see jat.alg.estimators.ProcessModel#numberOfStates()
	 */
	public int numberOfStates() {
		int n = 19 + this.nsv;
		return n;
	}
	
	/**
	 * @see jat.alg.integrators.Derivatives#derivs(double, double[])
	 */
	public double[] derivs(double t, double[] x, INS_Measurement measl_1, INS_Measurement measl, int sw) {
		
		int n = this.numberOfStates();
		
		VectorN out = new VectorN(x.length);

		// get the gyro measurements
		VectorN dthetal_1 = measl_1.omega;
		VectorN dthetal = measl.omega;
		
		// get accelerometer measurements
		VectorN dvl_1 = measl_1.f;
		VectorN dvl = measl.f;
		
		// strip out the incoming data
		VectorN r = new VectorN(x[0], x[1], x[2]);
		VectorN v = new VectorN(x[3], x[4], x[5]);
		Quaternion q = new Quaternion(x[6], x[7], x[8], x[9]);
		
		q.unitize();
		
		EstSTM stm = new EstSTM(x, n);
		Matrix phi = stm.phi();
		
		VectorN bg = new VectorN(x[10], x[11], x[12]);
		VectorN ba = new VectorN(x[13], x[14], x[15]);
		
		// strip off clock states		
		VectorN clock1 = new VectorN(2);
		clock1.set(0, x[16]);
		clock1.set(1, x[17]);

		double del_iono = x[18];  // iono state
		
		// strip off incoming ure states
		VectorN urevec = new VectorN(this.nsv);
		for (int i = 0; i < this.nsv; i++) {
			urevec.set(i, x[i+19]);
		}
								
		// Body to inertial DCM
		Matrix cib = q.quat2DCM();
		
		
		// form the appropriate specific force and omega vectors
		VectorN f;
		VectorN rate;

		switch(sw) {
			case 0: 
				f = this.dvl_2(dvl_1, dvl);
				rate = this.wl_2(dthetal_1, dthetal);	 
				break;
			case 1: 
				f = this.dvl_1(dvl_1, dvl); 
				rate = this.wl_1(dthetal_1, dthetal);	 
				break;
			case 2: 
				f = this.dvl_0(dvl_1, dvl); 
				rate = this.wl_0(dthetal_1, dthetal);	 
				break;
			default: 
				f = new VectorN(3);
				rate = new VectorN(3);
				break;
		}
				
		// compensate for IMU errors
		f = f.plus(ba);
		rate = rate.plus(bg);
				
		// position derivatives
		out.set(0, v);
		
		// velocity derivatives
//		TwoBody orbit = new TwoBody(Constants.GM_Earth, r, v);
//		VectorN g = orbit.local_grav();	
		J2Gravity j2chaser = new J2Gravity(r);
		VectorN g = j2chaser.local_gravity();	

		double Mjd = this.t_mjd0 + t/86400.0;
        EarthRef ref = new EarthRef(Mjd);
//        ref.setIERS(3.3E-07, 1.17E-06, 0.649232);
//        Matrix E = ref.eci2ecef();
//
//        // Acceleration due to harmonic gravity field        
//        VectorN g = jgm3.gravity(r, E);


		VectorN sf = cib.times(f);
		VectorN vdot = sf.plus(g);
		out.set(3, vdot);
		
		// quaternion derivatives
		Matrix omega = Quaternion.omega(rate);
		VectorN qdot = omega.times(q);
		out.set(6, qdot);
		
		// gyro bias derivatives
		VectorN bgdot = this.gyro.biasProcess(bg);
		out.set(10, bgdot);
		
		// accelerometer bias derivatives
		VectorN badot = this.accel.biasProcess(ba);
		out.set(13, badot);
		
		// GPS clock model derivatives
		VectorN bcdot = rcvr.biasProcess(clock1);
		out.set(16, bcdot);
		
		// iono derivs
		double ionodot = iono.ionoProcess(del_iono);
		out.set(18, ionodot);
		
		//ure derivs
		VectorN uredot = ure.ureProcess(urevec);
		out.set(19, uredot);


		// A matrix
		Matrix A = new Matrix(n, n);
		
		// position rows
		Matrix eye = new Matrix(3);
		A.setMatrix(0, 3, eye);
		
		// velocity rows
		Matrix G = j2chaser.gravityGradient();
		A.setMatrix(3, 0, G);
		
		Matrix abar = this.abar(f);
        Matrix atrans = abar.transpose();
        Matrix rhat = q.rMatrix();
        Matrix rtrans = rhat.transpose();
        Matrix d = cib.times(atrans.times(rtrans));
        d = d.times(2.0);
        A.setMatrix(3, 6, d);
        A.setMatrix(3, 13, cib);
        
        // quaternion rows
        A.setMatrix(6, 6, omega);
        Matrix qmat = q.qMatrix();
        A.setMatrix(6, 10, qmat);
        
        // gyro bias rows
        Matrix taug = eye.times(-1.0/SIGIGyroFilterModel.correlationTime);
        A.setMatrix(10, 10, taug);
		
        // accel bias rows
        Matrix taua = eye.times(-1.0/SIGIAccelFilterModel.correlationTime);
        A.setMatrix(13, 13, taua);
        
        //clock drift row
        A.set(16, 17, 1.0);
        
        // iono row
        double tau_iono = -1.0/iono.correlationTime;
        A.set(18, 18, tau_iono);
        
        // ure part
        Matrix bigeye = new Matrix(this.nsv);
        Matrix tau_ure = eye.times(-1.0/URE_Model.correlationTime);
        A.setMatrix(19, 19, tau_ure);
        
        		
		// phi derivatives
		Matrix phidot = A.times(phi);
		
		// put phi derivatives into output array
		int k = n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				out.x[k] = phidot.A[i][j];
				k = k + 1;
			}
		}
		
		return out.x;
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
////		double sp = 1.0E-03;
//
////		double sp = 1.0E-14;
//
//
//
////		double sp = 1.0E-08;
//
//		double sp = 1.0E-05;
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
//		q.set(6, 6, 1.0E-25);
//		q.set(7, 7, 1.0E-25);
//		q.set(8, 8, 1.0E-25);
//		q.set(9, 9, 1.0E-25);
//		
//		q.set(10, 10, gyro.biasQ());
//		q.set(11, 11, gyro.biasQ());		
//		q.set(12, 12, gyro.biasQ());		
//		q.set(13, 13, accel.biasQ());
//		q.set(14, 14, accel.biasQ());		
//		q.set(15, 15, accel.biasQ());
//
//		q.setMatrix(16, 16, rcvr.biasQ(dt));
//		
//		q.set(18, 18, iono.ionoQ(dt));
//		
//		for (int i = 0; i < this.nsv; i++) {
//			q.set((i+19),(i+19), ure.biasQ());
//		}
//				
//		VectorN r = new VectorN(x.x[0], x.x[1], x.x[2]);
//		VectorN v = new VectorN(x.x[3], x.x[4], x.x[5]);
//		Quaternion w = new Quaternion(x.x[6], x.x[7], x.x[8], x.x[9]);
//		
//		RSW_Frame rsw = new RSW_Frame(r, v);
//		Matrix Cb2i = rsw.ECI2RSW().transpose();
//		Matrix W = w.qMatrix();
//		
//		Matrix g = new Matrix(n);
//		g.setMatrix(3, 13, Cb2i);
//		g.setMatrix(6, 10, W);
////		g.print("gmatrix");
//		
//		Matrix gT = g.transpose();
//		
//		Matrix out = g.times(q.times(gT));				
//		return out;
//	}

	public Matrix Q(double t, double dt, EstSTM stm) {
		int n = this.numberOfStates();
		
		VectorN x = stm.state();
		Matrix phi = stm.phi();
		Matrix phiT = phi.transpose();

		VectorN r = new VectorN(x.x[0], x.x[1], x.x[2]);
		VectorN v = new VectorN(x.x[3], x.x[4], x.x[5]);
		Quaternion w = new Quaternion(x.x[6], x.x[7], x.x[8], x.x[9]);
		
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
		Matrix W = w.qMatrix();


		
		Matrix q = new Matrix(n, n);

    	double accel_sigma = 0.00075;             // in m/s/rt-hr
    	double accel_q = accel_sigma * accel_sigma / 3600.0;    // in (m/s)^2/Hz
		VectorN acc = new VectorN(3);

//		double sp = 1.0E-08;

		acc.set(accel_q);
		VectorN accelx = Cb2i.times(acc);
		q.set(3, 3, accelx.x[0]);
		q.set(4, 4, accelx.x[1]);
		q.set(5, 5, accelx.x[2]);
//		q.set(3, 3, sp);
//		q.set(4, 4, sp);
//		q.set(5, 5, sp);

    	double gyro_sigma = 0.002 * MathUtils.DEG2RAD;             // in rad/rt-hr
    	double gyro_q = gyro_sigma * gyro_sigma / 3600.0;    // in (rad/s)^2/Hz
		VectorN gyro_sig = new VectorN(3);
		gyro_sig.set(gyro_q);
		VectorN gyros = W.times(gyro_sig);
		
		q.set(6, 6, gyros.x[0]);
		q.set(7, 7, gyros.x[1]);
		q.set(8, 8, gyros.x[2]);
		q.set(9, 9, gyros.x[3]);
		
		q.set(10, 10, gyro.Q());
		q.set(11, 11, gyro.Q());		
		q.set(12, 12, gyro.Q());		
		q.set(13, 13, accel.Q());
		q.set(14, 14, accel.Q());		
		q.set(15, 15, accel.Q());

		q.setMatrix(16, 16, rcvr.Q());
		
		q.set(18, 18, iono.Q());
		
		for (int i = 0; i < this.nsv; i++) {
			q.set((i+19),(i+19), ure.Q());
		}
				
		
		Matrix g = new Matrix(n);
		
		g.setMatrix(3, 3, Cb2i);
		
		g.setMatrix(3, 13, Cb2i);
		g.setMatrix(6, 10, W);
//		g.print("gmatrix");
		
		Matrix gT = g.transpose();
		
		Matrix temp = g.times(q.times(gT));	
		
		Matrix out = stm.phi().times(temp.times(phiT));
		out = out.times(dt);			
		return out;
	}


	/**
	 * @see jat.alg.estimators.ProcessModel#propagate(double, double[], double)
	 */
	public double[] propagate(double t0, double[] xin, double tf) {
		
		if ((tf - t0) != dt) {
			System.out.println("propagate step size messed up");
		}

		// fix the attitude quaternion
		EstSTM x = new EstSTM(xin, this.numberOfStates());
		VectorN q = this.getQuat(x);
		q.unitize();
		x.setState(6, q);
		double [] xold = x.longarray();
		double [] xnew = x.longarray();
		
		// propagate
		xnew = rk4.step(t0, xold, this);
		
		// fix the attitude quaternion again
		EstSTM y = new EstSTM(xnew, this.numberOfStates());
		q = this.getQuat(y);
		q.unitize();
		y.setState(6, q);
		double [] out = y.longarray();
		
		return out;
	}
	
	/**
	 * Get the quaternion out of the state vector
	 * @param x EstSTM containing the state
	 * @return VectorN containing the quaternion
	 */
	public VectorN getQuat(EstSTM x){
		VectorN q = x.get(6, 4);
		return q;
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
	 * Close the LinePrinters
	 */
	public void closeLinePrinter(){
		lp1.close();
		lp2.close();
	}
	
	private VectorN wl_0 (VectorN dthetal_1, VectorN dthetal_0){
		VectorN part1 = dthetal_0.times(3.0);
		VectorN w = part1.minus(dthetal_1);
		return w;
	}

	private VectorN wl_1 (VectorN dthetal_1, VectorN dthetal_0){
		VectorN w = dthetal_0.plus(dthetal_1);
		return w;
	}
	
	private VectorN wl_2 (VectorN dthetal_1, VectorN dthetal_0){
		VectorN part1 = dthetal_1.times(3.0);
		VectorN w = part1.minus(dthetal_0);
		return w;
	}
	
	private VectorN dvl_0 (VectorN dvl_1, VectorN dvl_0){
		VectorN part1 = dvl_0.times(3.0);
		VectorN out = part1.minus(dvl_1);
		return out;
	}

	private VectorN dvl_1 (VectorN dvl_1, VectorN dvl_0){
		VectorN out = dvl_0.plus(dvl_1);
		return out;
	}

	private VectorN dvl_2 (VectorN dvl_1, VectorN dvl_0){
		VectorN part1 = dvl_1.times(3.0);
		VectorN out = part1.minus(dvl_0);
		return out;
	}
			

}
