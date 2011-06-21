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
 * File Created on Sep 25, 2003
 */
 
package jat.gps.filters.relative;
import jat.matvec.data.*;
import jat.alg.estimators.*;
import jat.alg.integrators.*;
import jat.cm.*;
import jat.gps.*;
import jat.gps.filters.*;
//import jat.gps_ins.*;
import jat.forces.density.earth.*;
import jat.forces.gravity.earth.J2Gravity;
import jat.timeRef.*;

/**
 * The RGPS_Thruster_EOM provides the equations of motion for the 
 * RGPS-only EKF including a thruster model.
 * 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class RGPS_Thruster_EOM implements Derivatives {

	private IonoModel iono;
	private ReceiverFilterModel rcvr;
	private double t_mjd0 = 51969.0;
	private double issMass = 128990.0;
	private double issArea = 640.7;
	private double issCd = 2.35;
	private CIRA_ExponentialDrag iss_ced = new CIRA_ExponentialDrag(this.issCd, this.issArea, this.issMass);

	private double stsMass = 104328.0;
	private double stsArea = 454.4;
	private double stsCd = 2.0;
	private CIRA_ExponentialDrag sts_ced = new CIRA_ExponentialDrag(this.stsCd, this.stsArea, this.stsMass);

	private URE_Model ure;
	private int nsv;
	private int issIndex;
	private int numberOfStates;
	
	private VectorN fburn;
		
	/**
	 * Constructor
	 * @param nsat number of GPS SVs
	 * @param nstates number of states
	 * @param io IonoModel
	 * @param rc ReceiverFilterModel
	 * @param ur URE_Model
	 */
	public RGPS_Thruster_EOM(int nsat, int iss, int nstates, IonoModel io, ReceiverFilterModel rc, URE_Model ur) {
		this.nsv = nsat;
		this.issIndex = iss;
		this.numberOfStates = nstates;
		this.ure = ur;
		this.iono = io;
		this.rcvr = rc;
	}

	public void setBurn(FiniteBurn b) {
		double thrust = b.accel;
		VectorN unit = b.unitVector;
		this.fburn = unit.times(thrust);		
	}
	
	/**
	 * @see jat.alg.integrators.Derivatives#derivs(double, double[])
	 */
	public double[] derivs(double t, double[] x) {
		
		int n = this.numberOfStates;
		
		VectorN out = new VectorN(x.length);

		
		// strip out the incoming data
		VectorN r = new VectorN(x[0], x[1], x[2]);
		VectorN v = new VectorN(x[3], x[4], x[5]);
		EstSTM stm = new EstSTM(x, n);
		Matrix phi = stm.phi();
		
		
		// strip off clock states		
		VectorN clock1 = new VectorN(2);
		clock1.set(0, x[6]);
		clock1.set(1, x[7]);
		
		double stsdrag = x[8];

		double del_iono = x[9];  // iono state
		
		// strip off incoming ure states
		VectorN urevec = new VectorN(this.nsv);
		for (int i = 0; i < this.nsv; i++) {
			urevec.set(i, x[i+10]);
		}
		
		// strip off ISS states
		int kk = this.issIndex;
		VectorN rISS = new VectorN(x[kk], x[kk+1], x[kk+2]);
		VectorN vISS = new VectorN(x[kk+3], x[kk+4], x[kk+5]);
		VectorN clock2 = new VectorN(2);
		clock2.set(0, x[kk+6]);
		clock2.set(1, x[kk+7]);
		double ddrag = x[kk+8];
						
				
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

		sts_ced.compute(ref, r, v);
		VectorN sts_drag0 = sts_ced.dragAccel();
		double dragfactor1 = 1.0 + stsdrag;
		VectorN drag = sts_drag0.times(dragfactor1);

        // acceleration due to the burn
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
        
//        Matrix Cb2i = q.quat2DCM();
        VectorN thrust = Cb2i.times(this.fburn);
        
        VectorN sf = drag.plus(thrust);
		
		VectorN vdot = g.plus(sf);

		out.set(3, vdot);
		
		
		// GPS clock model derivatives
		VectorN bcdot = rcvr.biasProcess(clock1);
		out.set(6, bcdot);

		double dragdot = DragProcessModel.dragProcess(stsdrag);
		out.set(8, dragdot);
		
		// iono derivs
		double ionodot = iono.ionoProcess(del_iono);
		out.set(9, ionodot);
		
		//ure derivs
		VectorN uredot = ure.ureProcess(urevec);
		out.set(10, uredot);

		// position derivatives
		out.set(kk, vISS);
		
		// velocity derivatives
//		TwoBody orbit2 = new TwoBody(Constants.GM_Earth, rISS, vISS);
//		g = orbit2.local_grav();		
		J2Gravity j2iss = new J2Gravity(rISS);
		g = j2iss.local_gravity();
		
//		g = jgm3.gravity(rISS, E);

//		double Mjd = this.t_mjd0 + t/86400.0;
//        EarthRef ref = new EarthRef(Mjd);
		iss_ced.compute(ref, rISS, vISS);
		VectorN iss_drag0 = iss_ced.dragAccel();
		double dragfactor2 = 1.0 + ddrag;
		drag = iss_drag0.times(dragfactor2);
		
		vdot = g.plus(drag);	
		out.set((kk+3), vdot);
		
		// GPS clock model derivatives
		VectorN bcdot2 = rcvr.biasProcess(clock2);
		out.set((kk+6), bcdot2);
		
		dragdot = DragProcessModel.dragProcess(ddrag);
		out.set((kk+8), dragdot);
		
		// integer ambiguity derivs = 0		

		// A matrix
		Matrix A = new Matrix(n, n);
		
		// position rows
		Matrix eye = new Matrix(3);
		A.setMatrix(0, 3, eye);
		
		// velocity rows
		Matrix G = j2chaser.gravityGradient();
		Matrix D = sts_ced.partialR().times(dragfactor1);
		Matrix GD = G.plus(D);
		A.setMatrix(3, 0, GD);
		
		D = sts_ced.partialV().times(dragfactor1);
		A.setMatrix(3, 3, D);
		
		// partials of drag accel wrt drag state
		A.set(3, 8, sts_drag0.x[0]);
		A.set(4, 8, sts_drag0.x[1]);
		A.set(5, 8, sts_drag0.x[2]);
		
        
        //clock drift row
        A.set(6, 7, 1.0);
        
        // drag rows
        double tau_drag = -1.0/DragProcessModel.correlationTime;
        A.set(8, 8, tau_drag);
        
        
        // iono row
        double tau_iono = -1.0/iono.correlationTime;
        A.set(9, 9, tau_iono);
        
        // ure part
        Matrix bigeye = new Matrix(this.nsv);
        Matrix tau_ure = eye.times(-1.0/URE_Model.correlationTime);
        A.setMatrix(10, 10, tau_ure);
        
		// position rows
		A.setMatrix(kk, (kk+3), eye);
		
		// velocity rows
		G = j2iss.gravityGradient();
		D = iss_ced.partialR().times(dragfactor2);
		GD = G.plus(D);
		A.setMatrix((kk+3), kk, GD);
		
		D = iss_ced.partialV().times(dragfactor2);
		A.setMatrix((kk+3),(kk+3), D);
		
		// partials of drag accel wrt drag state
		A.set((kk+3),(kk+8), iss_drag0.x[0]);
		A.set((kk+4),(kk+8), iss_drag0.x[1]);
		A.set((kk+5),(kk+8), iss_drag0.x[2]);
        
        //clock drift row
        A.set((kk+6),(kk+7), 1.0);
        
        // drag rows
        A.set((kk+8), (kk+8), tau_drag);
        		
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
	

}
