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
 
package jat.gps.filters.absolute;
import jat.matvec.data.*;
import jat.alg.estimators.*;
import jat.alg.integrators.*;
//import jat.cm.*;
import jat.gps.*;
import jat.gps.filters.*;
//import jat.gps_ins.*;
import jat.forces.*;
import jat.forces.density.earth.*;
import jat.forces.gravity.earth.J2Gravity;
import jat.timeRef.*;

/**
 * The GPS_EOM provides the equations of motion for the 
 * GPS-only EKF without a thruster model.
 * 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GPS_EOM implements Derivatives {

	private IonoModel iono;
	private ReceiverFilterModel rcvr;
	private double t_mjd0 = 51969.0;

	private double stsMass = 104328.0;
	private double stsArea = 454.4;
	private double stsCd = 2.0;
	private CIRA_ExponentialDrag sts_ced = new CIRA_ExponentialDrag(this.stsCd, this.stsArea, this.stsMass);

	private URE_Model ure;
	private int nsv;
	private int numberOfStates;
	
	/**
	 * Constructor
	 * @param nsat number of GPS SVs
	 * @param nstates number of states
	 * @param io IonoModel
	 * @param rc ReceiverFilterModel
	 * @param ur URE_Model
	 */
	public GPS_EOM(int nsat, int nstates, IonoModel io, ReceiverFilterModel rc, URE_Model ur) {
		this.nsv = nsat;
		this.numberOfStates = nstates;
		this.ure = ur;
		this.iono = io;
		this.rcvr = rc;
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
		
						
				
		// position derivatives
		out.set(0, v);
		
		// velocity derivatives
		J2Gravity j2chaser = new J2Gravity(r);
		VectorN g = j2chaser.local_gravity();	

		double Mjd = this.t_mjd0 + t/86400.0;
        EarthRef ref = new EarthRef(Mjd);

		sts_ced.compute(ref, r, v);
		VectorN sts_drag0 = sts_ced.dragAccel();
		double dragfactor1 = 1.0 + stsdrag;
		VectorN drag = sts_drag0.times(dragfactor1);
		
		VectorN vdot = g.plus(drag);

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
