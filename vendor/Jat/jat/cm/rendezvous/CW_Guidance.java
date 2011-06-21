package jat.cm.rendezvous;

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
 * File Created on Aug 26, 2003
 */
 
import jat.cm.ClohessyWiltshire;
import jat.cm.Constants;
import jat.cm.TwoBody;
import jat.matvec.data.*;
import jat.timeRef.*;

/**
 * <P>
 * The CW_Guidance Class provides an rendezvous guidance scheme using the CW and glideslope equations.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
 
public class CW_Guidance {
	
	/** List of targets */
	public CW_TargetList tgtList = new CW_TargetList();
	
	private ClohessyWiltshire cw;
	
	private double tburn1;
	private double tburn2;
	private double tburn3;
	
	/** array containing burn times */
	public double[] burntimes;
	
	/**
	 * Constructor
	 * @param ref TwoBody object containing the orbit
	 * @param tcoast the initial coast time
	 */
	public CW_Guidance(TwoBody ref, double tcoast){
		this.cw = new ClohessyWiltshire(ref);
		this.tburn1 = tcoast;
	}
	
	/**
	 * Compute targets using the glideslope algorithm
	 * @param mc4 Vector containing the mc4 point
	 * @param tof time of flight from Ti to MC4
	 * @param rbar glideslope initiation point
	 * @param tof1 time of flight from MC4 to GS1
	 * @param rtgt target position
	 * @param rhoDot0 initial glideslope speed
	 * @param rhoDotT final safe relative speed
	 * @param tof2 time of flight from GS1 to dock 
	 * @param n number of glideslope burns
	 */
	public void computeTargets(VectorN mc4, double tof, VectorN rbar, double tof1, VectorN rtgt, double rhoDot0, double rhoDotT, double tof2, int n){
		this.tburn2 = this.tburn1 + tof;
		this.tburn3 = this.tburn2 + tof1;
		this.burntimes = new double[n+2];
		burntimes[0] = this.tburn1;
		burntimes[1] = this.tburn2;
		VectorN rmc4 = mc4.copy();
		CW_Target tgt0 = new CW_Target(this.tburn2, rmc4);
		this.tgtList.add(tgt0);
		
		
//		CW_Target tgt0 = new CW_Target(this.tburn2, gs0);
//		this.tgtList.add(tgt0);
		
		// get glideslope initiation point
		VectorN r0 = rbar.copy();
		
		// set up glideslope targetting
		double dt = tof2 / n;
		VectorN rhoZero = r0.minus(rtgt);
		double rho_0 = rhoZero.mag();
		double rho_T = rtgt.mag();
		VectorN u = rhoZero.unitVector();
		double a = (rhoDot0 - rhoDotT)/rho_0;
		
		// get state transition matrix
		Matrix phi = this.cw.phiMatrix(dt);
		Matrix M = phi.getMatrix(0, 2, 0, 2);
		Matrix N = phi.getMatrix(0, 2, 3, 5);
		Matrix S = phi.getMatrix(3, 5, 0, 2);
		Matrix T = phi.getMatrix(3, 5, 3, 5);
		Matrix Ninv = N.inverse();
		
		int i = 2;
		// glideslope targetting		
		for (int m = 0; m < n; m++){
			double tm = m * dt;
			double eatm = Math.exp(a*tm);
			double rho = rho_0 * eatm + (rhoDotT/a)*(eatm - 1.0);
			VectorN temp1 = u.times(rho);
			VectorN rmvec = rtgt.plus(temp1);
			double tburn = this.tburn3 + tm;
			burntimes[i] = tburn;
			i = i + 1;
			CW_Target tgt = new CW_Target(tburn, rmvec);
			this.tgtList.add(tgt);
		}
		double tfinal = this.tburn3 + tof2;
		CW_Target tgt = new CW_Target(tfinal, rtgt);
		this.tgtList.add(tgt);
	}
	
	/**
	 * Compute the delta-v
	 * @param t time
	 * @param r position vector
	 * @param v velocity vector
	 * @param rref reference (for CW frame) position
	 * @param vref reference (for CW frame) velocity
	 * @param index burn index
	 */
	public VectorN computeBurn(double t, VectorN r, VectorN v, VectorN rref, VectorN vref, int index){
		CW_Target tgt = this.tgtList.get(index);
		RSW_Frame rsw = new RSW_Frame(rref, vref);
		VectorN dr_eci = r.minus(rref);
		VectorN dv_eci = v.minus(vref);
		VectorN drv = rsw.transform(dr_eci, dv_eci);
		VectorN dr0 = drv.get(0, 3);
		VectorN dv0 = drv.get(3, 3);
		double tof = tgt.t - t;
		VectorN drtgt = tgt.rtgt;
		VectorN out = cw.intercept(tof, dr0, dv0, drtgt);
		return out;
	}
	
	public static void main(String[] args){
		TwoBody sat =
			new TwoBody(Constants.GM_Earth, 6765500.0, 0.0, 51.8, 0.0, 0.0, 0.0);
			CW_Guidance guid = new CW_Guidance(sat, 1000.0);
			
			VectorN rtgt = new VectorN(3);
			VectorN gs0 = new VectorN(0.0, -183.0, 0.0);
//			VectorN gs0 = new VectorN(-183.0, 0.0, 0.0);
			VectorN mc4 = new VectorN(-549.0, -274.0, 0.0);
			double tof = 4620.0;
			double tof1 = 780.0;
			double tof2 = 3000.0;
			double rhoDot0 = -0.2;
			double rhoDotT = 0.0;
			guid.computeTargets(mc4, tof, gs0, tof1, rtgt, rhoDot0, rhoDotT, tof2, 4);
			guid.tgtList.sendToFile("C:\\Jat\\jat\\output\\vbar_targets.txt");		
	}
	
	

}
