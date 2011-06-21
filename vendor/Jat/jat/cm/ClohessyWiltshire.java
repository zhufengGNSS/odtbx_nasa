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

package jat.cm;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.plot.*;

/**
 * <P>
 * The ClohessyWiltshire Class provides the ability to propagate orbits in the CW system.
 * It also provides a method to compute two impulse rendezvous.
 * This model is only valid for two-body, near circular orbits.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class ClohessyWiltshire implements Derivatives, Printable {

	public TwoBody target;
	private double n;
	public VectorN deltav1;
	public VectorN deltav2;
	public SinglePlot traj_plot = new SinglePlot();
	public VectorList dvlist = new VectorList();

	/**
	 * Constructor for ClohessyWiltshire.
	 */
	public ClohessyWiltshire(TwoBody tgt) {
		this.target = tgt.copy();
		this.n = this.target.meanMotion();
	}

	public double[] derivs(double t, double[] y) {
		double out[] = new double[6];

		double nsq = this.n * this.n;

		out[0] = y[3];
		out[1] = y[4];
		out[2] = y[5];
		out[3] = 3.0 * nsq * y[0] + 2.0 * this.n * y[4];
		out[4] = -2.0 * this.n * y[3];
		out[5] = -1.0 * nsq * y[2];
		//        VectorN x = new VectorN(out);
		//        x.print("x");

		return out;
	}

	/** Derivatives interface for OD Toolbox
	 * 
	 * @param t - time
	 * @param y - input state in CW frame [xyz position & velocity]
	 * @param n - mean orbit rate (radians)
	 * @return time derivatives, units determined by the inputs t and y
	 */
	public static double[] derivs(double t, double[] y, double n) {
		double out[] = new double[6];

		double nsq = n * n;

		out[0] = y[3];
		out[1] = y[4];
		out[2] = y[5];
		out[3] = 3.0 * nsq * y[0] + 2.0 * n * y[4];
		out[4] = -2.0 * n * y[3];
		out[5] = -1.0 * nsq * y[2];

		return out;
	}

	/** Implements the Printable interface to get the data out of the propagator and pass it to the plot.
	 *  This method is executed by the propagator at each integration step.
	 * @param t Time.
	 * @param y Data array.
	 */
	public void print(double t, double[] y) {

		// handle the first variable for plotting - this is a little mystery but it works
		boolean first = true;
		if (t == 0.0)
			first = false;

		// add data point to the plot
		this.traj_plot.plot.addPoint(0, y[1], y[0], first);

		// also print to the screen for warm fuzzy feeling
		System.out.println(t + " " + y[0] + " " + y[1]);
	}

	/**
	 * Compute the two impulse rendezvous
	 * @param tof time of flight
	 * @param tf final time
	 * @param dr0 VectorN containing initial position
	 * @param dv0 VectorN containing initial velocity
	 * @return total delta-v
	 */
	public double twoImpulseRendezvous(
		double tof,
		double tf,
		VectorN dr0,
		VectorN dv0) {
		double nt = this.n * tf;
		double cosnt = Math.cos(nt);
		double sinnt = Math.sin(nt);
		double oon = 1.0 / this.n;
		double omcosnt = 1.0 - cosnt;

		Matrix M = new Matrix(3, 3);
		Matrix N = new Matrix(3, 3);
		Matrix S = new Matrix(3, 3);
		Matrix T = new Matrix(3, 3);

		M.set(0, 0, 4.0 - 3.0 * cosnt);
		M.set(1, 0, 6.0 * (sinnt - nt));
		M.set(1, 1, 1.0);
		M.set(2, 2, cosnt);

		N.set(0, 0, oon * sinnt);
		N.set(0, 1, 2.0 * oon * omcosnt);
		N.set(1, 0, -2.0 * oon * omcosnt);
		N.set(1, 1, oon * (4.0 * sinnt - 3.0 * nt));
		N.set(2, 2, oon * sinnt);

		S.set(0, 0, 3.0 * this.n * sinnt);
		S.set(1, 0, -6.0 * this.n * omcosnt);
		S.set(2, 2, -1.0 * this.n * sinnt);

		T.set(0, 0, cosnt);
		T.set(0, 1, 2.0 * sinnt);
		T.set(1, 0, -2.0 * sinnt);
		T.set(1, 1, 4.0 * cosnt - 3.0);
		T.set(2, 2, cosnt);

		Matrix Ninv = N.inverse();
		VectorN dv0plus = Ninv.times(M.times(dr0));
		dv0plus = dv0plus.times(-1.0);
		this.deltav1 = dv0plus.minus(dv0);

		VectorN temp1 = S.times(dr0);
		VectorN temp2 = T.times(dv0plus);
		VectorN dvfminus = temp1.plus(temp2);
		this.deltav2 = dvfminus.times(-1.0);

		double dvtotal = deltav1.mag() + deltav2.mag();
		return dvtotal;
	}
	
	/**
	 * Compute the one impulse intercept
	 * @param tof time of flight
	 * @param dr0 VectorN containing initial position
	 * @param dv0 VectorN containing initial velocity
	 * @param drtgt VectorN containing final target position
	 * @return intercept delta-v vector
	 */
	public VectorN intercept(double tof, VectorN dr0, VectorN dv0, VectorN drtgt) {
		Matrix phi = this.phiMatrix(tof);
		Matrix M = phi.getMatrix(0, 2, 0, 2);
		Matrix N = phi.getMatrix(0, 2, 3, 5);
		Matrix S = phi.getMatrix(3, 5, 0, 2);
		Matrix T = phi.getMatrix(3, 5, 3, 5);
		Matrix Ninv = N.inverse();
		VectorN temp1 = M.times(dr0);
		VectorN temp2 = drtgt.minus(temp1);
		VectorN dv0plus = Ninv.times(temp2);
		VectorN deltav = dv0plus.minus(dv0);
		
		deltav.print("intercept deltav");

		return deltav;
	}
	
	/**
	 * Return the state transition matrix
	 * @param tof time of flight
	 * @return the state transition matrix
	 */	
	public Matrix phiMatrix (double tof){
		double nt = this.n * tof;
		double cosnt = Math.cos(nt);
		double sinnt = Math.sin(nt);
		double oon = 1.0 / this.n;
		double omcosnt = 1.0 - cosnt;

		Matrix M = new Matrix(3, 3);
		Matrix N = new Matrix(3, 3);
		Matrix S = new Matrix(3, 3);
		Matrix T = new Matrix(3, 3);

		M.set(0, 0, 4.0 - 3.0 * cosnt);
		M.set(1, 0, 6.0 * (sinnt - nt));
		M.set(1, 1, 1.0);
		M.set(2, 2, cosnt);

		N.set(0, 0, oon * sinnt);
		N.set(0, 1, 2.0 * oon * omcosnt);
		N.set(1, 0, -2.0 * oon * omcosnt);
		N.set(1, 1, oon * (4.0 * sinnt - 3.0 * nt));
		N.set(2, 2, oon * sinnt);

		S.set(0, 0, 3.0 * this.n * sinnt);
		S.set(1, 0, -6.0 * this.n * omcosnt);
		S.set(2, 2, -1.0 * this.n * sinnt);

		T.set(0, 0, cosnt);
		T.set(0, 1, 2.0 * sinnt);
		T.set(1, 0, -2.0 * sinnt);
		T.set(1, 1, 4.0 * cosnt - 3.0);
		T.set(2, 2, cosnt);
		
		Matrix out = new Matrix(6,6);
		out.setMatrix(0, 0, M);
		out.setMatrix(0, 3, N);
		out.setMatrix(3, 0, S);
		out.setMatrix(3, 3, T);
		return out;
		
	}
	
	/**
	 * Return the state transition matrix for OD Toolbox.  The time units
	 * are specified by tof and n.
	 * @param tof time of flight (time units)
	 * @param n  mean orbit rate (radians/time unit)
	 * @return the state transition matrix
	 */	
	public static double[][] phiMatrix (double tof, double n){
		double nt = n * tof;
		double cosnt = Math.cos(nt);
		double sinnt = Math.sin(nt);
		double oon = 1.0 / n;
		double omcosnt = 1.0 - cosnt;

		Matrix M = new Matrix(3, 3);
		Matrix N = new Matrix(3, 3);
		Matrix S = new Matrix(3, 3);
		Matrix T = new Matrix(3, 3);

		M.set(0, 0, 4.0 - 3.0 * cosnt);
		M.set(1, 0, 6.0 * (sinnt - nt));
		M.set(1, 1, 1.0);
		M.set(2, 2, cosnt);

		N.set(0, 0, oon * sinnt);
		N.set(0, 1, 2.0 * oon * omcosnt);
		N.set(1, 0, -2.0 * oon * omcosnt);
		N.set(1, 1, oon * (4.0 * sinnt - 3.0 * nt));
		N.set(2, 2, oon * sinnt);

		S.set(0, 0, 3.0 * n * sinnt);
		S.set(1, 0, -6.0 * n * omcosnt);
		S.set(2, 2, -1.0 * n * sinnt);

		T.set(0, 0, cosnt);
		T.set(0, 1, 2.0 * sinnt);
		T.set(1, 0, -2.0 * sinnt);
		T.set(1, 1, 4.0 * cosnt - 3.0);
		T.set(2, 2, cosnt);
		
		double[][] out = new double[6][6];
		
		out[0][0] = M.get(0,0);
		out[0][1] = M.get(0,1);
		out[0][2] = M.get(0,2);
		out[0][3] = N.get(0,0);
		out[0][4] = N.get(0,1);
		out[0][5] = N.get(0,2);
		
		out[1][0] = M.get(1,0);
		out[1][1] = M.get(1,1);
		out[1][2] = M.get(1,2);
		out[1][3] = N.get(1,0);
		out[1][4] = N.get(1,1);
		out[1][5] = N.get(1,2);
		
		out[2][0] = M.get(2,0);
		out[2][1] = M.get(2,1);
		out[2][2] = M.get(2,2);
		out[2][3] = N.get(2,0);
		out[2][4] = N.get(2,1);
		out[2][5] = N.get(2,2);

		out[3][0] = S.get(0,0);
		out[3][1] = S.get(0,1);
		out[3][2] = S.get(0,2);
		out[3][3] = T.get(0,0);
		out[3][4] = T.get(0,1);
		out[3][5] = T.get(0,2);
		
		out[4][0] = S.get(1,0);
		out[4][1] = S.get(1,1);
		out[4][2] = S.get(1,2);
		out[4][3] = T.get(1,0);
		out[4][4] = T.get(1,1);
		out[4][5] = T.get(1,2);
		
		out[5][0] = S.get(2,0);
		out[5][1] = S.get(2,1);
		out[5][2] = S.get(2,2);
		out[5][3] = T.get(2,0);
		out[5][4] = T.get(2,1);
		out[5][5] = T.get(2,2);

		return out;
		
	}
	
	/**
	 * Compute the glideslope
	 * Reference: Hablani, et al
	 * @param r0 initial postion vector
	 * @param v0 initial velocity vector
	 * @param rtgt final target position vector
	 * @param rhoDot0 initial relative velocity
	 * @param rhoDotT safe final relative velocity
	 * @param tof time of flight
	 * @param n number of burns
	 */
	public void glideslope(VectorN r0, VectorN v0, VectorN rtgt, double rhoDot0, double rhoDotT, double tof, int n){
		double dt = tof / n;
		VectorN rhoZero = r0.minus(rtgt);
		double rho_0 = rhoZero.mag();
		double rho_T = rtgt.mag();
		VectorN u = rhoZero.unitVector();
		double a = (rhoDot0 - rhoDotT)/rho_0;
		
		// get state transition matrix
		Matrix phi = this.phiMatrix(dt);
		Matrix M = phi.getMatrix(0, 2, 0, 2);
		Matrix N = phi.getMatrix(0, 2, 3, 5);
		Matrix S = phi.getMatrix(3, 5, 0, 2);
		Matrix T = phi.getMatrix(3, 5, 3, 5);
		
		Matrix Ninv = N.inverse();
		
		VectorList rm = new VectorList();
		
		for (int m = 0; m < n; m++){
			double tm = m * dt;
			double eatm = Math.exp(a*tm);
			double rho = rho_0 * eatm + (rhoDotT/a)*(eatm - 1.0);
			VectorN temp1 = u.times(rho);
			VectorN rmvec = rtgt.plus(temp1);
			rm.add(rmvec);
		}
		rm.add(rtgt);
		
		VectorList vplus = new VectorList();
		VectorList vminus = new VectorList();
		vminus.add(v0);		
		for (int m = 0; m < n; m++){
			VectorN r_m = rm.get(m);
			VectorN temp2 = M.times(r_m);
			VectorN rmp1 = rm.get(m+1);
			VectorN temp3 = rmp1.minus(temp2);
			VectorN vmplus = Ninv.times(temp3);
			vplus.add(vmplus);
			
			VectorN temp4 = S.times(r_m);
			VectorN temp5 = T.times(vmplus);
			VectorN vmp1 = temp4.plus(temp5);
			vminus.add(vmp1);
		}
		
		for (int m = 0; m < n; m++) {
			VectorN vp = vplus.get(m);
			VectorN vm = vminus.get(m);
			VectorN dv = vp.minus(vm);
			dv.print("delta-v");
			this.dvlist.add(dv);
		}
		
	}
	

	public static void main(String[] args) {
		
		// set up the target orbit
		double a = 6765500.0;
		double mu = Constants.GM_Earth;
		TwoBody tgt = new TwoBody( mu, a, 0.0, 0.0, 0.0, 0.0, 0.0);

		// create the CW model		
		ClohessyWiltshire cw = new ClohessyWiltshire(tgt);

		// set up the trajectory plot
//		cw.traj_plot.setTitle("Two Body Trajectory");
//		cw.traj_plot.plot.setXLabel("y (km)");
//		cw.traj_plot.plot.setYLabel("x (km)");

		// create an RungeKutta8 integrator with step-size of 10 sec
		RungeKutta8 rk8 = new RungeKutta8(1.0);
		
		// set up a LinePrinter for printing results to a file
//		LinePrinter lp = new LinePrinter("c:\\temp\\", "sts-iss.dat");


		// set up the initial conditions
		VectorN dr0 = new VectorN(-16.628475931473076, -14999.987710825975, 0.0);
		VectorN dv0 = new VectorN(3);
		
//		VectorN drtgt = new VectorN(-183.0, 0.0, 0.0);
		VectorN drtgt = new VectorN(0.0, -183.0, 0.0);
		
		// set up the times
		double tof = 5000.0;
		double tcoast = 0.0;
		double t0 = 0.0;
		double tf = tcoast + tof;


		// initialize the variables
		double[] x0 = new double[6];
		x0[0] = dr0.x[0];
		x0[1] = dr0.x[1];
		x0[2] = dr0.x[2];
		x0[3] = dv0.x[0];
		x0[4] = dv0.x[1];
		x0[5] = dv0.x[2];

		// integrate the equations to tcoast
		double[] rv = rk8.integrate(t0, x0, tcoast, cw, cw, true);
				
		VectorN drminus = new VectorN(rv[0], rv[1], rv[2]);
		VectorN dvminus = new VectorN(rv[3], rv[4], rv[5]);
		
		// compute the delta-v's		
//		double dvtot = cw.twoImpulseRendezvous(tof, tof, drminus, dvminus);
//		System.out.println("dvtot = " + dvtot);
//		cw.deltav1.print("delta-v1");
//		cw.deltav2.print("delta-v2");
//		
//		// add in the effects of the first deltav
//		VectorN dv0plus = dvminus.plus(cw.deltav1);
//		x0[0] = drminus.x[0];
//		x0[1] = drminus.x[1];
//		x0[2] = drminus.x[2];
//		x0[3] = dv0plus.x[0];
//		x0[4] = dv0plus.x[1];
//		x0[5] = dv0plus.x[2];

		VectorN dv = cw.intercept(tof, dr0, dv0, drtgt);
		VectorN dv0plus = dvminus.plus(dv);
		x0[0] = drminus.x[0];
		x0[1] = drminus.x[1];
		x0[2] = drminus.x[2];
		x0[3] = dv0plus.x[0];
		x0[4] = dv0plus.x[1];
		x0[5] = dv0plus.x[2];
		

		// integrate the equations from tcoast to tf       
		double [] xnew = rk8.integrate(tcoast, x0, tf, cw, cw, true);
		
		VectorN r0 = new VectorN(xnew[0], xnew[1], xnew[2]);
		VectorN v0 = new VectorN(xnew[3], xnew[4], xnew[5]);
		VectorN rtgt = new VectorN(0.0, 0.0, 0.0);
		double rhoDot0 = -0.2;
		double rhoDotT = 0.0;
		double tglide = 3000.0;
		
		cw.glideslope(r0, v0, rtgt, rhoDot0, rhoDotT, tglide, 4);
		
		int index = 0; 
		t0 = tf;
		tf = tf + 750.0;
		while (cw.dvlist.hasNext(index)){
			VectorN deltav = cw.dvlist.get(index);
			x0[0] = xnew[0];
			x0[1] = xnew[1];
			x0[2] = xnew[2];
			x0[3] = xnew[3] + deltav.x[0];
			x0[4] = xnew[4] + deltav.x[1];
			x0[5] = xnew[5] + deltav.x[2];
			
			xnew = rk8.integrate(t0, x0, tf, cw, cw, true);
			t0 = t0 + 750.0;
			tf = tf + 750.0;
			index = index + 1;
		}

		// make the plot visible
		cw.traj_plot.setVisible(true);

		// close the LinePrinter
//		lp.close();

	}

}
