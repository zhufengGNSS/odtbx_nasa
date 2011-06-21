package jat.alg.estimators;

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
 * File Created on May 7, 2003
 */
import jat.matvec.data.*;
//import jat.audio.*;

/**
 * The ExtendedKalmanFilter Class processes measurements using an EKF algorith,
 * given the measurements, measurement model and dynamics (or process) model.
* Assumes a scalar measurement update.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class ExtendedKalmanFilter {

	/** Measurement Model */
	public MeasurementModel measModel;

	/** Measurement Data */
	public MeasurementData meas;

	/** Dynamics or Process Model */
	public ProcessModel process;

	/** Number of states or unknowns */
	public int n;

	/** Nominal time step */
	private double dtNominal = 1.0;

	/**
	 * Constructor.
	 * @param mm MeasurementModel
	 * @param pm ProcessModel
	 * @param rp LinePrinter
	 */
	public ExtendedKalmanFilter(
		MeasurementModel mm,
		MeasurementData md,
		ProcessModel pm) {
		this.measModel = mm;
		this.meas = md;
		this.process = pm;
		this.n = pm.numberOfStates();

	}

	private Matrix updateCov(VectorN k, VectorN h, Matrix p) {
		Matrix eye = new Matrix(this.n);
		Matrix kh = k.outerProduct(h);
		Matrix i_kh = eye.minus(kh);
		Matrix out = i_kh.times(p);
		return out;
	}

	private Matrix updateCovJoseph(VectorN k, VectorN h, Matrix p) {
		Matrix eye = new Matrix(this.n);
		Matrix kh = k.outerProduct(h);
		Matrix i_kh = eye.minus(kh);
		Matrix i_khT = i_kh.transpose();
		double r = measModel.R();
		Matrix kkT = k.outerProduct(k);
		Matrix krkT = kkT.times(r);

		Matrix part1 = i_kh.times(p);
		part1 = part1.times(i_khT);
		Matrix out = part1.plus(krkT);
		return out;
	}

	private VectorN updateState(VectorN k, double z) {
		VectorN xhat = k.times(z);
		return xhat;
	}

	private Matrix propCov(Matrix p, Matrix phi, Matrix q) {
		Matrix phitrans = phi.transpose();
		Matrix temp1 = p.times(phitrans);
		Matrix temp2 = phi.times(temp1);
		Matrix out = temp2.plus(q);
		return out;
	}

	private VectorN kalmanGain(Matrix p, VectorN h, double r) {
		VectorN ph = p.times(h);
		VectorN hp = h.times(p);
		double hph = hp.dotProduct(h);
		double hphr_inv = 1.0 / (hph + r);
		VectorN out = ph.times(hphr_inv);
		
		return out;
	}

	/** Process the measurements
	 */
	public void process() {

		System.out.println("EKF Processing ...");

		// initialize
		double tprev = 0.0;
		Matrix pold = process.P0();
		Matrix pnew = pold.copy();
		EstSTM xref = new EstSTM(process.xref0());
		double[] xprev = xref.longarray();
		VectorN k = new VectorN(process.numberOfStates());

		//		xref.print("ic");

		int i = 0;
		//		for (int kk = 0; kk < 12; kk++){
//		while (meas.hasNext(i) && (tprev < 80.0)) {
		while (meas.hasNext(i)) {

			//            System.out.println("Processing Measurement "+(i+1));

			// read in next measurement
			double t = meas.time(i);
			double z = meas.z(i);

			// move to new time
			double dt = t - tprev;

			// detect backwards time jump
			if (dt < 0.0) {
				System.out.println("backwards time jump");
				System.exit(1);
			}

			// propagate state and covariance to new time
			if (dt > 0.0) {

				while (tprev < t) {
					
					// print output to line printer
					process.print(tprev, xref.state(), pold);

					double tnext = tprev + this.dtNominal;

					if (tnext > t) {
						tnext = t;
						System.out.println(
							"measurement gap not an integer number of time steps");
					}

					double[] xnew = process.propagate(tprev, xprev, tnext);
					xref = new EstSTM(xnew, this.n);
					Matrix phi = xref.phi();
					//				phi.print("phi");
					Matrix q = process.Q(tnext, this.dtNominal, xref);
					//				q.print("q matrix");
					//				pold.print("p before prop");
					pnew = this.propCov(pold, phi, q);
					//				pnew.print("p after prop");
					tprev = tnext;
					xref.resetPhi();
					xprev = xref.longarray();
					pold = pnew.copy();
				}
			} else {
				pnew = pold.copy(); // dt = 0, no change to covariance
			}

			//			xref.state.print("state before processing");
			//			pnew.print("covariance before update");

			// measurement update
			double zpred = measModel.zPred(i, t, xref.state());
			double y = z - zpred;
			double r = measModel.R();
			
			// resid check at 9 sigma level
//			if (Math.abs(y) < (9.0*Math.sqrt(r))) {
			
				//			System.out.println("predicted range: "+zpred+" residual = "+y);
	
				VectorN h = measModel.H(xref.state());
				//			h.print("h");
	
				// compute the Kalman gain
				k = this.kalmanGain(pnew, h, r);
				//			k.print("kalman gains");
	
				// compute new best estimate
				VectorN xhat = k.times(y);
				//			xhat.print("corrections");
				
				// update state and covariance
				xref.update(xhat); // extended Kalman filter
				//            xref.print("updated state");
				//			pnew.print("p before update");
				//System.out.println(xref.state().toString());
				pold = this.updateCov(k, h, pnew);
				//			pold.print("p after update");
				//			VectorN sigmas = pold.diagonal();
				//			sigmas = sigmas.ebeSqrt();
	
				// check the update
				double zafter = measModel.zPred(i, t, xref.state());
				double yafter = z - zafter;
				process.printResiduals(t, y, yafter);

	//			VectorN printvec = k.append(y);
	//			printvec = printvec.append(yafter);
	//			lp2.print(t, printvec.x);
	//			if (t >= 1798.0){
	//				System.out.println(t+"\t"+y+"\t"+yafter);
	//			}
	
				// print output to line printer
				//			VectorN printvector = new VectorN(xref.state(), sigmas);
				//			lp.print(t, printvector.x);
				//			System.out.println(t+"\t"+h.toString());

//			}
			// get ready for next data point
			
			i = i + 1;
			xref.resetPhi(); // re-linearize
			tprev = t;
			xprev = xref.longarray();
		}
		
		process.closeLinePrinter();

		System.out.println("EKF Processing Completed");
//		SoundPlayer.play("C:\\Jat\\jat\\jat\\audio\\sounds\\itisdone.wav");

	}

}
