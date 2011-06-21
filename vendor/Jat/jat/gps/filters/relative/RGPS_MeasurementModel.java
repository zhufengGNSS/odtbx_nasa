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
//import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.gps.*;
//import jat.gps_ins.*;
//import jat.gps_ins.*;
//import jat.cm.*;

/**
* The RGPS_MeasurementModel.java Class provides ...
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RGPS_MeasurementModel implements MeasurementModel {

	private VectorN htilde;

	private RGPS_MeasurementList gpsData;
	
	private GPS_Constellation constell;
			
	private static final double rangeSigma = 4.0;
	
	private static final double cpSigma = 0.02;
			
	private int index = 0;
	
	private int numberOfStates;
				
	private IonoModel iono = new IonoModel();
	
	private int nsv;
	
	private int type = 0;
	
	private double rangeA = 0.0;
	
	private double rangeB = 0.0;
	
	private static final int stateIndexA = 0;
	private static final int clockIndexA = 6;
	private int stateIndexB;
	private int clockIndexB;
	private int iaIndex;
	
	private boolean[] firstTime;

	
	public RGPS_MeasurementModel(RGPS_MeasurementList gps, GPS_Constellation con){
		this.gpsData = gps;
		this.constell = con;
		this.nsv = con.size();
		this.numberOfStates = 19 + 2 * this.nsv;
		this.stateIndexB = this.nsv + 10;
		this.clockIndexB = this.stateIndexB + 6;
		this.iaIndex = this.nsv + 19;
		this.firstTime = new boolean[this.nsv];
		for (int i = 0; i < this.nsv; i++) {
			this.firstTime[i] = true;
		}
		
	}

	/**
	 * @see jat.alg.estimators.MeasurementModel#H(VectorN)
	 */
	public VectorN H(VectorN xref) {
		return this.htilde;
	}

	/**
	 * @see jat.alg.estimators.MeasurementModel#R()
	 */
	public double R() {
		double r = rangeSigma*rangeSigma;
		if (type == 2) r = 2.0*cpSigma*cpSigma;
		return r;
	}

	/**
	 * @see jat.alg.estimators.MeasurementModel#zPred(int, double, VectorN)
	 */
	public double zPred(int index, double t, VectorN xref) {
		
		double out = 0.0;
		
		// get the SVID and MJD time of measurement
		RGPS_Measurement meas = gpsData.get(index);
		int prn = meas.svid;
		int svindex = this.constell.getIndex(prn);
		this.type = meas.type();
		double t_mjd = meas.t_mjd();
		double range = meas.range();
		
		//time check
		double tgps = meas.t();
		if (tgps != t) {
			System.out.println("zPred: time sync error, meas time = "+tgps+" tsim = "+t);
		}
		
		if (type == 0) {
			out = rangePred(t_mjd, xref, prn, type);
			this.rangeA = out;
		}
		
		if (type == 1) {
			out = rangePred(t_mjd, xref, prn, type);
			this.rangeB = out;
		}
		
		if (type == 2) {
			out = cpPred(t_mjd, xref, prn, range);
		}
				
		return out;
	}
	
	private double rangePred(double t_mjd, VectorN xref, int prn, int type ){
		
		// set indices
		int svindex = this.constell.getIndex(prn);
		int stateIndex = stateIndexA;
		int clockIndex = clockIndexA;
		if (type == 1) {
			stateIndex = this.stateIndexB;
			clockIndex = this.clockIndexB;
		}

		// extract the current spacecraft position vector		
		VectorN r = xref.get(stateIndex, 3);
		VectorN v = xref.get((stateIndex+3), 3);
				
		// get the GPS SV from the constellation
		GPS_SV sv = constell.getPRN(prn);
		
		// compute the time of transmission
		double ts_mjd = GPS_Utils.transmitTime(t_mjd, sv, r);
		
		// compute the GPS SV position vector at transmission time
		VectorN rvGPS = sv.rvECI(ts_mjd);
		VectorN rGPS = new VectorN(rvGPS.x[0], rvGPS.x[1], rvGPS.x[2]);
		VectorN vGPS = new VectorN(rvGPS.x[3], rvGPS.x[4], rvGPS.x[5]);
		
		// compute the LOS vector
		VectorN los = rGPS.minus(r);
		VectorN losu = los.unitVector();
		
		// compute the expected range
		double range = los.mag();
				
		// compute the range rate
//		VectorN vrel = vGPS.minus(v);
		double range_rate = GPS_Utils.rangeRate(los, v, vGPS);
		
		// add clock bias contribution to range
		double bc = xref.x[clockIndex];
		double omrroc = 1.0 - (range_rate/GPS_Utils.c);
		double bcpart = omrroc*bc;
		range = range + bcpart;
		
		// add iono contribution to range
		int iono_index = 9;
		double diono = xref.x[iono_index];
		double iv = (1.0 + diono)*iono.Ivbar;
		double elev = GPS_Utils.elevation(r, rGPS);
		double ionopart = IonoModel.del_iono(iv, elev);
		range = range + ionopart;
		
		// add ure contribution to range
		int ure_index = 10 + svindex;
		double dure = xref.x[ure_index];
		range = range + dure;
		
		// compute the H vector
		this.htilde = new VectorN(this.numberOfStates);
		VectorN drrdr = this.drrdr(los, vGPS, range_rate);
		double term1 = -1.0 * bc / GPS_Utils.c;
		drrdr = drrdr.times(term1);
		VectorN drrdv = this.drrdv(los, vGPS);
		VectorN drhodv = drrdv.times(term1);
		VectorN drhodr = drrdr.minus(losu);
		htilde.set(stateIndex, drhodr);
		htilde.set((stateIndex+3), drhodv);
		double drdbc = 1.0 - range_rate/GPS_Utils.c;
		htilde.set(clockIndex, drdbc); 

		double drdiono = IonoModel.del_iono(iono.Ivbar, elev);
		htilde.set(iono_index, drdiono);
		htilde.set(ure_index, 1.0);
		
		return range;
	}
		
	private double cpPred(double t_mjd, VectorN xref, int prn, double meas){
		
		// set indices
		int svindex = this.constell.getIndex(prn);

		// extract the current spacecraft position vector		
		VectorN r = xref.get(stateIndexA, 3);
		VectorN v = xref.get((stateIndexA+3), 3);
		VectorN rISS = xref.get(this.stateIndexB, 3);
		VectorN vISS = xref.get((this.stateIndexB+3), 3);
				
		// get the GPS SV from the constellation
		GPS_SV sv = constell.getPRN(prn);
		
		// compute the time of transmission
		double ts_mjd1 = GPS_Utils.transmitTime(t_mjd, sv, r);
		double ts_mjd2 = GPS_Utils.transmitTime(t_mjd, sv, rISS);
		
		// compute the GPS SV position vector at transmission time
		VectorN rvGPS1 = sv.rvECI(ts_mjd1);
		VectorN rGPS1 = new VectorN(rvGPS1.x[0], rvGPS1.x[1], rvGPS1.x[2]);
		VectorN vGPS1 = new VectorN(rvGPS1.x[3], rvGPS1.x[4], rvGPS1.x[5]);
		VectorN rvGPS2 = sv.rvECI(ts_mjd2);
		VectorN rGPS2 = new VectorN(rvGPS2.x[0], rvGPS2.x[1], rvGPS2.x[2]);
		VectorN vGPS2 = new VectorN(rvGPS2.x[3], rvGPS2.x[4], rvGPS2.x[5]);
		
		// compute the LOS vector
		VectorN los1 = rGPS1.minus(r);
		VectorN losu1 = los1.unitVector();
		VectorN los2 = rGPS2.minus(rISS);
		VectorN losu2 = los2.unitVector();
		
		// compute the expected range
		double range1 = los1.mag();
		double range2 = los2.mag();
				
		// compute the range rate
//		VectorN vrel1 = vGPS1.minus(v);
		double range_rate1 = GPS_Utils.rangeRate(los1, v, vGPS1);
//		VectorN vrel2 = vGPS2.minus(vISS);
		double range_rate2 = GPS_Utils.rangeRate(los2, vISS, vGPS2);
		
		// add clock bias contribution to range
		double bc1 = xref.x[clockIndexA];
		double omrroc1 = 1.0 - (range_rate1/GPS_Utils.c);
		double bcpart1 = omrroc1*bc1;
		range1 = range1 + bcpart1;
		
		double bc2 = xref.x[clockIndexB];
		double omrroc2 = 1.0 - (range_rate2/GPS_Utils.c);
		double bcpart2 = omrroc2*bc2;
		range2 = range2 + bcpart2; 
		
		double drange = range1 - range2;
		
		// add integer ambiguity to range
		int ia_index = this.iaIndex + svindex;
		double dia = xref.x[ia_index];
		double ia_est = meas - this.rangeA + this.rangeB;
		
		if (firstTime[svindex]){                  // initialize integer amb state and cov
			xref.x[ia_index] = ia_est;
			dia = ia_est;
			firstTime[svindex] = false;
		}
		
//		double ia_err = Math.abs(dia - ia_est);
//		if (ia_err > 10.0) {
//			xref.x[ia_index] = ia_est;
//			dia = ia_est;
//			System.out.println("reinit integers");
//		}
		
		drange = drange + dia;
		
		
		// compute the H vector
		this.htilde = new VectorN(this.numberOfStates);
				
		VectorN drrdr1 = this.drrdr(los1, vGPS1, range_rate1);
		VectorN drrdr2 = this.drrdr(los2, vGPS2, range_rate2);
		
		double term1 = -1.0 * bc1 / GPS_Utils.c;
		double term2 = -1.0 * bc2 / GPS_Utils.c;
		
		drrdr1 = drrdr1.times(term1);
		drrdr2 = drrdr2.times(term2);		

		VectorN drrdv1 = this.drrdv(los1, vGPS1);
		VectorN drrdv2 = this.drrdv(los2, vGPS2);		

		VectorN drhodr1 = drrdr1.minus(losu1);
		VectorN drhodr2 = drrdr2.minus(losu2);
		drhodr2 = drhodr2.times(-1.0);		

		htilde.set(stateIndexA, drhodr1);
		htilde.set(stateIndexB, drhodr2);
		
		VectorN drhodv1 = drrdv1.times(term1);
		VectorN drhodv2 = drrdv2.times(term2);
		drhodv2 = drhodv2.times(-1.0);
		htilde.set((stateIndexA+3), drhodv1);
		htilde.set((this.stateIndexB+3), drhodv2);
		
		double drdbc1 = 1.0 - range_rate1/GPS_Utils.c;
		double drdbc2 = range_rate2/GPS_Utils.c - 1.0;  // subtract it
		
		htilde.set(clockIndexA, drdbc1); 
		htilde.set(this.clockIndexB, drdbc2);
		
		htilde.set(ia_index, 1.0);
		
		
		return drange;
	}


	
	private VectorN drrdr (VectorN los, VectorN vGPS, double rr){
		VectorN vc = vGPS.divide(GPS_Utils.c);
		double range = los.mag();
		double corr = los.dotProduct(vc);
		double factor = 1.0 / (range + corr);
		VectorN losu = los.unitVector();
		VectorN part2 = losu.plus(vc);
		part2 = part2.times(rr);
		VectorN numerator = part2.minus(los);
		VectorN out = numerator.times(factor);
		return out;
	}

	private VectorN drrdv (VectorN los, VectorN vGPS){
		VectorN vc = vGPS.divide(GPS_Utils.c);
		double range = los.mag();
		double corr = los.dotProduct(vc);
		double factor = -1.0 / (range + corr);
		VectorN out = los.times(factor);
		return out;
	}
		
		

}
