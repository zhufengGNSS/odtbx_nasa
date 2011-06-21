/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2009 United States Government as represented by the
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
 */
package jat.cm;

import jat.alg.integrators.Printable;
import jat.alg.integrators.SecondDerivatives;
import jat.matvec.data.*;
import jat.alg.integrators.*;
import jat.eph.DE405;
import jat.eph.DE405_Body;
import jat.forces.gravity.*;
import jat.constants.*;
import jat.spacetime.*;


public class Encke implements SecondDerivatives {


	private VectorN r_osc = null;
	private VectorN v_osc = null;
	private VectorN rv = null;
	private double _mu = 0.0;
	private TwoBody tb = null;
	private double step_size = 60.0;
	private Nystrom4 nystrom = new Nystrom4(step_size);
	private DE405 de405 = new DE405();
	private ThirdBodyGravity tbg = new ThirdBodyGravity(de405, DE405_Body.EARTH, DE405_Body.MOON, IERS_1996.GM_Moon);
	private Time time0 = null;
	public int rectify_counter = 0;
	private boolean rectify = false;
	private boolean once = true;
	public double devmag = 0.0;
	
	public Encke(double mu, Time t0, VectorN r, VectorN v){
		r_osc = r.copy();
		v_osc = v.copy();
		rv = r_osc.append(v_osc);
		_mu = mu;
		tb = new TwoBody(mu, r_osc, v_osc);
		time0 = t0;
		
	}
	public double[] derivs(double t, double[] x, double[] xdot) {
		// get the deviations
		VectorN del = new VectorN(x);
		VectorN nu = new VectorN(xdot);
		// compute current orbit
		VectorN r = r_osc.plus(del);
		VectorN v = v_osc.plus(nu);
				
		// form first term
		double rho2 = r_osc.mag() * r_osc.mag();
		double rho3 = rho2 * r_osc.mag();		
		double muorho3 = _mu/rho3;
		
		double r2 = r.mag() * r.mag();
		double r3 = r2 * r.mag();
		double omr3 = 1.0 - (rho3/r3);
		
		VectorN term1 = r.times(omr3);
		VectorN term2 = (term1.minus(del)).times(muorho3);
		

		// compute the perturbations
		VectorN accel_d = perts(t, r.getArray(), v.getArray());
		double pertmag = accel_d.mag();
		// make the output
//		VectorN out = sum.plus(accel_d);
		VectorN out = accel_d.plus(term2);
		// do rectification check
		if (term2.mag() > 0.01*pertmag) rectify = true;
		
		return out.getArray();
	}	
	
//	public double[] derivs(double t, double[] x, double[] xdot) {
//		// get the deviations
//		VectorN del = new VectorN(x);
//		VectorN nu = new VectorN(xdot);
//		// compute current orbit
//		VectorN r = r_osc.plus(del);
//		VectorN v = v_osc.plus(nu);
//		// form first term
//		double rho2 = r_osc.mag() * r_osc.mag();
//		double rho3 = rho2 * r_osc.mag();		
//		double muor3 = -1.0*_mu/rho3;
//		// form first term
//		VectorN term1 = del.times(muor3);
//		// form second term
//		double num = del.dotProduct(del.minus(r.times(2.0)));
//		double r2 = Math.pow(r.mag(),2.0);
//		double q = num/r2;
//		double f = q*(3.0 + 3.0*q + q*q)/(1.0 + Math.pow(1.0+q, 1.5));
//		double fmuor3 = f*muor3;		
//		VectorN term2 = r.times(fmuor3);
//		// compute the perturbations
//		VectorN accel_d = perts(t, r.getArray(), v.getArray());
//		double pertmag = accel_d.mag();
//		// make the output
//		VectorN sum = term1.plus(term2);
////		VectorN out = sum.plus(accel_d);
//		VectorN out = accel_d.copy();
//		// do rectification check
//		if (sum.mag() > pertmag) rectify = true;
//		
//		return out.getArray();
//	}
	
	public VectorN perts (double t, double[] r, double[] v){
		Time time = time0.plus(t);
		VectorN rin = new VectorN(r);
		VectorN out = tbg.acceleration(time, rin);
//		VectorN out = new VectorN(3);
		return out;
	}
	
	private void updateOsculatingOrbit(double dt){
		tb.step(dt);

		// update osculating position
		r_osc = tb.getR();
				
		// update osculating velocity
		v_osc = tb.getV();

	}
	
	private void rectify(){
		r_osc = rv.get(0, 3);
		v_osc = rv.get(3, 3);
//		r0 = r_osc.copy();
//		v0 = v_osc.copy();

		tb.setRV( r_osc, v_osc);
		
	}
	
    /** Integrate from t0 to tf.
     * @param t0    initial time or independent variable
     * @param x0    double[] containing the initial state.
     * @param tf    final time or independent variable.
     * @param derivs   Object containing the Equations of Motion
     * @param pr    Object containing the print method.
     * @param print_switch  Boolean (true = print, false = don't print).
     * @return      double[] containing the new state
     */
    public VectorN runloop(double t0, VectorN rv0, double tf, Printable pr, boolean print_switch) {
        int neqns = 6;

        double dt = this.step_size;
        double t = t0;
        

        // print initial state
        if (print_switch) {
            pr.print(t, rv0.getArray());
        }
        VectorN dev = new VectorN(neqns);
//        dev = new VectorN(nystrom.step(t+dt, dev.getArray(), this));

        // main integration loop

        while (t < tf) {
            // update time
            if ((t + dt) > tf){
            	dt = tf - t;
            	this.step_size = dt;
            }
            t = t + dt;
        	
//        	System.out.println("loop "+t);
            // step osculating orbit forward
            updateOsculatingOrbit(dt);        	
        	// step the deviation
            dev = new VectorN(nystrom.step(t, dev.getArray(), this));
            VectorN tmp = dev.get(0, 3);
            devmag = tmp.mag();
            
            // compute new state
            rv = (r_osc.append(v_osc)).plus(dev);            

            // check for rectification condition
            VectorN rdev = dev.get(0, 3);
//            if ((rdev.mag()/r_osc.mag())> 0.01) {
//            if ((t > 432000.0)&&once) {
//            if (rectify) {
//            	rectify_counter++;
//            	rectify();
//            	dev = new VectorN(neqns);
//            	rectify = false;
//            	once = false;
//            }
            // print new state
            if (print_switch) {
            	pr.print(t, rv.getArray());
            }
//            System.out.print("t = "+t+" ");
           
        }
        
        return rv;
    }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Time t0 = new Time(2007, 7, 1, 12, 0, 0.0);
		VectorN r0 = new VectorN(6.678137E+06, 0.0, 0.0);
		VectorN v0 = new VectorN(0.0, 6.78953030027267E+03, 3.68641417440091E+03);
		VectorN rv0 = r0.append(v0);
		LinePrinter lp = new LinePrinter();
		lp.setIsAdaptive(true);
		Encke encke = new Encke(IERS_1996.GM_Earth, t0, r0, v0);
		VectorN rf = encke.runloop(0.0, rv0, 864000.0, lp, true);
		System.out.println("rectify counter = "+ encke.rectify_counter);
		System.out.println("dev mag = "+encke.devmag);
		
		lp.close();
	}

}
