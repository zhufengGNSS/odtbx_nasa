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
 * File Created on Aug 28, 2003
 */
 
package jat.traj;

import jat.alg.integrators.*;
import jat.matvec.data.*;
import jat.cm.*;
import jat.timeRef.*;
import jat.forces.density.earth.*;
import jat.forces.gravity.earth.JGM3;
import jat.spacetime.Time;

/**
 * <P>
 * The ChaserBurnEOM Class implements the EquationsOfMotion interface. It provides
 * the equations of motion for a body in a JGM3 gravity field with Harris-Priester drag
 * and finite duration thrust.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
 
public class ChaserBurnEOM implements EquationsOfMotion {
	private Trajectory traj;
    private HarrisPriester hp; // earth atmosphere model
	private double t_mjd0;
	private double mass;
	private double area;
	private double cd;
	private JGM3 jgm3 = new JGM3(12,12);
	
	private VectorN fburn;
	private int index = 0;
	private boolean used = false;
	
	/**
	 * Constructor - see HarrisPriester for parameter details
	 * @param t Trajectory
	 * @param mjd0 initial MJD
	 * @param m mass
	 * @param a area
	 * @param c Cd
	 */
	public ChaserBurnEOM(Trajectory t, double mjd0, double m, double a, double c){
		this.traj = t;
		this.t_mjd0 = mjd0;
		this.mass = m;
		this.area = a;
		this.cd = c;
		this.hp = new HarrisPriester(this.cd, this.area, this.mass);
	}
	
	/**
	 * Set the next finite burn
	 * @param b FiniteBurn
	 */
	public void setBurn(FiniteBurn b) {
		double thrust = b.accel;
		VectorN unit = b.unitVector;
		this.fburn = unit.times(thrust);		
	}
	
	/** Implements the Printable interface to get the data out of the propagator and pass it to the trajectory.
	 *  This method is executed by the propagator at each integration step.
	 * @param t Time.
	 * @param y Data array.
	 */
	public void print(double t, double[] y) {
		traj.add(t, y);
	}
	
	/** Computes the derivatives for integration.
	 * @param t Time (not used).
	 * @param y State array.
	 * @return Array containing the derivatives.
	 */
	public double[] derivs(double t, double[] y) {
		
		VectorN out = new VectorN(y.length);
		
		// strip out incoming data    
		VectorN r = new VectorN(y[0], y[1], y[2]);
		VectorN v = new VectorN(y[3], y[4], y[5]);
		Quaternion q = new Quaternion(y[6], y[7], y[8], y[9]);
		q.unitize();

//		// construct the true orbit
//		TwoBody true_orbit = new TwoBody(Constants.GM_Earth, r, v);
//
//		// compute local gravity
//		VectorN g = true_orbit.local_grav();
		
		double Mjd = this.t_mjd0 + t/86400.0;
        EarthRef ref = new EarthRef(Mjd);
        ref.setIERS(3.3E-07, 1.17E-06, 0.649232);
        Matrix E = ref.eci2ecef();

        // Acceleration due to harmonic gravity field        
        VectorN g = jgm3.gravity(r, E);
        
        // acceleration due to the burn
		RSW_Frame rsw = new RSW_Frame(r, v);
		Matrix Cb2i = rsw.ECI2RSW().transpose();
        
//        Matrix Cb2i = q.quat2DCM();
        VectorN thrust = Cb2i.times(this.fburn);
        
        Time time = new Time(Mjd);
        hp.compute(time,ref, r, v);
        VectorN drag = hp.dragAccel();
        VectorN sf = thrust.plus(drag);
        
        
        VectorN accel = sf.plus(g);
        
		
		// compute attitude rate
		
		VectorN w = rsw.omega();
		Matrix omega_true = Quaternion.omega(w);
		VectorN q_deriv = omega_true.times(q);

		// derivatives for true orbit
		out.set(0, v);
		out.set(3, accel);
		out.set(6, q_deriv);
		
		// accelerometer state
		Matrix eci2rsw = Cb2i.transpose();
		VectorN df = eci2rsw.times(sf);
		out.set(10, df);
		
//		System.out.println("df = "+df);
		
		// gyro state
		out.set(13, w);
		

		return out.x;
	}
	
	

}
