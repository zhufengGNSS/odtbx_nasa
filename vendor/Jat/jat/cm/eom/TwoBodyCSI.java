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

package jat.cm.eom;

import jat.alg.*;
import jat.alg.integrators.*;
import jat.cm.*;
import jat.matvec.data.*;

/**
 * Equations of motion for a powered spacecraft with constant specific impulse
 * in the force field of the two-body problem  
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public class TwoBodyCSI extends TwoBody implements Derivatives /*, Printable*/
{
	int n = 7;
	double T; // Thrust magnitude (parameter but always max when on)
	double Isp; // Specific Impulse
	double c; // Exhaust velocity
	double g0 = 0.00981;
	//double b=-T/c;
	VectorTimeFunction u; // thrust direction
	public VectorN rvm; // vector containing position, velocity, and mass

	/** Create the equations of motion for a spacecraft with a constant specific 
	 * impulse engine in a two-body force field. The initial state of the spacecraft 
	 * is given by a state on a two-body orbit. The thrust direction u is a callback
	 * function that the user must implement. 
	 * @param k Orbital elements of the initial state
	 * @param mu 
	 * @param m	initial mass in kg
	 * @param u thrust direction at a given time
	 * @param T maximum thrust magnitude in kg*m/s^2
	 * @param Isp Specific impulse in seconds
	 */
	public TwoBodyCSI(KeplerElements k, double mu, double m, VectorTimeFunction u, double T, double Isp)
	{
		super(mu, k);
		this.u = u;
		this.Isp = Isp;
		this.c = Isp * g0;
		this.T = T;
		rvm = new VectorN(7);
		rvm.x[0] = rv.x[0];
		rvm.x[1] = rv.x[1];
		rvm.x[2] = rv.x[2];
		rvm.x[3] = rv.x[3];
		rvm.x[4] = rv.x[4];
		rvm.x[5] = rv.x[5];
		rvm.x[6] = m;
	}

	public double[] derivs(double t, double[] x_)
	{
		double x = x_[0], y = x_[1], z = x_[2];
		double vx = x_[3], vy = x_[4], vz = x_[5];
		double m = x_[6];
		double dxdt[] = new double[n];
		double ux, uy, uz; // components of thrust direction vector
		double r = Math.sqrt(x * x + y * y + z * z);
		double rcubed = r * r * r;
		double muorc = -1.0 * this.mu / rcubed;

		// Terms
		//r = Math.sqrt(x * x + y * y + z * z);
		//r3 = r * r * r;
		//T = Math.sqrt(Tx * Tx + Ty * Ty + Tz * Tz);

		// get thrust direction for given state and time
		VectorN nxf = new VectorN(x_);
		VectorN cf = u.evaluate(nxf, t);
		//cf.print();
		ux = cf.x[0];
		uy = cf.x[1];
		uz = cf.x[2];

		// Derivatives
		dxdt[0] = vx;
		dxdt[1] = vy;
		dxdt[2] = vz;
		dxdt[3] = muorc * x + T/m * ux ;
		dxdt[4] = muorc * y + T/m * uy ;
		dxdt[5] = muorc * z + T/m * uz ;
		dxdt[6] = -T / c;

		return dxdt;
	}
}

/*
public double[] integrate_TCT(double[] x0, double t_0, double t_s1, double t_s2, double t,  double T)
{
	double[] x;

	// Thrust arc
	if ( t<=t_s1 )
	{
		x=integrate(x0, t_0, t, T);
	}
	else
	{
		x=integrate(x0, t_0, t_s1, T);
		// Coast arc
		if ( t<=t_s2 )
			x=integrate(x, t_s1, t, 0.);
		else
		{
			x=integrate(x, t_s1, t_s2, 0.);
			// Thrust arc
			x=integrate(x, t_s2, t, T);
		}
	}
	//myutil.print_Matrix(new Matrix(x,1).transpose(),"state x");
	return x;
}


double[] integrate(double[] x0, double t_start, double t_final,  double T)
{
	this.T=T;
	return integrate(x0, t_start, t_final);
}

double[] integrate(double[] x0, double t_start, double t_f)
{
	double[] x,tmp_x0=new double[n];

	myutil.copy(x0, tmp_x0, n);
	x=super.integrate(tmp_x0, t_start, t_f);
	return x;
}

public double switching_function(double[] x0, double t_0, double t_s1, double t_s2, double t_s, double T)
{
	double[] x_=integrate_TCT(x0, t_0, t_s1, t_s2, t_s, T);
	double x=x_[0],y=x_[1],z=x_[2],vx=x_[3],vy=x_[4],vz=x_[5],m=x_[6];
	double lx=x_[7],ly=x_[8],lz=x_[9],lvx=x_[10],lvy=x_[11],lvz=x_[12],lm=x_[13];
	double lv=Math.sqrt(lvx*lvx+lvy*lvy+lvz*lvz);
	double c=Isp*g0;

	return lv/m-lm/c;
}
*/
