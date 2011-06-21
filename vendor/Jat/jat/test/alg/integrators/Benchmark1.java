/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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

package jat.test.alg.integrators;


import jat.alg.integrators.*;
import jat.cm.*;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;


/**
 * Benchmark and Validation of numerical integrators by propagating a two-body orbit 
 * over varying final times. Reference for validation is the Kepler solution for 
 * two-body orbits.
 * 
 * @author Tobias Berthold
 *
 */
public class Benchmark1 implements Printable, Derivatives
{
	protected double mu = 398600.4415; // GM in km^3/s^2 (value from JGM-3	
	static int dcount=0;
	
	public void print(double t, double [] y)
    {

        // handle the first variable for plotting - this is a little mystery but it works
        boolean first = true;
        if (t == 0.0) first = false;

        // also print to the screen for warm fuzzy feeling
        System.out.println(t+" "+y[0]+" "+y[1]);
    }

    
	public double[] derivs(double t, double[] y)
	{
		double out[] = new double[6];
		VectorN r = new VectorN(y[0], y[1], y[2]);
		double rmag = r.mag();
		double rcubed = rmag * rmag * rmag;
		double muorc = -1.0 * this.mu / rcubed;
	
		//System.out.println("I am being called! "+dcount++);
		dcount++;
		out[0] = y[3];
		out[1] = y[4];
		out[2] = y[5];
		out[3] = muorc * y[0];
		out[4] = muorc * y[1];
		out[5] = muorc * y[2];

		return out;
	}    
    
	
	public static void main(String[] args)
	{
		double[] x0, xf;
		Benchmark1 dpt2 = new Benchmark1();

        // create a TwoBody orbit using orbit elements
        TwoBody orbit = new TwoBody(7000.0, 0.5, 0.0, 0.0, 0.0, 0.0);
        orbit.setSteps(2);

        // initial time and final time, initial state
        double t0 = 0.0;
        double period=orbit.period();
        double tf = 9.*period;
        x0 = orbit.randv();
		new VectorN(x0).print("initial state");
        

        // create an Dormand Prince integrator
        DormandPrince dp=new DormandPrince();
        dp.setIntegration_error(1.e-7);

        // create an RungeKuttaFehlberg78 integrator
        RungeKuttaFehlberg78 rk78f = new RungeKuttaFehlberg78();
        rk78f.setAccuracy(1.e-14);

        
        // warm-up: propagate the orbit with each method
        // Kepler solution
        orbit.propagate(t0, tf);
        orbit.rv.print("Kepler equation");

        // Dormand Prince
        dcount=0;
        xf = dp.integrate(t0, x0, tf, dpt2);  
		new VectorN(xf).print("DP");
      
		// Runge-Kutta
		dcount=0;
        xf = rk78f.integrate(t0, x0, tf, dpt2, dpt2, false);
		new VectorN(xf).print("RKF78");        

		// create the table
		Matrix M=new Matrix(10,6);
		for(int i=0;i<10;i+=1)
		{
			int j=2*i;
	        tf = j*period;
			M.set(i,0 , j);
	        // Kepler solution
			orbit.propagate(t0, tf);
			M.set(i,1 , orbit.rv.x[1]);
	        // Dormand Prince
	        dcount=0;
	        xf = dp.integrate(t0, x0, tf, dpt2);  
			M.set(i,2 , xf[1]);
			M.set(i,3 , dcount);
			// Runge-Kutta
	        dcount=0;
	        xf = rk78f.integrate(t0, x0, tf, dpt2, dpt2, false);
			M.set(i,4 , xf[1]);
			M.set(i,5 , dcount);
		}	
		//new VectorN(x0).print("initial state");
		M.print("periods  Kepler DP feval  RK  feval");
	}
}
