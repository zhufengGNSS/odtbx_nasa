/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.demo.OST;

import jat.alg.integrators.*;
import jat.cm.*;
import jat.plot.*;

/**
 * @author Tobias Berthold
 *
 */
public class Orbit implements Printable
{
	SinglePlot plot;
	TwoBody twobody;
	
	public Orbit(KeplerElements k_initial)
	{

		twobody = new TwoBody( cm.mu, k_initial);

	}

	public void plot_orbit(SinglePlot plot)
	{
		this.plot=plot;
		// create an RungeKuttaFehlberg78 integrator
		RungeKutta8 rk = new RungeKutta8();
		//RungeKuttaFehlberg78 rk=new RungeKuttaFehlberg78();
		//rk.setNonAdaptive();
		rk.setStepSize(3.e2);
		//rk.setAdaptive();
		double tfin = twobody.period();

		double[] x0initial = twobody.randv();
		rk.integrate(0., x0initial, tfin, twobody, this, true);
		
	}

	public void print(double t, double[] y)
	{
		boolean first = true;
		if (t == 0.0)
			first = false;

		// print to the screen for warm fuzzy feeling
		//System.out.println(t + " " + y[0] + " " + y[1] + " " + first);

		// add data point to the plot
		plot.plot.addPoint(0, y[0], y[1], first);
	}

}

//KeplerElements k_mid = new KeplerElements(11000.0, 0., 0.0, 0.0, 0.0, 0.0);
//Orbit mid=new Orbit(k_mid);
//mid.plot_orbit(plot);
