/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2004 United States Government as represented by the
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
 * Created on Jun 22, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package jat.demo.OST;

import jat.alg.VectorTimeFunction;
import jat.alg.integrators.Printable;
import jat.math.Interpolator;
import jat.matvec.data.VectorN;
import jat.plot.SinglePlot;
import java.text.*;

/**
 * @author Tobias Berthold
 *
 */
public class ThrustDirection implements Printable, VectorTimeFunction
{
	// warning: u must not be null vector
	double ux, uy, uz;
	double[] x_search;
	SinglePlot plot;
	boolean plot_yes_no;
	VectorN u = new VectorN(3);
	DecimalFormat df; // formatting for double  
	NumberFormat inf;  // formatting for int

	public ThrustDirection(SinglePlot plot)
	{
		this.plot = plot;
		df = (DecimalFormat)NumberFormat.getInstance();
		df.applyPattern(" ###.##;-###.##");
		df.setMinimumFractionDigits(2);
		df.setMinimumIntegerDigits(4);
		// Integers
		inf = NumberFormat.getInstance();
		inf.setMinimumIntegerDigits(6);

	}

	void set_x_search(double[] x_search)
	{
		this.x_search = x_search;
	}

	public void print(double t, double[] y)
	{
		if (plot_yes_no)
		{
			boolean first = true;
			if (t == 0.0)
				first = false;

			// print to the screen for warm fuzzy feeling
			//System.out.println(t + " " + y[0] + " " + y[1] + " " + y[6]);

			// add data point to the plot
			plot.plot.addPoint(1, y[0], y[1], first);
		}
	}

	/*
	void set_thrust(double theta)
	{
		ux = Math.cos(theta);
		uy = Math.sin(theta);
		uz = 0.;
		u.x[0] = Math.cos(x_search[1]);
		u.x[1] = Math.sin(x_search[1]);
	}
	*/
	void set_thrust(double t)
	{
		double theta;
		double tf;
		int n=x_search.length;
		int segs=n-1;
		
		tf = x_search[0];
		//int index=(int)Math.floor(segs*t/tf);
		int index=1;
		double r=t/tf;
		if (r <= 1/3.)
			index=1;
		if (r >= 1 / 3. && r <= 2 / 3.)
			index=2;
		if (r >= 2/3.)
			index=3;

		//System.out.println(df.format(t)+" "+df.format(tf)+" "+index);

		theta = x_search[index];
		//System.out.println(t);//	set_thrust(x_search[2]);
		u.x[0] = Math.cos(theta);
		u.x[1] = Math.sin(theta);
		u.x[2] = 0.;
	}

	public VectorN evaluate(VectorN x, double t)
	{
		set_thrust(t);

		return u;
		//System.out.println(tf.format(t)+"  "+df.format(ux.get_value(t))+" "+uy.get_value(t));
		//return u.unitVector();
	}

	public void thrust_interpolator(double tf_CSI)
	{
		double[] uxdata, uydata, uzdata, t;
		//		double[] uxdata={1.,1.,1.};
		//		double[] uydata={1.,1.,1.};
		//		double[] uzdata={1.,1.,1.};
		//		double[] t={0.,1.,1.};
		Interpolator ux, uy, uz;

		int segments = 10;
		uxdata = new double[segments + 1];
		uydata = new double[segments + 1];
		uzdata = new double[segments + 1];
		t = new double[segments + 1];
		for (int i = 0; i < segments + 1; i++)
		{
			t[i] = tf_CSI * i / segments;
			uxdata[i] = 0.1;
			uydata[i] = 0.1;
			uzdata[i] = 0.;
			//System.out.println(tf.format(t[i]) + " " + df.format(uxdata[i]) + " " + uydata[i] + " " + uzdata[i]);
		}
		ux = new Interpolator(t, uxdata);
		uy = new Interpolator(t, uydata);
		uz = new Interpolator(t, uzdata);
	}

}

//u.x[0] = ux;
//u.x[1] = uy;
//u.x[2] = uz;
/*
if (t <= tf / 3.)
	theta = x_search[1];
else
{
	if (t >= tf / 3. && t <= 2 * tf / 3.)
		theta = x_search[2];
	else
	{
		if (t <= 2 * tf / 3.)
			theta = x_search[3];
	}
}
*/
