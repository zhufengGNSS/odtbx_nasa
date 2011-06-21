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
package jat.vr;

import jat.cm.*;
import jat.alg.integrators.*;
import javax.media.j3d.*;
import javax.vecmath.*;

/**
 * The Orbit class creates a Shape3D object in Java3D for a twobody orbit. 
 *
 * @author Tobias Berthold
 * @version 1.0
 */
// TODO: orbit color, number of segments as parameters
public class Orbit extends Shape3D implements Printable
{
	public double[] coords;
	public double[] t, x, y, z;
	int j = 0;
	Color3f Color = Colors.pink;
	private double a; // sma in km
	private double e; // eccentricity
	private double i; // inclination in radians
	private double raan; // right ascension of ascending node in radians
	private double w; // argument of perigee in radians
	private double ta; // true anomaly in radians
	private int steps=500;
	Constants c = new Constants();

	public Orbit(double[] coords)
	{
		this.coords = coords;
	}

	public Orbit()
	{
		draw_orbit();
	}

	/** Construct a TwoBody orbit from 6 orbit elements. Angles are input in degrees.
	 * @param a Semi-major axis in km.
	 * @param e Eccentricity.
	 * @param i Inclination in degrees.
	 * @param raan RAAN in degrees.
	 * @param w Argument of perigee in degrees.
	 * @param ta True anomaly in degrees.
	 */
	public Orbit(double a, double e, double i, double raan, double w, double ta)
	{
		this.a = a;
		this.e = e;
		this.i = i;
		this.raan = raan;
		this.w = w;
		this.ta = ta;
		draw_orbit();
	}

	public Orbit(double a, double e, double i, double raan, double w, double ta, Color3f Color, int steps)
	{
		this.a = a;
		this.e = e;
		this.i = i;
		this.raan = raan;
		this.w = w;
		this.ta = ta;
		this.Color = Color;
		this.steps=steps;
		draw_orbit();
	}

	public Orbit(KeplerElements el)
	{
		this.a = el.a;
		this.e = el.e;
		this.i = Math.toRadians(el.i);
		this.raan = Math.toRadians(el.raan);
		this.w = Math.toRadians(el.w) ;
		this.ta = Math.toRadians(el.ta);
		draw_orbit();
	}

	public Orbit(KeplerElements el, Color3f Color)
	{
		this.a = el.a;
		this.e = el.e;
		this.i = Math.toRadians(el.i);
		this.raan = Math.toRadians(el.raan);
		this.w = Math.toRadians(el.w) ;
		this.ta = Math.toRadians(el.ta);
		this.Color = Color;
		draw_orbit();
	}

	public void print(double time, double[] pos)
	{
		// also print to the screen for warm fuzzy feeling
		//		System.out.println(j+"  "+time + " " + pos[0] + " " + pos[1] + " " + pos[2]);
		t[j] = time;
		x[j] = pos[0];
		y[j] = pos[1];
		z[j] = pos[2];
		//		coords[j + 0] = pos[0];
		//		coords[j + 1] = pos[1];
		//		coords[j + 2] = pos[2];
		j++;
	}

	private void draw_orbit()
	{
		//int steps = 10000;
		// number of steps in propagate in TwoBody should be made a parameter
		coords = new double[3 * steps + 6];

		t = new double[steps + 2];
		x = new double[steps + 2];
		y = new double[steps + 2];
		z = new double[steps + 2];

		// create a TwoBody orbit using orbit elements
		TwoBody sat = new TwoBody(a, e, i, raan, w, ta);
		//sat.printElements("Orbit");

		// find out the period of the orbit
		double tf = sat.period();

		// propagate the orbit
		//sat.propagate(0., tf, x, true);
		sat.propagate(0., tf, this, true,steps);

		// Copy data into coords array
		coords = new double[steps * 3];
		for (int k = 0; k < steps; k++)
		{
			coords[k * 3 + 0] = x[k];
			coords[k * 3 + 1] = y[k];
			coords[k * 3 + 2] = z[k];
		}
		int num_vert = coords.length / 3;
		int[] stripLengths = { num_vert };

		LineStripArray myLines =
			new LineStripArray(num_vert, GeometryArray.COORDINATES | GeometryArray.COLOR_3, stripLengths);
		Color3f colors[] = new Color3f[num_vert];
		for (int i = 0; i < num_vert; i++)
			colors[i] = Color;
		myLines.setColors(0, colors);
		myLines.setCoordinates(0, coords);

		this.setGeometry(myLines);
	}
}
