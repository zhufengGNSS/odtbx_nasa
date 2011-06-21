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
import jat.math.*;
import jat.matvec.data.VectorN;
import javax.media.j3d.*;
import javax.vecmath.*;

/**
 * @author Tobias Berthold
 *
 */
public class PlanetGrid3D extends Body3D
{
	double radius = 1000.;
	double long_section = 10;
	double lat_section = 10;
	int long_inc = 1;
	int lat_inc = 1;
	Color3f Color = Colors.darkgray;
	Constants c = new Constants();
	//CoordTransform t = new CoordTransform();

	/** Creates a 3D panet grid object
	 */
	public PlanetGrid3D(double radius)
	{
		super();
		this.radius = radius;
		createGeometry();
	}

	public PlanetGrid3D(double radius, double long_section, double lat_section)
	{
		super();
		this.radius = radius;
		this.long_section = long_section;
		this.lat_section = lat_section;
		createGeometry();
	}

	private void createGeometry()
	{
		VectorN point;
		int max;
		double theta = 0., phi = 0.;

		// Longitude circles
		max = 360 / long_inc + 1;
		//		int num_vert = coords.length / 3;
		int num_vert = max;
		int[] stripLengths = { num_vert };
		for (phi = 0.; phi < 360.; phi += long_section)
		{
			// one circle
			LineStripArray myLines = new LineStripArray(num_vert, GeometryArray.COORDINATES | GeometryArray.COLOR_3, stripLengths);
			double[] coords = new double[max * 3];
			theta = 0.;
			int j = 0;
			while (j < max * 3)
			{
				//System.out.println("theta " + theta + " phi " + phi);
				point = CoordTransform.Spherical_to_Cartesian(radius, Math.toRadians(theta), Math.toRadians(phi));
				//point.print();
				coords[j + 0] = point.get(0);
				coords[j + 1] = point.get(1);
				coords[j + 2] = point.get(2);
				theta += long_inc;
				j += 3;
			}
			myLines.setCoordinates(0, coords);
			// Set colors
			Color3f colors[] = new Color3f[num_vert];
			for (int i = 0; i < num_vert; i++)
				colors[i] = Color;
			myLines.setColors(0, colors);
			Shape3D s = new Shape3D();
			s.setGeometry(myLines);
			addChild(s);
		}

		// Latitude circles
		max = 360 / lat_inc + 1;
		for (theta = 0.; theta <= 181.; theta += lat_section)
		{
			// one circle
			LineStripArray myLines = new LineStripArray(num_vert, GeometryArray.COORDINATES | GeometryArray.COLOR_3, stripLengths);
			double[] coords = new double[max * 3];
			phi = 0.;
			int j = 0;
			while (j < max * 3)
			{
				//System.out.println("theta " + theta + " phi " + phi);
				point = CoordTransform.Spherical_to_Cartesian(radius, Math.toRadians(theta), Math.toRadians(phi));
				//point.print();
				coords[j + 0] = point.get(0);
				coords[j + 1] = point.get(1);
				coords[j + 2] = point.get(2);
				phi += lat_inc;
				j += 3;
			}
			myLines.setCoordinates(0, coords);
			// Set colors
			Color3f colors[] = new Color3f[num_vert];
			for (int i = 0; i < num_vert; i++)
				colors[i] = Color;
			myLines.setColors(0, colors);
			Shape3D s = new Shape3D();
			s.setGeometry(myLines);
			addChild(s);
		}
	} // end of createGeometry()

	/*
		private void createGeometrynew()
		{
			VectorN point;
			double theta = 0., phi = 0.;
	
			// Longitude circles
			for (phi = 0.; phi < 360.; phi += long_section)
			{
				// one circle
				addCircle(phi, theta);
			}
			// Latitude circles
			for (theta = 0.; theta <= 181.; theta += lat_section)
			{
				// one circle
				addCircle(phi, theta);
			}
		}
	
		private void addCircle(double phi, double theta)
		{
			VectorN point;
			int max;
	
			// smoothness of circle
			max = 360 / long_inc + 1;
			//		int num_vert = coords.length / 3;
			int num_vert = max;
			int[] stripLengths = { num_vert };
			// one circle
			LineStripArray myLines =
				new LineStripArray(num_vert, GeometryArray.COORDINATES | GeometryArray.COLOR_3, stripLengths);
			double[] coords = new double[max * 3];
			int j = 0;
			while (j < max * 3)
			{
				//System.out.println("theta " + theta + " phi " + phi);
				point = t.Spherical_to_Cartesian(radius, theta * c.deg2rad, phi * c.deg2rad);
				//point.print();
				coords[j + 0] = point.get(0);
				coords[j + 1] = point.get(1);
				coords[j + 2] = point.get(2);
				j += 3;
			}
			myLines.setCoordinates(0, coords);
			// Set colors
			Color3f colors[] = new Color3f[num_vert];
			for (int i = 0; i < num_vert; i++)
				colors[i] = Colors.darkgray;
			myLines.setColors(0, colors);
			Shape3D s = new Shape3D();
			s.setGeometry(myLines);
			addChild(s);
		}
	*/
} // end of class
