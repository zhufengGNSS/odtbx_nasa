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

import jat.util.*;
import javax.media.j3d.*;

/**
 * Create Java3D trajectory from COPERNICUS data file
 *
 * @author Tobias Berthold
 */
public class Copernicus_Trajectory extends Shape3D
{
	int num = 0; // number of records in file
	public double[] coords;
	public double[] t_read, x_read, y_read, z_read; // values from file
	public double[] ux_read, uy_read, uz_read; // values from file
	public double[] t, x, y, z; // original or interpolated values
	public double[] ux, uy, uz; // original or interpolated values
	double tmp;

	public Copernicus_Trajectory(double[] coords)
	{
		this.coords = coords;
	}

	/**
	 * Create interpolated Java3D trajectory from COPERNICUS data file
	 *
	 * @param filename
	 * @param steps
	 */
	public Copernicus_Trajectory(String filename, int steps)
	{
		int i;

		get_data_from_file(filename); // read data into arrays t_read, etc.

		//int steps=10;
		double delta_time = (t_read[num - 1] - t_read[0]) / steps;
		jat.math.Interpolator in_x = new jat.math.Interpolator(t_read, x_read);
		jat.math.Interpolator in_y = new jat.math.Interpolator(t_read, y_read);
		jat.math.Interpolator in_z = new jat.math.Interpolator(t_read, z_read);
		jat.math.Interpolator in_ux = new jat.math.Interpolator(t_read, ux_read);
		jat.math.Interpolator in_uy = new jat.math.Interpolator(t_read, uy_read);
		jat.math.Interpolator in_uz = new jat.math.Interpolator(t_read, uz_read);
		// Create new arrays for interpolated values
		t = new double[steps + 1];
		x = new double[steps + 1];
		y = new double[steps + 1];
		z = new double[steps + 1];
		ux = new double[steps + 1];
		uy = new double[steps + 1];
		uz = new double[steps + 1];

		for (i = 0; i < steps + 1; i++)
		{
			t[i] = t_read[0] + i * delta_time;
			x[i] = in_x.get_value(t[i]);
			y[i] = in_y.get_value(t[i]);
			z[i] = in_z.get_value(t[i]);
			ux[i] = in_ux.get_value(t[i]);
			uy[i] = in_uy.get_value(t[i]);
			uz[i] = in_uz.get_value(t[i]);
			//System.out.println(t[i]+" "+x[i]);
		}

		// Copy data into coords array
		coords = new double[steps * 3];
		for (i = 0; i < steps; i++)
		{
			coords[i * 3 + 0] = x[i];
			coords[i * 3 + 1] = y[i];
			coords[i * 3 + 2] = z[i];
		}

		create_trajectory_lines();

	}

	/**
	 * Create Java3D trajectory from COPERNICUS data file
	 *
	 * @param filename
	 */
	public Copernicus_Trajectory(String filename)
	{

		get_data_from_file(filename); // read data into arrays t_read, etc.

		// Copy data into coords array
		//coords=new double[num*3];

		//draw_trajectory();
	}

	public Copernicus_Trajectory()
	{
		//draw_trajectory();
	}

	public void get_data_from_file(String filename)
	{
		EasyReader inFile = new EasyReader(filename);
		if (inFile.bad())
		{
			System.err.println("Can't open " + filename);
			System.exit(1);
		}

		// find out how many numbers in file
		while (!inFile.eof())
		{
			inFile.readLine();
			//System.out.println(inFile.readLine());
			num++;
		}
		inFile.close();
		num -= 3;
		System.out.println("lines " + num);

		// read the data from file
		inFile = new EasyReader(filename);
		t_read = new double[num];
		x_read = new double[num];
		y_read = new double[num];
		z_read = new double[num];
		ux_read = new double[num];
		uy_read = new double[num];
		uz_read = new double[num];

		int i, line = 0;
		inFile.readLine(); // read over first line

		for (i = 0; i < num; i++)
		{
			t_read[i] = inFile.readFORTRANDouble(); // Julian date
			tmp = inFile.readFORTRANDouble(); // sim_time
			tmp = inFile.readFORTRANDouble(); // seg_time
			x_read[i] = inFile.readFORTRANDouble(); // x(km)
			y_read[i] = inFile.readFORTRANDouble(); // y(km)
			z_read[i] = inFile.readFORTRANDouble(); // z(km)
			tmp = inFile.readFORTRANDouble(); // xd(km/s)
			tmp = inFile.readFORTRANDouble(); // yd(km/s)
			tmp = inFile.readFORTRANDouble(); // zd(km/s)
			tmp = inFile.readFORTRANDouble(); // mass(kg)
			ux_read[i] = inFile.readFORTRANDouble(); // t/m ux
			uy_read[i] = inFile.readFORTRANDouble(); // t/m uy
			uz_read[i] = inFile.readFORTRANDouble(); // t/m uz
			tmp = inFile.readFORTRANDouble(); // lx
			tmp = inFile.readFORTRANDouble(); // ly
			tmp = inFile.readFORTRANDouble(); // lz
			tmp = inFile.readFORTRANDouble(); // lxd
			tmp = inFile.readFORTRANDouble(); // lyd
			tmp = inFile.readFORTRANDouble(); // lzd
			tmp = inFile.readFORTRANDouble(); // lm
			tmp = inFile.readFORTRANDouble(); // lvmag
			tmp = inFile.readFORTRANDouble(); // lvdmag
			tmp = inFile.readFORTRANDouble(); // lrvlvd
			tmp = inFile.readFORTRANDouble(); // alpha(deg)
			tmp = inFile.readFORTRANDouble(); // beta(deg)
			tmp = inFile.readFORTRANDouble(); // u_dot_mag(deg/day)
			tmp = inFile.readFORTRANDouble(); // Hamiltonian
			tmp = inFile.readFORTRANDouble(); // switch_func
			tmp = inFile.readFORTRANDouble(); // thrust(N)
			tmp = inFile.readFORTRANDouble(); // Isp(s)
			tmp = inFile.readFORTRANDouble(); // c_exhaust(km/s)
			tmp = inFile.readFORTRANDouble(); // power(watts)
			//System.out.println(""+coords[line+2]);
			inFile.readLine();
			line += 3;
		}
		inFile.close();
	}

	private void create_trajectory_lines()
	{

		int num_vert = coords.length / 3;
		//int[] stripLengths = { 200};
		int[] stripLengths = { num_vert };

		LineStripArray myLines = new LineStripArray(num_vert, GeometryArray.COORDINATES, stripLengths);
		myLines.setCoordinates(0, coords);

		this.setGeometry(myLines);
	}

	private void create_test_lines()
	{
		double[] coords = { 0., 0., 0., 100000., 0., 0., 1000., 1000., 0., 1000., 1000., 1000. };

		int[] stripLengths = { 4 };

		LineStripArray myLines = new LineStripArray(4, GeometryArray.COORDINATES, stripLengths);
		myLines.setCoordinates(0, coords);

		this.setGeometry(myLines);
	}
}
