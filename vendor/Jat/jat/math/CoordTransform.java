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

package jat.math;

import javax.vecmath.Vector3d;

import jat.matvec.data.VectorN;

/**
 * Transformations between different coordinate systems
 * 
 * @author Tobias Berthold
 * 
 */
public class CoordTransform
{
	/**
	 * Convert spherical coordiantes to Cartesian coordinates
	 * @param r radius
	 * @param theta angle between z-axis and point in radians
	 * @param phi angle between x-axis and projection of point onto x-y plane in radians
	 * @return VectorN
	 */
	public static VectorN Spherical_to_Cartesian(double r, double theta, double phi)
	{
		VectorN out = new VectorN(3);
		out.set(0,r*Math.sin(theta)*Math.cos(phi));
		out.set(1,r*Math.sin(theta)*Math.sin(phi));
		out.set(2,r*Math.cos(theta));
		return out;
	}

	public static Vector3d Spherical_to_Cartesian(Vector3d coord)
	{
		Vector3d out = new Vector3d();
		double r=coord.x;
		double theta=coord.y;
		double phi=coord.z;
		out.x=r*Math.sin(theta)*Math.cos(phi);
		out.y=r*Math.sin(theta)*Math.sin(phi);
		out.z=r*Math.cos(theta);
		return out;
	}
}
