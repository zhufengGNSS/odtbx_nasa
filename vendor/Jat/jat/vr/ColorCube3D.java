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

import javax.vecmath.Color3f;

import com.sun.j3d.utils.geometry.*;

/** Sphere class
 * @author Tobias Berthold
 */
public class ColorCube3D extends Body3D
{
	//float size=100.f;

	/**
	 * Constructor Sphere3D.
	 * @param size size
	 */
	public ColorCube3D(float size)
	{
//		this.size=size;
		addChild(new ColorCube(size));
	}

	public ColorCube3D(double size)
	{
//		this.size=size;
		addChild(new ColorCube(size));
	}

	public ColorCube3D(float size, Color3f color, float x, float y, float z)
	{
		addChild(new ColorCube(size));
		set_position(x, y, z);
	}
}
