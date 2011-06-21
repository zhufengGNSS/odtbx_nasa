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

import javax.media.j3d.*;
import javax.vecmath.*;

public class Surface3D extends Shape3D
{

	////////////////////////////////////////////
	//
	// create Shape3D with geometry and appearance
	// the geometry is created in method yoyoGeometry
	// the appearance is created in method yoyoAppearance
	//
	public Surface3D()
	{

		this.setGeometry(yoyoGeometry());
		this.setAppearance(yoyoAppearance());

	} // end of Yoyo constructor

	////////////////////////////////////////////
	//
	// create yoyo geometry
	// four triangle fans represent the yoyo
	// strip   indicies_______________
	//   0     0N+0 to 1N+0 ( 0 to N )
	//   1     1N+1 to 2N+1
	//   2     2N+2 to 3N+2
	//   3     3N+4 to 4N+3
	//
	private Geometry yoyoGeometry()
	{

		TriangleFanArray tfa;
		int N = 17;
		int totalN = 4 * (N + 1);
		Point3f coords[] = new Point3f[totalN];
		Color3f colors[] = new Color3f[totalN];
		Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		Color3f yellow = new Color3f(0.7f, 0.5f, 0.0f);
		int stripCounts[] = { N + 1, N + 1, N + 1, N + 1 };
		float r = 6000.f;
		float w = 4000.f;
		int n;
		double a;
		float x, y;

		// set the central points for four triangle fan strips
		coords[0 * (N + 1)] = new Point3f(0.0f, 0.0f, w);
		coords[1 * (N + 1)] = new Point3f(0.0f, 0.0f, 0.0f);
		coords[2 * (N + 1)] = new Point3f(0.0f, 0.0f, 0.0f);
		coords[3 * (N + 1)] = new Point3f(0.0f, 0.0f, -w);

		colors[0 * (N + 1)] = red;
		colors[1 * (N + 1)] = yellow;
		colors[2 * (N + 1)] = yellow;
		colors[3 * (N + 1)] = red;

		for (a = 0, n = 0; n < N; a = 2.0 * Math.PI / (N - 1) * ++n)
		{
			x = (float) (r * Math.cos(a));
			y = (float) (r * Math.sin(a));
			coords[0 * (N + 1) + n + 1] = new Point3f(x, y, w);
			coords[1 * (N + 1) + N - n] = new Point3f(x, y, w);
			coords[2 * (N + 1) + n + 1] = new Point3f(x, y, -w);
			coords[3 * (N + 1) + N - n] = new Point3f(x, y, -w);

			colors[0 * (N + 1) + N - n] = red;
			colors[1 * (N + 1) + n + 1] = yellow;
			colors[2 * (N + 1) + N - n] = yellow;
			colors[3 * (N + 1) + n + 1] = red;
		}

		tfa = new TriangleFanArray(totalN, TriangleFanArray.COORDINATES | TriangleFanArray.COLOR_3, stripCounts);

		tfa.setCoordinates(0, coords);
		tfa.setColors(0, colors);

		return tfa;

	}

	private Appearance yoyoAppearance()
	{

		Appearance appearance = new Appearance();

//		PolygonAttributes polyAttrib = new PolygonAttributes();
//		polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_LINE);
//		appearance.setPolygonAttributes(polyAttrib);

		TransparencyAttributes ta=new TransparencyAttributes();
		ta.setTransparency(0.5f);
		ta.setTransparencyMode(TransparencyAttributes.BLENDED);
		appearance.setTransparencyAttributes(ta);		

		return appearance;

	} // end of method yoyoAppearance of class Yoyo

} // end of class
