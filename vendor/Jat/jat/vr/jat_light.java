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
/**
 * Lighting related functions for JAT
 * @author Tobias Berthold
 * @version 1.0
 */
public class jat_light
{
	public static DirectionalLight lightD;

	public static AmbientLight AmbientLight(BoundingSphere bounds)
	{
		AmbientLight lightA = new AmbientLight(Colors.white);
		lightA.setInfluencingBounds(bounds);
		//BG_root.addChild(lightA);
		Background background = new Background();
		//background.setApplicationBounds(bounds);
		//background.setColor(Colors.darkgray);
		background.setColor(Colors.blue);
		//BG_root.addChild(background);
		return lightA;
	}

	public static DirectionalLight DirectionalLight(BoundingSphere bounds)
	{
		// *** Lights
		// Set up the directional lights
		Color3f lightColor = new Color3f(.8f, .8f, .8f);
		Vector3f lightDirection = new Vector3f(10.0f, 10.f, 0.f);
		lightD = new DirectionalLight(lightColor, lightDirection);
		lightD.setInfluencingBounds(bounds);

		/*
		Vector3f direction = new Vector3f(1.0f, 0.f, 0.f);
		direction.normalize();
		lightD1.setDirection(direction);
		lightD1.setColor(new Color3f(1.0f, 1.0f, 1.0f));
		//BG_root.addChild(lightD1);
		 */

		Background background = new Background();
		//background.setApplicationBounds(bounds);
		background.setColor(0.1f, 0.1f, 0.1f);
		//BG_root.addChild(background);
		return lightD;
	}

	public static void setDirection(float x, float y, float z)
	{
		Vector3f direction = new Vector3f(x, y, z);
		lightD.setDirection(direction);
	}

	public SpotLight SpotLight(Bounds bounds, Point3f pos, float spread, float concentration)
	{
		SpotLight sl = new SpotLight();
		sl.setInfluencingBounds(bounds);
		sl.setPosition(pos);
		//sl.setSpreadAngle(spread);
		sl.setConcentration(concentration);
		return sl;
	}

	public static PointLight PointLight(Bounds bounds)
	{
		PointLight pl = new PointLight();
		pl.setInfluencingBounds(bounds);
		return pl;
	}

}
