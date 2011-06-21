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
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.Color3f;
import com.sun.j3d.utils.image.TextureLoader;
import com.sun.j3d.utils.geometry.*;

/** Planet class
 * @author Tobias Berthold
 */
public class Planet3D extends Body3D
{
	public static final int MERCURY = 1, VENUS = 2, EARTH = 3, MARS = 4, JUPITER = 5, MOON = 11;
	float radius;
	String Texturefilename;
	Appearance app;
	Color3f Planetcolor; // planet color if texture not found
	int divisions = 60; // number of divisions for sphere

	// TODO : make divisions a parameter

	/**
	 * @param myapplet Applet displaying this body (this parameter unfortunatley
	 * necessary for textures)
	 * @param planet_number
	 * @param scale
	 */
	public Planet3D(Applet myapplet, int planet_number, float scale)
	{
		super(myapplet);
		super.scale = scale;
		CreatePlanet(myapplet, planet_number);
	}

	public Planet3D(Applet myapplet, int planet_number)
	{
		super(myapplet);
		CreatePlanet(myapplet, planet_number);

	}

	// had to create this to have code that is common to different constructors.
	// Is there a better way?
	private void CreatePlanet(Applet myapplet, int planet_number)
	{

		switch (planet_number)
		{
			case MERCURY :
				Texturefilename = images_path + "mercury.jpg";
				radius = (float)cm.mercury_radius;
				Planetcolor = Colors.red;
				break;
			case VENUS :
				Texturefilename = images_path + "venus.jpg";
				radius = (float)cm.venus_radius;
				Planetcolor = Colors.green;
				break;
			case EARTH :
				Texturefilename = images_path + "earth.jpg";
				radius = (float)cm.earth_radius;
				Planetcolor = Colors.blue;
				break;
			case MARS :
				Texturefilename = images_path + "mars.jpg";
				radius = (float)cm.mars_radius;
				Planetcolor = Colors.blue;
				break;
			case JUPITER :
				Texturefilename = images_path + "jupiter.jpg";
				radius = (float)cm.jupiter_radius;
				Planetcolor = Colors.blue;
				break;
			case MOON :
				Texturefilename = images_path + "moon.jpg";
				radius = (float)cm.moon_radius;
				Planetcolor = Colors.blue;
				break;
		}

		if (Texturefilename == null)
		{
			app = createMatAppear_planet(Planetcolor, Colors.white, 10.0f);
		} else
		{
			TextureLoader tex = new TextureLoader(Texturefilename, myapplet);
			TextureAttributes ta = new TextureAttributes();
			ta.setTextureMode(TextureAttributes.MODULATE);
			app = createMatAppear_planet(Colors.white, Colors.white, 10.0f);
			app.setTextureAttributes(ta);
			app.setTexture(tex.getTexture());
		}

		//TG_plan.addChild( createLabel( szName, 1.2f, 1.2f, 0 ) );
		addChild(
			new Sphere(
				radius,
				Sphere.GENERATE_NORMALS | Sphere.GENERATE_TEXTURE_COORDS,
				divisions,
				app));
		set_scale(scale);				
	}

	protected static Appearance createMatAppear_planet(Color3f dColor, Color3f sColor, float shine)
	{
		Appearance appear = new Appearance();
		Material material = new Material();
		material.setDiffuseColor(dColor);
		material.setSpecularColor(sColor);
		material.setShininess(shine);
		material.setEmissiveColor(0.1f, 0.1f, 0.1f);
		appear.setMaterial(material);
		return appear;
	}
}
