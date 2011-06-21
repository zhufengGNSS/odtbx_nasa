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

import java.applet.Applet;
import com.sun.j3d.loaders.lw3d.Lw3dLoader;
import com.sun.j3d.loaders.Loader;
import com.sun.j3d.loaders.Scene;

/** LightWaveObject class
 * @author Tobias Berthold
 */
public class LightWaveObject extends Body3D
{

	public LightWaveObject(Applet myapplet, String filename, float scale)
	{
		super(myapplet);
		this.scale = scale;
		String fullname;
		Scene s = null;

		// Construct the Lw3d loader and load the file
		fullname = Lightwave_path + filename;
		System.out.println(fullname);
		Loader lw3dLoader = new Lw3dLoader(Loader.LOAD_ALL);
		try
		{
			s = lw3dLoader.load(fullname);
		} catch (Exception e)
		{
			System.err.println("Exception loading file: " + e);
		}
		addChild(s.getSceneGroup());
	}
}
