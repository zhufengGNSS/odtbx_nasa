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

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.loaders.Scene;
import com.mnstarfire.loaders3d.Inspector3DS;

/** 3DStudio Object class
 * @author Tobias Berthold
 */
public class ThreeDStudioObject extends Body3D
{
	Sphere3D left, right;
	Cylinder3D lc,rc;

	public ThreeDStudioObject(Applet myapplet, String filename, float scale)
	{
		super(myapplet);
		super.scale = scale;
		String fullname;
		Scene s = null;

		// Construct the Lw3d loader and load the file
		fullname = ThreeDStudio_path + filename;
		System.out.println(fullname);

		Inspector3DS loader = new Inspector3DS(fullname); // constructor
		loader.parseIt(); // process the file
		TransformGroup theModel = loader.getModel();

		//Orient vehicle such that its main axis is lined up with x-axis
		Transform3D Rotate = new Transform3D();
		Transform3D R1 = new Transform3D();
		Transform3D R2 = new Transform3D();
		Rotate.setIdentity();
		R1.setIdentity();
		R2.setIdentity();
		R1.rotZ(Math.PI / 2);
		R2.rotY(Math.PI / 2);
		R2.mul(R1); // R2 = R2 . R1
		Rotate = R2;

		// For MarsOrbit class
		//Rotate.rotY(0.6 * Math.PI);

		// These values compensate for the fact that carrier.3ds is not centered
		// at the center of mass. Different for every model- make parameter later
		Vector3d V = new Vector3d();
		Transform3D T_3D = new Transform3D();
		V.x = 875.f;
		V.y = -222.f;
		V.z = -209.f;
		T_3D.set(V);
		T_3D.mul(Rotate, T_3D);
		theModel.setTransform(T_3D);

		super.set_scale(scale);

		setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		setCapability(TransformGroup.ALLOW_TRANSFORM_READ);


		// Thrust plumes
		// values for sphere
//		Point3d lt = new Point3d(-380.f, -42.f, 0.f); //right thruster
//		Point3d rt = new Point3d(-380.f, 35.f, 0.f); //left thruster
		//addChild(left = new Sphere3D( 18, Colors.white,Colors.turquoise)); 
//		addChild(left  = new Sphere3D( 18, Colors.white,Colors.turquoise,lt.x, lt.y, lt.z)); 
//		addChild(right = new Sphere3D( 18, Colors.white,Colors.turquoise,rt.x, rt.y, rt.z));
//		left.set_scale(0);
//		right.set_scale(0);
		
 
		// values for cylinder
		Point3d rt = new Point3d(-485.f, -42.f, 0.f); //left thruster
		Point3d lt = new Point3d(-485.f, 35.f, 0.f); //right thruster
		addChild(rc = new Cylinder3D( 18,200, Colors.turquoise));
		addChild(lc = new Cylinder3D( 18,200, Colors.turquoise));
		rc.set_attitude(0,0,0.5*Math.PI);
		lc.set_attitude(0,0,0.5*Math.PI);
		rc.set_position(rt); 
		lc.set_position(lt); 
		rc.set_scale(0.01);
		lc.set_scale(0.01);

		addChild(theModel);

		//BoundingSphere bounds = new BoundingSphere(rt, 40.f);//100000000.0
		//addChild(jat_light.PointLight(bounds));
	}
	
	public void thrusters_on()
	{
		lc.set_scale(1.);
		rc.set_scale(1.);
	}
}
