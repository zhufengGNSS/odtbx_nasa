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
import jat.util.*;
import jat.matvec.data.*;

/** Base class for 3D objects. Planet, moons, and spacecraft extend this class.
 * To be extended by specific object
 * @author Tobias Berthold
 */
// Temporarily renamed class Body3D to ABody3D because of difficulties with CVS

//abstract public class body3D extends TransformGroup
public class Body3D extends TransformGroup
{
	Vector3f Vf = new Vector3f();
	Vector3d Vd = new Vector3d();
	Vector3d VRot = new Vector3d();
	Transform3D Trans = new Transform3D();
	double scale = 1.0; // scale factor for 3D objects
	Applet myapplet;
	static String images_path, Lightwave_path, Wavefront_path, ThreeDStudio_path;

	public Body3D()
	{
		setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
	}

	public Body3D(Applet myapplet)
	{
		this.myapplet = myapplet;
		String thisClassName="Body3D";
		images_path = FileUtil.getClassFilePath("jat.vr", thisClassName) + "images_hires/";
		Wavefront_path = FileUtil.getClassFilePath("jat.vr", thisClassName) + "Wavefront\\";
		Lightwave_path = FileUtil.getClassFilePath("jat.vr", thisClassName) + "Lightwave\\";
		ThreeDStudio_path = FileUtil.getClassFilePath("jat.vr", thisClassName) + "3DStudio\\";
		//System.out.println(images_path);

		setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		setCapability(TransformGroup.ALLOW_TRANSFORM_READ);

		//add();
	}

	// Methods to be implemented in subclasses
	//abstract public void add();

	public void set_scale(double scale)
	{
		getTransform(Trans);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	/** Set the position of the body in km
	 * @param x x position
	 * @param y y position
	 * @param z z position
	 */
	public void set_position(double x, double y, double z)
	{
		getTransform(Trans);
		Vd.x = x;
		Vd.y = y;
		Vd.z = z;
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	/**
	 * Set new body position.
	 * @param rv new position
	 */
	public void set_position(double[] rv)
	{
		getTransform(Trans);
		Vd.x = rv[0];
		Vd.y = rv[1];
		Vd.z = rv[2];
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	/**
	 * Set new body position.
	 * @param r new position
	 */
	public void set_position(VectorN r)
	{
		getTransform(Trans);
		Vd.x = r.x[0];
		Vd.y = r.x[1];
		Vd.z = r.x[2];
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	public void set_position(Point3d r)
	{
		getTransform(Trans);
		Vd.x = r.x;
		Vd.y = r.y;
		Vd.z = r.z;
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	public void set_position(Vector3d r)
	{
		getTransform(Trans);
		Vd.x = r.x;
		Vd.y = r.y;
		Vd.z = r.z;
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	/**
	 * Set_body attitude.using quaternion
	 * @param quatObject quaternion
	 */
	public void set_attitude(Transform3D quatObject)
	{
		setTransform(quatObject);
	}

	/** Set body attitude without changing position or scale using Euler angles
	 * @param alpha x angle
	 * @param beta y angle
	 * @param gamma z angle
	 */
	public void set_attitude(double alpha, double beta, double gamma)
	{
		getTransform(Trans);
		Trans.get(Vf);
		VRot.x = alpha;
		VRot.y = beta;
		VRot.z = gamma;
		Trans.setEuler(VRot);
		Trans.setTranslation(Vf);
		Trans.setScale(scale);
		setTransform(Trans);
	}

	/**
	 * Set body position and attitude
	 * @param x
	 * @param y
	 * @param z
	 * @param alpha
	 * @param beta
	 * @param gamma
	 */
	public void set_pos_attitude(double x, double y, double z, double alpha, double beta, double gamma)
	{
		getTransform(Trans);
		Trans.get(Vf);
		VRot.x = alpha;
		VRot.y = beta;
		VRot.z = gamma;
		Trans.setEuler(VRot);
		Vd.x = x;
		Vd.y = y;
		Vd.z = z;
		Trans.setTranslation(Vd);
		Trans.setScale(scale);
		setTransform(Trans);
	}
}


//Transform3D T_3D = new Transform3D();
//Transform3D RotX = new Transform3D();
//Transform3D RotY = new Transform3D();
//Transform3D RotZ = new Transform3D();
	/*
		public void set_attitude(double alpha, double beta, double gamma)
		{
			getTransform(T_3D);
			T_3D.get(Vf);
	//		translate.set(Vf);

			RotX.setIdentity();
			RotY.setIdentity();
			RotZ.setIdentity();
			RotX.rotX(alpha);
			RotY.rotY(beta);
			RotZ.rotZ(gamma);
			RotX.mul(RotY); // RotX=RotX . RotY
			RotX.mul(RotZ); // RotX=RotX . RotZ

			RotX.setTranslation(Vf);
			setTransform(RotX);
		}


		public void set_attitude(double alpha, double beta, double gamma)
		{
			getTransform(RotX);
			RotX.get(Vf);
	//		scale=RotX.getScale();
			RotY.setIdentity();
			RotZ.setIdentity();
			RotX.rotX(alpha);
			RotY.rotY(beta);
			RotZ.rotZ(gamma);
			RotX.mul(RotY); // RotX=RotX . RotY
			RotX.mul(RotZ); // RotX=RotX . RotZ

			RotX.setTranslation(Vf);
			RotX.setScale(scale);
			setTransform(RotX);
		}
	*/


//	public void set_earth_rotation(double angle)
//	{
//		earthRotate.setIdentity();
//		earthRotate.rotZ(angle);
//	}
//
//	public void rotate_earth()
//	{
//		TG_earth.getTransform(T_3D);
//		T_3D.mul(earthRotate, T_3D);
//		TG_earth.setTransform(T_3D);
//	}
