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
import jat.cm.*;
import jat.matvec.data.VectorN;

/**
 * The jat_view class contains methods related to the viewing platform in Java3D,
 * specific to JAT.
 *
 * @author Tobias Berthold
 * @version 1.0
*  Date        :   10-7-2002
 */
public class jat_view
{
	public static TransformGroup TG_vp;
//	static Matrix3d m1 = new Matrix3d();
//	static Matrix3d m2 = new Matrix3d();
	static Transform3D Trans = new Transform3D();
	public static View view;
	public static float FrontClipDistance = 1.e3f;
	public static float BackClipDistance = 1.e9f;
	public static Vector3d V_initial_view_location = new Vector3d(-50000.0f, 0.0f, 0.0f);
	public static Vector3d V = new Vector3d();
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);

	public static BranchGroup view(double x, double y, double z, Canvas3D c)
	{
		V_initial_view_location.x = x;
		V_initial_view_location.y = y;
		V_initial_view_location.z = z;
		return view(c);
	}

	public static BranchGroup view(double x, double y, double z, Canvas3D c, float f, float b)
	{
		V_initial_view_location.x = x;
		V_initial_view_location.y = y;
		V_initial_view_location.z = z;
		FrontClipDistance = f;
		BackClipDistance = b;
		return view(c);
	}

	public static BranchGroup view(Canvas3D c)
	{
		// View
		PhysicalBody body = new PhysicalBody();
		PhysicalEnvironment environment = new PhysicalEnvironment();
		view = new View();
		view.addCanvas3D(c);
		view.setPhysicalBody(body);
		view.setPhysicalEnvironment(environment);
		view.setProjectionPolicy(View.PERSPECTIVE_PROJECTION);
		view.setFrontClipDistance(FrontClipDistance);
		view.setBackClipDistance(BackClipDistance);

		//System.out.println(" "+Float.MAX_VALUE);
		//view.setFieldOfView(cm.Rad(45.));
		view.setFieldOfView(cm.Rad(30.));

		//Transform3D t = new Transform3D();
		//t.set(new Vector3f(0.0f,0.0f,0.0f));
		//Vector3d view_location=new Vector3d(-50000.0f,0.0f,0.0f);
		Vector3d look_to = new Vector3d(0., 0., 0.);
		Transform3D t = jat_view.lookat_T(V_initial_view_location, look_to);
		//Transform3D t = new Transform3D();
		ViewPlatform vp = new ViewPlatform();
		TG_vp = new TransformGroup(t);
		TG_vp.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		TG_vp.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		TG_vp.addChild(vp);
		BranchGroup BG_vp = new BranchGroup();
		BG_vp.addChild(TG_vp);
		view.attachViewPlatform(vp);

		return BG_vp;
	}

	public static void set_FrontBackClipDistance(double FrontClipDistance, double BackClipDistance)
	{
		jat_view.FrontClipDistance = (float)FrontClipDistance;
		jat_view.BackClipDistance = (float)BackClipDistance;
		view.setFrontClipDistance(FrontClipDistance);
		view.setBackClipDistance(BackClipDistance);
	}

	/**
	 * Set the coordinates of the viewers position without changing the viewing direction
	 *
	 * @param position
	 */
	public static void set_view_position(Vector3d position)
	{
		TG_vp.getTransform(Trans);
		Trans.setTranslation(position);
		TG_vp.setTransform(Trans);
	}

	/**
	 * Set the coordinates of the viewers position without changing the viewing direction
	 *
	 * @param position
	 */
	public static void set_view_position(VectorN position)
	{
		set_view_position(position.get(0), position.get(1), position.get(2));
	}

	/**
	 * Set the coordinates of the viewers position without changing the viewing direction
	 *
	 * @param
	 */
	public static void set_view_position(float x, float y, float z)
	{
		TG_vp.getTransform(Trans);
		V.x = x;
		V.y = y;
		V.z = z;
		Trans.setTranslation(V);
		TG_vp.setTransform(Trans);
	}

	/**
	 * Set the coordinates of the viewers position without changing the viewing direction
	 *
	 * @param
	 */
	public static void set_view_position(double x, double y, double z)
	{
		TG_vp.getTransform(Trans);
		V.x = x;
		V.y = y;
		V.z = z;
		Trans.setTranslation(V);
		TG_vp.setTransform(Trans);
	}

	public static void set_view_rotation(Vector3d rotation)
	{	
		TG_vp.getTransform(Trans);
		Trans.get(V);	// save translational component
		Trans.setEuler(rotation);
		Trans.setTranslation(V);	// restore translational component
		TG_vp.setTransform(Trans);
	}

	public static void set_view_direction(double x, double y, double z)
	{
		TG_vp.getTransform(Trans);
		Trans.get(V);
		Vector3d look_to = new Vector3d(x, y, z);
		Trans = jat_view.lookat_T(V, look_to);
		TG_vp.setTransform(Trans);
	}

	/** Change view direction while keeping z-axis up
	 * @param look_to
	 */
	public static void set_view_direction(Vector3d position, Vector3d look_to)
	{
		Trans = jat_view.lookat_T(position, look_to);
		TG_vp.setTransform(Trans);
	}

	public static void set_pos_direction(Vector3d position, Vector3d look_to)
	{
		Trans = jat_view.lookat_T(position, look_to);
		Trans.setTranslation(position);
		TG_vp.setTransform(Trans);
	}


	public static Transform3D lookat_T(Vector3d position, Vector3d look_to)
	{
		Vector3d diff=new Vector3d();
		diff.sub(look_to,position);
		double x=diff.x;
		double y=diff.y;
		double z=diff.z;
//		System.out.println("Look to : x "+x+"  y "+y+"  z "+z);
		Vector3d VRot = new Vector3d();
		Vector3d zaxis=new Vector3d(0.,0.,-1.);
		Vector3d project_yz=new Vector3d(0.,Math.sqrt(y*y+z*z),z);
//		System.out.println("Projection : x "+project_yz.x+"  y "+project_yz.y+"  z "+project_yz.z);
		// Angle between negative z-axis and projection of vector onto y-z-plane
		//double xrot=Math.atan2(z,Math.sqrt(y*y+z*z));
		double xrot=zaxis.angle(diff);
//		double xrot=zaxis.angle(project_yz);
		//double xrot_angle=project_yz.angle(zaxis);
		// Angle between vector and x-z-plane	
		//double zrot=zaxis.angle(project_yz); 
		double zrot=-Math.atan2(x,y); // Angle between vector and x-z-plane 
//		System.out.println("xrot angle "+Constants.rad2deg*xrot+" Degrees");
//		System.out.println("zrot angle "+Constants.rad2deg*zrot+" Degrees");
		VRot.x = xrot;
		VRot.y = 0;
		VRot.z = zrot;
		Trans.setEuler(VRot);
		position.get(V);			// save translational component
		Trans.setTranslation(V);	// restore translational component
		return Trans;
	}




	public static Vector3d get_view()
	{
		TG_vp.getTransform(Trans);
		Trans.get(V);
		return V;
	}




	/**
	 * Set the view direction without changing the viewers position
	 *
	 * @param position First point for constructing a viewing direction
	 * @param look_to Second point to which to set the viewing direction
	 */
/*
	public static void set_view_direction(Vector3d position, Point3d look_to)
	{
		Trans = jat_view.lookat_T(position, look_to);
		TG_vp.setTransform(Trans);
	}

	public static Vector3d get_vp()
	{
		TG_vp.getTransform(Trans);
		Trans.get(V);
		return V;
	}
*/
	/**
	 * Method set_view_direction.
	 * Set the view direction without changing the viewers position
	 *
	 * @param look_to point to which to set the viewing direction
	 */
/*	
	public static void set_view_direction(Point3d look_to)
	{
		//TG_vp.getTransform().;
		//Vector3d position=TG_vp.getTransform().;
		Trans = jat_view.lookat_T(get_vp(), look_to);
		TG_vp.setTransform(Trans);
	}

	public static void set_view_direction(double x, double y, double z)
	{
		Point3d look_to = new Point3d(x, y, z);
		Trans = jat_view.lookat_T(get_vp(), look_to);
		TG_vp.setTransform(Trans);
	}
*/

	/**
	 * Method lookat_T.
	 * @param position
	 * @param look_to
	 * @return Transform3D
	 */
/*
	public static Transform3D lookat_T(Vector3d position, Point3d look_to)
	{

		Matrix3d m1 = mymat(position, look_to);
		Transform3D T_3d = new Transform3D(m1, position, 1.0);
		return T_3d;
	}
	private static Matrix3d mymat(Vector3d location, Point3d look_to)
	{
		Vector3d negzaxis = new Vector3d(0.0, 0.0, -1.0);
		Vector3d yaxis = new Vector3d(0.0, 1.0, 0.0);
		Vector3d lk_vector = new Vector3d();
		Vector3d tmp = new Vector3d();
		double lk_angle;
		double tilt;

		// Determin the angle and normal from the normal look angle (straight
		// down the Z axis, to the look_to point
		//
		lk_vector.sub(look_to, location);
		//System.out.println("Vector " + lk_vector);

		// Get
		tmp.set(lk_vector.x, 0.0, lk_vector.z);
		lk_angle = negzaxis.angle(tmp);
		if (tmp.x < 0)
			lk_angle = -lk_angle;
		m1.setIdentity();
		if (!Double.isNaN(lk_angle))
			m1.rotY(lk_angle);
		m1.transform(lk_vector);
		tilt = lk_vector.angle(negzaxis);
		if (lk_vector.y < 0)
			tilt = -tilt;

		m1.setIdentity();
		if (!Double.isNaN(lk_angle))
			m1.rotY(-lk_angle);
		m2.rotX(tilt);
		m1.mul(m2);

		return m1;
	}
*/	
	
}
