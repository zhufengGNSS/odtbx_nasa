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

package jat.demo.vr.Copernicus1;

import jat.vr.*;
import jat.util.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; // For Timer
import java.awt.event.*;

/**
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public class QuaternionTest extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Transform3D Trans = new Transform3D();
	ThreeDStudioObject apollo;
	ColorCube3D cc;
	Planet3D earth, mars;
	RGBAxis3D mainaxis;
	Orbit apolloorbit = null;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	double alpha;
	int i = 0;

	public QuaternionTest()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.Copernicus1", "EarthOrbit");
		System.out.println(b);

		//Applet window
		setLayout(new BorderLayout());
		c = createCanvas(b + "frames/");
		add("Center", c);
		panel = new ControlPanel(BG_root);
		add("South", panel);

		// 3D Objects
		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);
		TG_scene = new TransformGroup();
		//TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH));
		TG_scene.addChild(mainaxis = new RGBAxis3D(500.0f));
		//TG_scene.addChild(cc = new ColorCube3D(100.0f));

		// Try Quaternions
		//		double theta=0.0 * Math.PI;
		//		Quat4d quat=new Quat4d(1*Math.cos(theta/2),0.,0., Math.sin(theta/2));
		//		Trans.setRotation(quat);
		//		cc.setTransform(Trans);
		//Transform3D Trans = new Transform3D();
		//Trans.setRotation(new Quat4d(0, 0, 1, 0.0 * Math.PI));
		//apollo.setTransform(Trans);

		TG_scene.addChild(apollo = new ThreeDStudioObject(this, "carrier.3DS", 1.f));
		apollo.set_position(0., 200., 0.);
		//apollo.set_position(cm.earth_radius + 1000., -2000., 1000.);

		// Try Euler angles
		apollo.set_attitude(0.,-0.25 * Math.PI, 0.5 * Math.PI);
		BG_root.addChild(TG_scene);

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(1.f, 1.f, -2.f);
		//BG_root.addChild(jat_light.AmbientLight(bounds));

		// View
		BG_vp = jat_view.view(1000., -1000., 2000., c, 1.f, 100000.f);
		//BG_vp = jat_view.view(0., 0., 1000., c, 1.f, 100000.f);

		//		Transform3D Trans = new Transform3D();
		//		Vector3f Vf = new Vector3f();
		//		jat_view.TG_vp.getTransform(Trans);
		//		Trans.get(Vf);
		//satTrans.setRotation(new Quat4d(1, 0, 0,Math.PI / 2. ));
		//Trans.setRotation(new Quat4d(1, 0, 0,0));
		//		jat_view.TG_vp.setTransform(Trans);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(10.f);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		// Use Timer for animation
		int delayValue = 50; // milliseconds
		animControl = new Timer(delayValue, this);
		//animControl.start();
	}

	// This method is called each time a timer event occurs
	public void actionPerformed(ActionEvent e)
	{
		//		if (i > 30 && i < 300)
		//		{
		//			animControl.stop();
		//			c.writeJPEG_ = true;
		//			animControl.start();
		//			//a.c.repaint();
		//		}
		//System.out.println("" + i);
		i++;
		if (i > 50)
		{
			apollo.thrusters_on();
			apollo.set_position(i - 50, 0, 0);
		}

		alpha += 0.002;

		// Try Euler angles
		//apollo.set_attitude(alpha,alpha,0);
		//apollo.set_attitude(0.,0.,alpha * Math.PI);

		// Try Quaternions
		//		alpha += 0.1 * Math.PI;
		//		System.out.println("" + alpha);
		//		Transform3D satTrans = new Transform3D();
		//		satTrans.setRotation(new Quat4d(0, 0, 1, alpha));
		//		apollo.setTransform(satTrans);
		//apollo.set_attitude(satTrans);

		//		mars.set_attitude(Math.PI / 2., alpha, 0);

		Vector3d V_view = new Vector3d(0.f, 0.f, 0.0f);
		Transform3D T3D = new Transform3D();
		jat_view.TG_vp.getTransform(T3D);
		T3D.get(V_view);
		panel.label.setText("  x " + (long)V_view.x + "  y " + (long)V_view.y + "  z " + (long)V_view.z);

		//jat_view.set_view_position(0.,50.-alpha*10,50.);
		//jat_view.set_view_direction(origin);
		//panel.label.setText(i + "  Time " + (long)orb1.t[i] + "  x " + (long)orb1.x[i] + "  y " + (long)orb1.y[i]);

		//apollo.set_attitude(alpha,alpha,0);
		//apollo.set_position(apolloorbit.x[i], apolloorbit.y[i], apolloorbit.z[i]);
		//apollo.set_attitude(0. * Math.PI, 0.6 * Math.PI, 0. * Math.PI);
		//jat_view.set_view_direction(apolloorbit.x[i], apolloorbit.y[i], apolloorbit.z[i]);
	}

	private CapturingCanvas3D createCanvas(String frames_path)
	{
		GraphicsConfigTemplate3D template = new GraphicsConfigTemplate3D();
		GraphicsConfiguration gc1 =
			GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(template);
		return new CapturingCanvas3D(gc1, frames_path);
	}

	public void init()
	{
	}

	public static void main(String[] args)
	{
		QuaternionTest sh = new QuaternionTest(); // Applet
		new MainFrame(sh, 800, 600);
	}

}
