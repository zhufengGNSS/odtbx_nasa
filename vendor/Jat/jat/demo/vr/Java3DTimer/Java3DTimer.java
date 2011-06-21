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

package jat.demo.vr.Java3DTimer;

import jat.vr.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;

/**
 * @author Tobias Berthold
 *
 * This example shows a VR animation using a Java3D timer for the animation.
 * It seems that Java timers based on threads are more stable, so it is 
 * recommended to use that method instead.
 * 
 */
public class Java3DTimer extends Applet
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ColorCube3D cube;
	Axis axis;
	ControlPanel panel;

	public Java3DTimer()
	{
		// Applet window
		setLayout(new BorderLayout());
		Canvas3D c = createCanvas();
		add("Center", c);
		panel = new ControlPanel(BG_root);
		add("South", panel);

		// 3D Objects
		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);
		TG_scene = new TransformGroup();
		TG_scene.addChild(cube = new ColorCube3D(1000.f));
		TG_scene.addChild(axis = new Axis(1500.0f));
		BG_root.addChild(TG_scene);

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(0.f, 1.f, -1.f);

		// View
		BG_vp = jat_view.view(-10000.f, 3000.f, -4000.f, c);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(500);

		// Animation
		SimulationClock SimClock = new SimulationClock(40, this, panel);
		SimClock.setSchedulingBounds(bounds);
		BG_vp.addChild(SimClock);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();
	}

	private Canvas3D createCanvas()
	{
		GraphicsConfigTemplate3D template = new GraphicsConfigTemplate3D();
		GraphicsConfiguration gc1 =
			GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(template);
		return new Canvas3D(gc1);
	}

	public void init()
	{
	}

	public static void main(String[] args)
	{
		Java3DTimer JT = new Java3DTimer(); // Applet
		new MainFrame(JT, 800, 600);
	}

}
