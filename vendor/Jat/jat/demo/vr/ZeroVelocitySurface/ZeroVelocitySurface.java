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

package jat.demo.vr.ZeroVelocitySurface;

import jat.vr.*;
import jat.util.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; 
import java.awt.event.*;

/**
 * Zero Velocity surfaces
 * 
 * @author Tobias Berthold
 * @version 1.0
 */
public class ZeroVelocitySurface extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Planet3D earth, moon;
	Star3D sun;
	Sphere3D s;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	double alpha;

	private Surface3D SF;

	public ZeroVelocitySurface()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.ZeroVelocitySurface", "ZeroVelocitySurface");
		System.out.println(b);

		setLayout(new BorderLayout());
		c = createCanvas(b + "frames/");
		add("Center", c);

		panel = new ControlPanel(BG_root);
		add("South", panel);

		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);

		// 3D Objects
		TG_scene = new TransformGroup();
		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH,0.5f));
		TG_scene.addChild(moon  = new Planet3D(this, Planet3D.MOON));
		TG_scene.addChild(sun  = new Star3D(this, 0.1f));
		TG_scene.addChild(s = new Sphere3D( 5000.f, Colors.turquoise, 5.f, 0.f, 0.f));
		moon.set_position(20000,0,0);
		sun.set_position(0,0,20000);
		s.set_position(0,20000,0);

		BG_root.addChild(TG_scene);
		BG_root.addChild(new Axis(10000.0f));
        BG_root.addChild(SF=new Surface3D());

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		//light.setDirection(0.f, -100.f, -100.f);
		//BG_root.addChild(jat_light.AmbientLight(bounds));

		// View
		BG_vp = jat_view.view(0, 0., 20000., c, 1.f, 1000000.f);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(10000.f);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		// Use Swing Timer
		int delayValue = 50; // milliseconds
		animControl = new Timer(delayValue, this);
		animControl.start();

	}

	public void actionPerformed(ActionEvent e)
	{
		// Perform the task
		// Update Transform3D object
		alpha += 0.01;
		earth.set_attitude( Math.PI / 2.,0.,alpha);
		moon.set_position(alpha*1000,10000,0);
		moon.set_attitude( alpha, 0,0);

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

		// create an instance
		ZeroVelocitySurface ZVS = new ZeroVelocitySurface();
		MainFrame m=new MainFrame(ZVS, 800, 600);
		m.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);
	}

}
