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

package jat.demo.vr.attitude;

import jat.vr.*;
import jat.util.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; 
import java.awt.event.*;

public class Flexible extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ColorCube3D spacecraft1, spacecraft2;
	Axis axis;
	ControlPanel panel;
	public CapturingCanvas3D c;
	double x_inc;
	Timer animControl;

	public Flexible()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.attitude", "attitude");
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
		TG_scene.addChild(spacecraft1=new ColorCube3D(10000.f));
		TG_scene.addChild(spacecraft2=new ColorCube3D(10000.f));
		TG_scene.addChild(axis = new Axis(20000.0f));
		BG_root.addChild(TG_scene);

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(0.f, 1.f, -1.f);

		// View
		BG_vp = jat_view.view(-180000, 0., -1.e4, c);
		jat_view.set_view_position(-180000, 20000., -20000);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);

		// Use Swing Timer
		int delayValue = 100; // millisecondes = 0.1 sec
		animControl = new Timer(delayValue, this);
		animControl.start();

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();
	}

	public void actionPerformed(ActionEvent e)
	{
		//Transform3D satTrans = new Transform3D();
		// Perform the task
		// Update Transform3D object

		x_inc += 100.;
		spacecraft2.set_position(0., x_inc, 0.);

		//Stop animation at end
		if (x_inc >= 10000)
			animControl.stop();
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
		Flexible at = new Flexible(); // Applet
		new MainFrame(at, 800, 600);
	}
}
