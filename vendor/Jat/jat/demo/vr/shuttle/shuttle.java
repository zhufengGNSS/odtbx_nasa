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

package jat.demo.vr.shuttle;

import jat.vr.*;
import jat.cm.*;
import jat.util.*;

import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; 
import java.awt.event.*;

/**
 * @author Tobias Berthold
 *
 */
public class shuttle extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	Planet3D earth;
	RGBAxis3D axis;
	LightWaveObject shuttle;
	Orbit orb1;
	Axis mainaxis, shuttleaxis;
	public CapturingCanvas3D c;
	Timer animControl;
	int i;
	int steps = 2000; // steps in orbit

	public shuttle()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.shuttle", "shuttle");
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
		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH));
		TG_scene.addChild(shuttle = new LightWaveObject(this, "SpaceShuttle.lws", 200.f));
		//TG_scene.addChild(axis = new RGBAxis3D(15000.0f));
		TG_scene.addChild(shuttleaxis = new Axis(10.0f));
		BG_root.addChild(TG_scene);
		BG_root.addChild(orb1 = new Orbit(cm.earth_radius + 1000., 0., 22.0, 0.0, 0.0, 0.0, Colors.pink, steps));
		i=1000;
		shuttle.set_position(orb1.x[i], orb1.y[i], orb1.z[i]);
		shuttle.set_scale(200.);
		double alpha = 1.1;
		shuttle.set_attitude(alpha, alpha, 0.);
		//shuttleaxis.set_position(cm.earth_radius + 500., 0., 0.);

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(100.f, 1.f, -1.f);
		Background background = new Background();
		background.setApplicationBounds(bounds);
		background.setColor(Colors.gray);
		BG_root.addChild(background);

		// View
		BG_vp = jat_view.view(-180000, 0., -1.e4, c);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(500);


		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		// Use Swing Timer for animation
		int delayValue = 50; // milliseconds
		animControl = new Timer(delayValue, this);
		animControl.start();

	}

	// This method is called each time a timer event occurs
	public void actionPerformed(ActionEvent e)
	{
		Transform3D T_3D = new Transform3D();
		if (i > 10 && i < 30)
		{
			System.out.println("" + i);
			//c.writeJPEG_ = true;
			//a.c.repaint();
		}

		i++;
		if (i > steps)
		{
			i = 0;
		}
		panel.label.setText(i + "  Time " + (long)orb1.t[i] + "  x " + (long)orb1.x[i] + "  y " + (long)orb1.y[i]);
		shuttle.set_position(orb1.x[i], orb1.y[i], orb1.z[i]);

		// Earth rotation
		earth.set_attitude(Math.PI / 2., 0.,i * 0.001f);
	}

	//private Canvas3D createCanvas()
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
		shuttle em = new shuttle(); // Applet
		new MainFrame(em, 800, 600);
	}

}
