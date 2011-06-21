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

package jat.demo.vr.PlanetGrid;

import jat.vr.*;
import jat.cm.*;
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
public class PlanetGrid extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Planet3D mars;
	PlanetGrid3D PG;
	Axis mainaxis;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	double alpha;

	public PlanetGrid()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		Orbit orb1 = null;
		String b = FileUtil.getClassFilePath("jat.demo.vr.PlanetGrid", "PlanetGrid");
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
		TG_scene.addChild(mars = new Planet3D(this, Planet3D.MARS));
		mars.set_attitude(Math.PI / 2., 0, 0);
		TG_scene.addChild(mainaxis = new Axis(cm.mars_radius-1000.0f));
		//TG_scene.addChild(PG=  new PlanetGrid3D(cm.mars_radius+5.));
		TG_scene.addChild(PG=  new PlanetGrid3D(cm.mars_radius+1.,5.,5.));
		BG_root.addChild(TG_scene);

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(0.f, 1.f, -1.f);
		//BG_root.addChild(jat_light.AmbientLight(bounds));
//		AmbientLight lightA = new AmbientLight(Colors.white);
//		lightA.setInfluencingBounds(bounds);
//		BG_root.addChild(lightA);

		Background background = new Background();
		background.setApplicationBounds(bounds);
		background.setColor(Colors.white);
		BG_root.addChild(background);


		// View
		BG_vp = jat_view.view(cm.mars_radius+1000., 0., cm.mars_radius+1000., c, 1.f, 100000.f);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(100.f);

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
		Vector3d V_view = new Vector3d( 0.f,0.f,0.0f);
		Transform3D T3D = new Transform3D();
		jat_view.TG_vp.getTransform(T3D);
		T3D.get(V_view);
		panel.label.setText("  x " + (long) V_view.x + "  y " + (long) V_view.y + "  z " + (long) V_view.z);

		alpha += 0.0001;
		//mars.set_attitude(Math.PI / 2., alpha, 0);
//		PG.set_attitude( 0., 0.,alpha);
		//jat_view.set_view_position(0.,50.-alpha*10,50.);
		//jat_view.set_view_direction(origin);
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
		PlanetGrid pg = new PlanetGrid(); // Applet
		new MainFrame(pg, 800, 600);
	}

}
