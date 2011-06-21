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
public class MarsOrbit extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	ThreeDStudioObject carrier;
	Planet3D earth, mars;
	RGBAxis3D mainaxis;
	Orbit carrierorbit = null;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	double alpha;
	int i = 0;
	int steps = 8500;

	public MarsOrbit()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.Copernicus1", "MarsOrbit");
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
		TG_scene.addChild(mars = new Planet3D(this, Planet3D.MARS));
		TG_scene.addChild(mainaxis = new RGBAxis3D(4000.0f));
		TG_scene.addChild(carrier = new ThreeDStudioObject(this, "carrier.3DS", 1.f));
		carrier.set_attitude(0. * Math.PI, 0.0 * Math.PI, 0.1 * Math.PI);
		//carrier.set_position(cm.earth_radius + 1000., -2000., 1000.);
		BG_root.addChild(TG_scene);
		BG_root.addChild(carrierorbit = new Orbit(cm.mars_radius + 500., 0.02, 5.0, 0.0, 0.0, 0.0, Colors.pink, steps));
		BG_root.addChild(carrierorbit = new Orbit(8000., 0.02, -45.0, 0.0, 0.0, 300.0, Colors.pink, steps));

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(0.f, 1.f, -0.5f);
		//BG_root.addChild(jat_light.AmbientLight(bounds));

		// View
		//BG_vp = jat_view.view(5418., -5239., 5304., c, 1.f, 100000.f);
		BG_vp = jat_view.view(5418., -5239., 5304., c, 1.f, 100000.f);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(100.f);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		// Use Timer for animation
		int delayValue = 50; // milliseconds
		animControl = new Timer(delayValue, this);
		animControl.start();
	}

	// This method is called each time a timer event occurs
	public void actionPerformed(ActionEvent e)
	{
		System.out.println("" + i);
		i++;
		if (i > steps)
			i = 0;

		if (i < 300)
		{
			alpha += 0.001;
			mars.set_attitude(Math.PI / 2., alpha, 0);

			Vector3d V_view = new Vector3d(0.f, 0.f, 0.0f);
			Transform3D T3D = new Transform3D();
			jat_view.TG_vp.getTransform(T3D);
			T3D.get(V_view);
			panel.label.setText("  x " + (long)V_view.x + "  y " + (long)V_view.y + "  z " + (long)V_view.z);

			//jat_view.set_view_position(0.,50.-alpha*10,50.);
			//jat_view.set_view_direction(origin);
			//panel.label.setText(i + "  Time " + (long)orb1.t[i] + "  x " + (long)orb1.x[i] + "  y " + (long)orb1.y[i]);
			//		panel.label.setText(i + "  Time " + (long)orb1.t[i] + "  x " + (long)orb1.x[i] + "  y " + (long)orb1.y[i]);

			//carrier.set_attitude(alpha,alpha,0);
			carrier.set_position(carrierorbit.x[i], carrierorbit.y[i], carrierorbit.z[i]);
			jat_view.set_view_direction(carrierorbit.x[i], carrierorbit.y[i], carrierorbit.z[i]);

			// Take frame screenshot
			try
			{
				Thread.sleep(10);
				//System.out.println("waiting..");
			} catch (Exception f)	{	};
			//c.takePicture();					
		}
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
		MarsOrbit sh = new MarsOrbit(); // Applet
		MainFrame m = new MainFrame(sh, 800, 600);
		m.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);
	}
}
