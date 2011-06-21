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
import jat.attitude.eom.*;
import jat.alg.integrators.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; 
import java.awt.event.*;

public class AttitudeTimer extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	ColorCube3D spacecraft;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	float quat_values[][];
	int curr_pt = 0;
	Timer animControl;
	int numberOfPts;

	public AttitudeTimer()
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
		TG_scene.addChild(spacecraft = new ColorCube3D(10000.f));
		BG_root.addChild(TG_scene);
		BG_root.addChild(new Axis(20000.0f));

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		jat_light.setDirection(100.f, 0.f, 0.f);

		// View
		BG_vp = jat_view.view(-180000, 0., -1.e4, c);

		// Behaviors
		//jat_behavior.behavior(this, BG_root, BG_vp, bounds);
		jat_behavior.behavior(BG_root, BG_vp, bounds);

		// Attitude Simulation Numerical Integration
		double Ixx = 10.42;
		double Iyy = 35.42;
		double Izz = 41.67;
		double timeDuration = 20;
		double time_step = 0.1;

		double[] x0 = new double[11];
		// Initial Conditions
		x0[0] = 0;
		x0[1] = 0;
		x0[2] = 1;
		x0[3] = 0;
		x0[4] = 0;
		x0[5] = 0;
		x0[6] = 1;
		x0[7] = 0;
		x0[8] = 0;
		x0[9] = 0;
		x0[10] = 0;

		double J = 10;
		double angle = 20;
		double psi = 30;
		double phi = 50;
		double theta = 60;
		numberOfPts = (int) (timeDuration / time_step) + 1;
		double tf = timeDuration;
		double t0 = 0;
		quat_values = new float[5][numberOfPts + 1]; // Don't forget to do the "Plus 1"!!
		// Start Simulation
		FourRWManeuver eomObject4 = new FourRWManeuver(time_step, Ixx, Iyy, Izz, J, angle, psi, phi, theta, quat_values);

		RungeKutta8 rk8 = new RungeKutta8(time_step);
		rk8.integrate(t0, x0, tf, eomObject4, true);
		quat_values = eomObject4.getQuaternion();

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
		Transform3D satTrans = new Transform3D();
		// Perform the task
		// Update Transform3D object

		satTrans.setRotation(
			new Quat4f(quat_values[1][curr_pt], quat_values[2][curr_pt], quat_values[3][curr_pt], quat_values[4][curr_pt]));
		spacecraft.set_attitude(satTrans);
		//Increment to next point
		curr_pt = curr_pt + 1;
		//Stop animation at end
		if (curr_pt >= numberOfPts)
			animControl.stop();
	}
	//private Canvas3D createCanvas()
	private CapturingCanvas3D createCanvas(String frames_path)
	{
		GraphicsConfigTemplate3D template = new GraphicsConfigTemplate3D();
		GraphicsConfiguration gc1 =
			GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(template);
		//return new Canvas3D(gc1);
		return new CapturingCanvas3D(gc1, frames_path);
	}

	public void init()
	{

	}

	public static void main(String[] args)
	{
		AttitudeTimer at = new AttitudeTimer(); // Applet
		new MainFrame(at, 800, 600);
	}

}
