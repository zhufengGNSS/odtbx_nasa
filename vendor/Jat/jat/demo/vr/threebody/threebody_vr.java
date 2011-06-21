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

package jat.demo.vr.threebody;

import jat.vr.*;
import jat.cm.*;
import jat.util.*;
import jat.alg.integrators.*;
import jat.matvec.data.*;
import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*; 
import java.awt.event.*;

/**
 * Three body gravity simulation
 * 
 * @author Tobias Berthold
 * @version 1.0
 */
public class threebody_vr extends Applet implements ActionListener, Printable
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Sphere3D mass1, mass2, mass3;
	Planet3D earth, moon;
	Star3D sun;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	static double[][] store_states;
	int timer_counter = 0, counter = 0;
	static double m1, m2, m3; // masses
	static int number_of_points; // how many points
	static ThreeBody tb; // Threebody object, holds derivatives
	double alpha;

	public threebody_vr()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.threebody", "threebody_vr");
		System.out.println(b);

		setLayout(new BorderLayout());
		//Canvas3D c = createCanvas();
		c = createCanvas(b + "frames/");
		add("Center", c);

		panel = new ControlPanel(BG_root);
		add("South", panel);

		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);

		// 3D Objects
		TG_scene = new TransformGroup();
		TG_scene.addChild(mass1 = new Sphere3D( (float) m1 / 10, Colors.red, 5.f, 0.f, 0.f));
		TG_scene.addChild(mass2 = new Sphere3D( (float) m2 / 10, Colors.green, 0.f, 5.f, 0.f));
		TG_scene.addChild(mass3 = new Sphere3D( (float) m3 / 10, Colors.blue, 0.f, 0.f, 5.f));

		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH,0.2f));
		TG_scene.addChild(moon  = new Planet3D(this, Planet3D.MOON));
		TG_scene.addChild(sun  = new Star3D(this, 0.1f));
		moon.set_position(20000,0,0);
		sun.set_position(0,0,20000);

		BG_root.addChild(TG_scene);
		BG_root.addChild(new Axis(10000.0f));

		// Lights
		BG_root.addChild(jat_light.DirectionalLight(bounds));
		//light.setDirection(0.f, -100.f, -100.f);
		//BG_root.addChild(jat_light.AmbientLight(bounds));

		// View
		BG_vp = jat_view.view(0, 0., 200., c, 1.f, 100000.f);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(1.f);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		// Use Swing Timer
		int delayValue = 50; // millisecondes = 0.1 sec
		animControl = new Timer(delayValue, this);
		animControl.start();

	}

	public void actionPerformed(ActionEvent e)
	{
		// Perform the task
		// Update Transform3D object
		alpha += 0.01;
		earth.set_attitude( Math.PI / 2.,alpha, 0);
		moon.set_position(alpha*1000,10000,0);
		moon.set_attitude( alpha, 0,0);


		//System.out.println("" + timer_counter);
		if (timer_counter < number_of_points)
		{
			VectorN V = new VectorN(store_states[timer_counter]);
			VectorN body1 = V.get(0, 6);
			VectorN body2 = V.get(6, 6);
			VectorN body3 = V.get(12, 6);
			//V.print();	body1.print();	body2.print();	body3.print();
			mass1.set_position(body1);
			mass2.set_position(body2);
			mass3.set_position(body3);
			VectorN cma = tb.center_of_mass(store_states[timer_counter]);
			//jat_view.set_view_position(cma.get(0), cma.get(1), cma.get(2) + 200);
			panel.label.setText("Energy: " + tb.Energy(store_states[timer_counter]));
		}
		timer_counter++;
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

	/** Implements the Printable interface to get the data out of the propagator
	 *  This method is executed by the propagator at each integration step.
	 * @param t Time.
	 * @param y Data array.
	 */
	public void print(double t, double[] y)
	{

		// print to the screen for warm fuzzy feeling
		//System.out.println(t + " " + y[0] + " " + y[1] + " " + y[2] + " " + counter);

		// Store states in time
		for (int i = 0; i < 18; i++)
			store_states[counter][i] = y[i];
		counter++;
	}

	public static void main(String[] args)
	{
		// Masses
		m1 = 5.;
		m2 = 10.;
		m3 = 30.;

		number_of_points = 300;
		store_states = new double[number_of_points + 2][18];

		// set the initial and final time
		double t0 = 0.0;
		double tf = 500.0;

		double step_size = (tf - t0) / number_of_points;

		// create an RungeKutta8 integrator with step-size
		RungeKutta8 rk8 = new RungeKutta8(step_size);

		// create an instance
		threebody_vr tbvr = new threebody_vr();

		// Threebody dynamics
		tb = new ThreeBody(1., m1, m2, m3);

		// Initial values
		double[] x0 = new double[18];
		// Mass 1
		x0[0] = 10.0;
		x0[1] = 0.0;
		x0[2] = 0.0;
		x0[3] = 0.0;
		x0[4] = 1.0;
		x0[5] = 0.0;
		// Mass 2
		x0[6] = -5.0;
		x0[7] = 8.66025;
		x0[8] = 0.0;
		x0[9] = -0.866025;
		x0[10] = -0.5;
		x0[11] = 0.0;
		// Mass 3
		x0[12] = -5.0;
		x0[13] = -8.66025;
		x0[14] = 0.0;
		x0[15] = 0.866025;
		x0[16] = -0.5;
		x0[17] = 0.0;

		// integrate the equations
		rk8.integrate(t0, x0, tf, tb, tbvr, true);

		new MainFrame(tbvr, 800, 600);
	}

}
