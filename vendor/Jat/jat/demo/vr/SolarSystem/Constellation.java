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

package jat.demo.vr.SolarSystem;

import jat.vr.*;
import jat.cm.*;
import jat.util.*;
import jat.eph.*;

import java.awt.*;
import java.util.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.media.j3d.Locale;
import javax.vecmath.*;
import java.text.SimpleDateFormat;
import com.sun.j3d.utils.applet.MainFrame;

/**
 * Solar system constellation simulation
 * Lets user use the keyboard to change time and view different planet constellations
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public class Constellation extends Applet //implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Star3D sun;
	Planet3D mercury, venus, earth, mars, moon;
	Axis axis;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e11);
	ControlPanel panel;
	public CapturingCanvas3D c;
	//Timer animControl;
	int timer_counter = 0, counter = 0;
	double alpha;
	Date dd;
	double jd;
	Calendar cal;
	SimpleDateFormat sdf;
	float scale_factor = 1000.f; // exaggerate planets by factor
	DE405 my_eph; // Ephemeris class

	public Constellation()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.SolarSystem", "Constellation");
		System.out.println(b);

		//Applet window
		setLayout(new BorderLayout());
		c = createCanvas(b + "frames/");
		add("Center", c);
		panel = new ControlPanel(BG_root);
		add("South", panel);

		// Ephemeris data
		String fs = FileUtil.file_separator();
		my_eph = new DE405(FileUtil.getClassFilePath("jat.eph","eph")+fs+"DE405"+fs);
		jd = cm.juliandate(2002, 2, 17, 12, 0, 0);

		// Date related functions
		sdf = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss");
		cal = new GregorianCalendar();
		dd = new Date(java.lang.System.currentTimeMillis());

		// 3D Objects
		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);
		TG_scene = new TransformGroup();
		TG_scene.addChild(axis = new Axis(1.2f * cm.mars_a));
		TG_scene.addChild(sun = new Star3D(this, 1.f));
		TG_scene.addChild(mercury = new Planet3D(this, Planet3D.MERCURY, scale_factor));
		TG_scene.addChild(venus = new Planet3D(this, Planet3D.VENUS, scale_factor));
		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH, scale_factor));
		TG_scene.addChild(mars = new Planet3D(this, Planet3D.MARS, scale_factor));
		//TG_scene.addChild(moon = new Planet3D(this, Planet3D.MOON));
		//moon.set_position(20000, 0, 0);
		BG_root.addChild(TG_scene);

		// Lights
		BG_root.addChild(jat_light.PointLight(bounds));

		// View
		BG_vp = jat_view.view(0., 0., 1.e9, c);
		jat_view.view.setFrontClipDistance(1.e5);
		jat_view.view.setBackClipDistance(1.e10);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(1.e7f);

		time_control_new tc = new time_control_new(this);
		tc.setSchedulingBounds(bounds);
		BG_vp.addChild(tc);

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Have Java 3D perform optimizations on this scene graph.
		//BG_root.compile();

		/*
				// Use Swing Timer
				int delayValue = 50; // millisecondes = 0.1 sec
				animControl = new Timer(delayValue, this);
				animControl.start();
		*/
	}

	/*
	 	public void actionPerformed(ActionEvent e)
		{
			// Perform the task
	
	        // Simulation time
	        dd.setTime(java.lang.System.currentTimeMillis());
	        cal.setTime(dd);
	        jd=cm.juliandate(cal);
	
	        // Update text in panel
	        //V=jat_view.get_vp();
	        String jds=pr.fwdts(jd,10,6);
	        //panel.label.setText("JD "+jds+ " Greg "+sdf.format( c.getTime())+"  x "+(long)V.x+ "  y "+(long)V.y+ "  z "+(long)V.z);
	        panel.label.setText("JD "+jds+ " Greg "+sdf.format( cal.getTime()));
	
			// Update Transform3D object
			alpha += 0.01;
			//earth.set_attitude(Math.PI / 2., alpha, 0);
			//moon.set_position(alpha * 1000, 10000, 0);
			//moon.set_attitude(alpha, 0, 0);
			mercury.set_position(my_eph.get_planet_posvel(eph.MERCURY, jd));
			venus.set_position(my_eph.get_planet_posvel(eph.VENUS, jd));
			earth.set_position(my_eph.get_planet_posvel(eph.EARTH, jd));
			mars.set_position(my_eph.get_planet_posvel(eph.MARS, jd));
		}
	*/

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
		// Create an instance
		Constellation tbvr = new Constellation();
		new MainFrame(tbvr, 800, 600);
	}

}
