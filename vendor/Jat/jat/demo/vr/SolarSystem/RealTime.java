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
import jat.matvec.data.*;
import jat.util.*;
import jat.eph.*;
import jat.spacetime.*;
import java.awt.*;
import java.util.Date;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import java.text.SimpleDateFormat;
import com.sun.j3d.utils.applet.MainFrame;
import javax.swing.*;
import java.awt.event.*;

/**
 * Realtime Solar system simulation
 * Uses new timer
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public class RealTime extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Star3D sun;
	Planet3D mercury, venus, earth, mars, moon;
	Ephemeris3D ephemeris_mercury, ephemeris_venus, ephemeris_earth, ephemeris_mars;
	RGBAxis3D axis;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e11);
	ControlPanel panel;
	public CapturingCanvas3D c;
	Timer animControl;
	int timer_counter = 0, counter = 0;
	double alpha;
	Date dd;
	double jd;
	Calendar cal;
	SimpleDateFormat sdf;
	float scale_factor = 1000.f; // exaggerate planets by factor
	DE405 my_eph; // Ephemeris class
	Matrix MRot;

	public RealTime()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.SolarSystem", "RealTime");
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
		MRot = new RotationMatrix(1, cm.Rad(Constants.eps));
//		jd = cm.juliandate(2002, 2, 17, 12, 0, 0);

		// Date related functions
		sdf = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss aaa");
		cal = new GregorianCalendar();
		dd = new Date(java.lang.System.currentTimeMillis());
		dd.setTime(java.lang.System.currentTimeMillis());
		cal.setTime(dd);
		jd = cm.juliandate(cal);

		// 3D Objects
		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);
		TG_scene = new TransformGroup();
		TG_scene.addChild(axis = new RGBAxis3D(1.2f * cm.mars_a));
		TG_scene.addChild(sun = new Star3D(this, 1.f));
		TG_scene.addChild(mercury = new Planet3D(this, Planet3D.MERCURY, scale_factor));
		TG_scene.addChild(venus = new Planet3D(this, Planet3D.VENUS, scale_factor));
		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH, scale_factor));
		TG_scene.addChild(mars = new Planet3D(this, Planet3D.MARS, scale_factor));
		BG_root.addChild(ephemeris_mercury = new Ephemeris3D(DE405_Body.MERCURY, jd, 100));
		BG_root.addChild(ephemeris_venus = new Ephemeris3D(DE405_Body.VENUS, jd, 300));
		BG_root.addChild(ephemeris_earth = new Ephemeris3D(DE405_Body.EARTH, jd, 365));
		BG_root.addChild(ephemeris_mars = new Ephemeris3D(DE405_Body.MARS, jd, 600));
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

	// TODO: timer freezes when system time set to earlier time
	public void actionPerformed(ActionEvent e)
	{
		// Simulation time
		dd.setTime(java.lang.System.currentTimeMillis());
		cal.setTime(dd);
		jd = cm.juliandate(cal);

		// Update text in panel
		String jds = pr.fwdts(jd, 10, 6);
		panel.label.setText("JD: " + jds + "   Greg: " + sdf.format(cal.getTime()));
		//panel.label.setText("JD "+jds+ " Greg "+sdf.format( c.getTime())+"  x "+(long)V.x+ "  y "+(long)V.y+ "  z "+(long)V.z);

		// Update Transform3D object
		alpha += 0.01;
		//earth.set_attitude(Math.PI / 2., alpha, 0);
		//moon.set_position(alpha * 1000, 10000, 0);
		//moon.set_attitude(alpha, 0, 0);
		double mjd_tt = TimeUtils.JDtoMJD(jd);
		mercury.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.MERCURY, mjd_tt)));
		venus.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.VENUS, mjd_tt)));
		earth.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.EARTH, mjd_tt)));
		mars.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.MARS, mjd_tt)));
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
		// Create an instance
		RealTime rt = new RealTime();
		new MainFrame(rt, 800, 600);
	}
}

/*
import java.util.*;
import java.applet.Applet;
import java.util.Date;
import java.util.Enumeration; //  09-17-01
import java.text.SimpleDateFormat;
import javax.media.j3d.*;
import javax.vecmath.*;
import jat.vr.*;
import jat.eph.*;
import jat.util.*;
import jat.cm.*;
import jat.env.*;
import jat.matvec.data.*;

public class SimulationClock extends Behavior
{
	// Animation time specific
	WakeupCriterion yawn;
	long SimTime;
	int DeltaT;
	Date dd;
	double jd;
	Calendar c;
	SimpleDateFormat sdf;
	// Java3D objects
	solar_real app;
	private Transform3D T_3D = new Transform3D();
	Vector3d V = new Vector3d(0.f, 0.f, 0.0f);
	Point3d P_origin = new Point3d(0.0, 0.0, 0.0);
	Vector3d V_origin = new Vector3d(0.0, 0.0, 0.0);
	Vector3d away = new Vector3d(100000.0, 0.0, 0.0);
	VectorN V_mercury, V_venus, V_earth, V_mars, V_jupiter;
	double[] origin = { 0., 0., 0. };
	ControlPanel panel;
	int i;
	eph my_eph; // Ephemeris class

	public SimulationClock(int ts, solar_real app, ControlPanel panel)
	{
		DeltaT = ts;
		//this.bod=bod;
		this.app = app;
		this.panel = panel;
		yawn = new WakeupOnElapsedTime(DeltaT);
		dd = new Date(java.lang.System.currentTimeMillis());

		c = new GregorianCalendar();
		c.setTime(dd);
		System.out.println(c.get(Calendar.YEAR));
		System.out.println(c.get(Calendar.MONTH));
		System.out.println(c.get(Calendar.DAY_OF_MONTH));
		System.out.println(c.get(Calendar.HOUR));
		System.out.println(c.get(Calendar.MINUTE));
		System.out.println(c.get(Calendar.SECOND));

		sdf = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss");
		System.out.println(sdf.format(c.getTime()));

		// initial position of sun
		app.sun.set_position(origin);

		jd = cm.juliandate(2002, 2, 17, 12, 0, 0);

		my_eph = new eph(Path.JATPATH + "\\DE405\\");

	}

	public void initialize()
	{
		//java.text.DateFormat DayTime;
		//DayTime = java.text.DateFormat.getDateTimeInstance(java.text.DateFormat.SHORT, java.text.DateFormat.SHORT);
		wakeupOn(yawn);
		SimTime = 0;
	}


} //  end class
*/
