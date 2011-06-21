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
import jat.eph.*;
import jat.util.*;
import jat.spacetime.*;
import jat.matvec.data.*;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import java.awt.*;
import javax.swing.*; // For Timer
import java.awt.event.*;
import com.sun.j3d.utils.applet.MainFrame;

/**
 * @author Tobias Berthold
 *  Date        :   11-18-2002
 *  Description :   Simulate Earth-Mars-Earth trajectory for Copernicus
 *
 */
public class Copernicus1 extends Applet implements ActionListener
{
	public CapturingCanvas3D c;
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e11); //100000000.0
	Planet3D mercury, venus, earth, mars, jupiter;
	Star3D sun;
	RGBAxis3D axis;
	ControlPanel panel;
	Copernicus_Trajectory ct1, ct2;
	//ColorCube3D spacecraft;
	ThreeDStudioObject spacecraft;
	Orbit orbit_mars, orbit_earth_moon;
	Ephemeris3D ephemeris_mercury, ephemeris_venus, ephemeris_earth, ephemeris_mars;
	Timer animControl;
	//int delayValue = 50; // milliseconds
	//Timer animControl = new Timer(delayValue, this);
	DE405 my_eph; // Ephemeris class
	double jd;
	Calendar cal;
	SimpleDateFormat sdf;
	Vector3d V = new Vector3d(0.f, 0.f, 0.0f);
	int i = 0;
	Matrix MRot;
	int steps = 300;

	public Copernicus1()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.Copernicus1", "Copernicus1");
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
		jd = cm.juliandate(2018, 6, 2, 12, 0, 0);
		MRot = new RotationMatrix(1, cm.Rad(Constants.eps));

		// Date related functions
		sdf = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss aaa");

		// 3D Objects
		BG_root = new BranchGroup();
		BG_root.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		BG_root.setBounds(bounds);
		float scale_factor = 1000.f; // exaggerate planets by factor
		TG_scene = new TransformGroup();
		TG_scene.addChild(axis = new RGBAxis3D(1.2f * cm.mars_a));
		TG_scene.addChild(sun = new Star3D(this, 1.0f));
		TG_scene.addChild(mercury = new Planet3D(this, Planet3D.MERCURY, scale_factor));
		TG_scene.addChild(venus = new Planet3D(this, Planet3D.VENUS, scale_factor));
		TG_scene.addChild(earth = new Planet3D(this, Planet3D.EARTH, scale_factor));
		TG_scene.addChild(mars = new Planet3D(this, Planet3D.MARS, scale_factor));
		TG_scene.addChild(spacecraft = new ThreeDStudioObject(this, "carrier.3DS", 50000.f));
		spacecraft.thrusters_on();
		//		TG_scene.addChild(spacecraft = new ColorCube3D(1.e7));
		BG_root.addChild(ephemeris_mercury = new Ephemeris3D(DE405_Body.MERCURY, jd, 100));
		BG_root.addChild(ephemeris_venus = new Ephemeris3D(DE405_Body.VENUS, jd, 300));
		BG_root.addChild(ephemeris_earth = new Ephemeris3D(DE405_Body.EARTH, jd, 365));
		BG_root.addChild(ephemeris_mars = new Ephemeris3D(DE405_Body.MARS, jd, 600));
		BG_root.addChild(TG_scene);
		// Copernicus Trajectory
		//ct1=new Copernicus_Trajectory(b+"Mission.txt",steps);
		ct1 = new Copernicus_Trajectory(b + "segment_1.txt", steps);
		BG_root.addChild(ct1);
		//ct2=new Copernicus_Trajectory(b+"segment_2.txt",steps);
		//BG_root.addChild(ct2);

		// Lights
		BG_root.addChild(jat_light.PointLight(bounds));
		//BG_root.addChild(jat_light.DirectionalLight(bounds));

		// View
		BG_vp = jat_view.view(0., 0., 1.e9, c);
		//BG_vp = jat_view.view(0., 0., 2.e8, c);
		//BG_vp = jat_view.view(2.e8,-3.e8, 1.e6, c);
		jat_view.set_FrontBackClipDistance(1.e5, 1.e10);

		// Behaviors
		jat_behavior.behavior(BG_root, BG_vp, bounds);
		jat_behavior.xyz_Behavior.set_translate(1.e7f);

		// Have Java 3D perform optimizations on this scene graph.
		BG_root.compile();

		VirtualUniverse universe = new VirtualUniverse();
		Locale locale = new Locale(universe);
		locale.addBranchGraph(BG_root);
		locale.addBranchGraph(BG_vp);

		// Use Timer for animation
		int delayValue = 20; // milliseconds
		animControl = new Timer(delayValue, this);
		animControl.start();
		//		System.out.println("i x y z xy_angle xz_angle");
	}

	// This method is called each time a timer event occurs
	public void actionPerformed(ActionEvent e)
	{
		//System.out.println("" + i);
		i++;

		//if (i > 300 && i < 580)
		//if(i > 10 && i<steps)
		if (i < steps)
		{

			jd = ct1.t[i];
			//jd+=ct1.t[1]-ct1.t[0];

			/*			
			// For Solarsystem view
			// set new viewing platform position
			Vector3d coord=new Vector3d(1.e9,0.5*Math.PI*i/steps,-0.25*Math.PI);
			Vector3d pos=CoordTransform.Spherical_to_Cartesian(coord);
			Vector3d look_to=new Vector3d(0.,0.,0.);
			jat_view.set_pos_direction(pos,look_to);
			*/

			// set new spacecraft position and attitude
			double x, y, z;
			x = ct1.ux[i];
			y = ct1.uy[i];
			z = ct1.uz[i];
			double xy_angle = -Math.atan2(z, Math.sqrt(x * x + y * y)); // Angle between vector and x-z-plane 
			double xz_angle = Math.atan2(y, x); // Angle between vector and x-z-plane 
			//	System.out.println(i+" "+x+" "+y+" "+z+" "+Constants.rad2deg*xy_angle+" "+Constants.rad2deg*xz_angle);
			spacecraft.set_pos_attitude(ct1.x[i], ct1.y[i], ct1.z[i], 0., xy_angle, xz_angle);
			
			// For above spacecraft view
			Vector3d look_to = new Vector3d(ct1.x[i], ct1.y[i], ct1.z[i]);
			Vector3d pos = new Vector3d(ct1.x[i], ct1.y[i], ct1.z[i] + 2.5e8);
			jat_view.set_pos_direction(pos, look_to);
			
		
			// Planet positions
			double mjd_tt = TimeUtils.JDtoMJD(jd);
			mercury.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.MERCURY, mjd_tt)));
			venus.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.VENUS, mjd_tt)));
			earth.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.EARTH, mjd_tt)));
			mars.set_position(MRot.times(my_eph.get_planet_pos(DE405_Body.MARS, mjd_tt)));

			// Take frame screenshot
			try
			{
				Thread.sleep(300);
				//System.out.println("waiting..");
			} catch (Exception f){};
			c.takePicture();			

		} //else
			//i = 0;

		// Update text in panel
		//V=jat_view.get_vp();
		cal = cm.Calendar(jd);
		String jds = pr.fwdts(jd, 10, 6);
		panel.label.setText(
			"JD: "
				+ jds
				+ "   Greg: "
				+ sdf.format(cal.getTime())
				+ "   Observer:  x:"
				+ (long)V.x
				+ "  y: "
				+ (long)V.y
				+ "  z: "
				+ (long)V.z);


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
		Copernicus1 c1 = new Copernicus1(); // Applet
		MainFrame m = new MainFrame(c1, 800, 600);
		m.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);
	}
}

/*
public class SimulationClock extends Behavior
{

	public void initialize()
	{
		//java.text.DateFormat DayTime;
		//DayTime = java.text.DateFormat.getDateTimeInstance(java.text.DateFormat.SHORT, java.text.DateFormat.SHORT);
	}

	public void processStimulus(Enumeration e)
	{
			//CalDate date = new CalDate( 2018, 3, 1, 0, 0, 0.0 );
			//double mjd = date.mjd();
			//EarthRef er=new EarthRef(mjd);
			//Matrix m=er.eci2ecef().inverse();
	}   //  end method

}   //  end class
*/
