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

package jat.demo.vr.EarthMoon;

import jat.vr.*;
import jat.cm.*;
import jat.util.*;

import java.awt.*;
import java.applet.Applet;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.image.*;
import com.sun.j3d.utils.geometry.*;

import javax.swing.*; 
import java.awt.event.*;

/**
 * @author Tobias Berthold
 *
 */
public class earth_moon_fly extends Applet implements ActionListener
{
	BranchGroup BG_root;
	BranchGroup BG_vp;
	TransformGroup TG_scene;
	Point3d origin = new Point3d(0.0f, 0.0f, 0.0f);
	BoundingSphere bounds = new BoundingSphere(origin, 1.e10); //100000000.0
	Planet3D earth, moon;
	ColorCube3D spacecraft;
	ThreeDStudioObject apollo;
	RGBAxis3D axis;
	ControlPanel panel;
	EricTrajectory tra;
	Orbit orb1;
	public CapturingCanvas3D c;
	Timer animControl;
	int i;
	int trapos=1390;

	public earth_moon_fly()
	{
		// Get path of this class, frames will be saved in subdirectory frames
		String b = FileUtil.getClassFilePath("jat.demo.vr.EarthMoon", "earth_moon_fly");
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
		TG_scene.addChild(moon = new Planet3D(this, Planet3D.MOON));
		//TG_scene.addChild(spacecraft = new ColorCube3D(1000.f));
		TG_scene.addChild(apollo = new ThreeDStudioObject(this, "Apollo/Apollo13.3ds", 10.f));
		//TG_scene.addChild(axis = new RGBAxis3D(15000.0f));
		BG_root.addChild(TG_scene);
		BG_root.addChild(tra = new EricTrajectory(b + "pab.txt"));
		//BG_root.addChild(orb1 = new Orbit(42162.0, 0.0, 0.0, 0.0, 0.0, 0.0));
		//float x_offset=800.f, y_offset=-880.f, z_offset=230.f;
		float x_offset=870.f, y_offset=620.f, z_offset=230.f;
		apollo.set_attitude(0., 0., Math.PI/2.);
		//apollo.set_position(x_offset, y_offset, z_offset);
		apollo.set_position(tra.coords[trapos * 3 + 0]+x_offset, tra.coords[trapos * 3 + 1]+y_offset, tra.coords[trapos * 3 + 2]+z_offset);
		
		// Star background
		Appearance app = new Appearance();
		// texturiertes bild mit Java 3D texture loader erstellen
		Texture tex = new TextureLoader(b+"stars3.jpg", this).getTexture();
		app.setTexture(tex);
		//Kugel wird erstellt um innen textur draufzuschmeissen(normalen nach innen!)
		Sphere sphere = new Sphere(5000000f, Primitive.GENERATE_TEXTURE_COORDS | Primitive.GENERATE_NORMALS_INWARD, app);
		BG_root.addChild(sphere);

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
		long SimTime=0;
		int DeltaT=40;
		//if(i>1350 && i<1450)
		if (i > 10 && i < 30)
		{
			System.out.println("" + i);
			//c.writeJPEG_ = true;
			//a.c.repaint();
		}
		i++;

		// general animation
		SimTime += DeltaT;

		// set new spacecraft position
		if (i < 1999)
		{
			SimTime = (long) tra.time[i];
			apollo.set_position(tra.coords[i * 3 + 0], tra.coords[i * 3 + 1], tra.coords[i * 3 + 2]);
		}

		// Update text in panel
		//		jat_view.TG_vp.getTransform(T_3D);
		//		T_3D.get(V_view);
		//panel.label.setText("Time " + SimTime + "  x " + (long) V_view.x + "  y " + (long) V_view.y + "  z " + (long) V_view.z);
		panel.label.setText("Time " + SimTime );
//
//		// set new viewing platform position
//		V_view.x += 100;
//		jat_view.set_view_position(V_view.x, V_view.y, V_view.z);

		// set a new viewing direction
		//jat_view.set_view_direction(V_view,origin);

		// Earth rotation
		earth.set_attitude(Math.PI / 2., 0.,i * 0.001f);

		// Moon position
		float inc_fac_cos = (float) Math.cos(cm.Rad(cm.moon_obl));
		float inc_fac_sin = (float) Math.sin(cm.Rad(cm.moon_obl));
		double x = cm.moon_a * (float) Math.cos(2 * Math.PI * SimTime / cm.moon_period);
		double y = cm.moon_a * (float) inc_fac_cos * Math.sin(2 * Math.PI * SimTime / cm.moon_period);
		double z = cm.moon_a * (float) inc_fac_sin * Math.sin(2 * Math.PI * SimTime / cm.moon_period);
		moon.set_position(x-2000., y-9500., z+1500.);
		//jat_view.set_view_position(0.,50.-i*10,50.);
		jat_view.set_view_position(-180000+i*100, 0., -1.e4);
		//jat_view.set_view_direction(origin);
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
		earth_moon_fly em = new earth_moon_fly(); // Applet
		MainFrame m=  new MainFrame(em, 800, 600);
		m.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);
	}

}
