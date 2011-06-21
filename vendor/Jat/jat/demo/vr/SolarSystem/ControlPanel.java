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

import java.awt.*;
import javax.media.j3d.*;

public class ControlPanel extends java.awt.Panel
{
    //public Label label;

    //private TextArea instructions;


    public ControlPanel(Group group)
    {
        super();

        setName("ControlPanel");
        setLayout(null);
        setBackground(java.awt.Color.lightGray);
        setSize(568, 150);
        /*
        label = new java.awt.Label();
        label.setName("sites");
        label.setText("Sites:");
        //label.setBounds(8, 2, 125, 30);
        add(label, label.getName());
        label.setAlignment(Label.LEFT);
        */

		//{{INIT_CONTROLS
		setLayout(new BorderLayout(0,0));
		setBackground(java.awt.Color.lightGray);
		setSize(0,0);
		label.setText("text");
		add(BorderLayout.CENTER, label);
		//}}
	}

   public void initialize()
   {
   }
	//{{DECLARE_CONTROLS
	java.awt.Label label = new java.awt.Label();
	//}}
}

