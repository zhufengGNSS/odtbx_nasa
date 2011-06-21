/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2003 United States Government as represented by the
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
 */

package jat.matvec.io.gui;


import javax.swing.JFrame;
import javax.swing.JPanel;

public class FrameView extends JFrame {

  public FrameView(JPanel panel) {
    setContentPane(panel);
    pack();
    setVisible(true);
  }

  public FrameView(String title,JPanel panel) {
    super(title);
    setContentPane(panel);
    pack();
    setVisible(true);
  }

   public FrameView(JPanel[] panels) {
    JPanel panel = new JPanel();
    for (int i=0;i<panels.length;i++) {
      panel.add(panels[i]);
    }
    setContentPane(panel);
    pack();
    setVisible(true);
  }


}

