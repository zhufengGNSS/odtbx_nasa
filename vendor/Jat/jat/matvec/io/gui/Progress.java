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
import javax.swing.JProgressBar;

public class Progress extends JFrame {

  private JProgressBar progress;
  private int min;
  private int max;
  private int val;
  private JPanel pane;

  public Progress(int m,int M) {

    min = m;
    max = M;

    val = min;

    this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    pane = new JPanel();

    progress = new JProgressBar(min,max);
    progress.setValue(val);
    progress.setString(null);

    pane.add(progress);
    this.setContentPane(pane);
    this.pack();
    this.setVisible(true);
  }

  public void setValue(int n) {
    val = n;
    progress.setValue(val);
    if (val>=max) setVisible(false);
  }

}