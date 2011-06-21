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

package jat.matvec.io.gui.plotTools;

//import javax.swing.*;
import java.awt.*;

public class NotedPoint2D {

  private int[] Sc = new int[2];;
  private double[] Pl = new double[2];
  private double[] Pl_0 = new double[2];

  private int[] Sc_X0 = new int[2];;
  private double[] Pl_X0 = new double[2];
  private int[] Sc_Y0 = new int[2];
  private double[] Pl_Y0 = new double[2];

  private boolean isVisible;

  private Axe2D axe;
  private Grid2D grid;

  public NotedPoint2D() {
    setVisible(false);
  }

  public NotedPoint2D(double x, double y,Axe2D ax,Grid2D g) {

    axe = ax;
    grid = g;

    Pl[0] = x;
    Pl[1] = y;

    Pl_0 = axe.getPl0();

    Pl_X0[0] = Pl_0[0];
    Pl_X0[1] = y;

    Pl_Y0[0] = x;
    Pl_Y0[1] = Pl_0[1];

    Pl2ScConvert();

    setVisible(true);
  }

  private void Pl2ScConvert() {
    Sc = axe.Pl2Sc(Pl);
    Sc_X0 = axe.Pl2Sc(Pl_X0);
    Sc_Y0 = axe.Pl2Sc(Pl_Y0);
  }

  public void setVisible(boolean b) {
    isVisible = b;
  }

  public void draw(Graphics2D comp2D) {
    if (isVisible) {
      //Pl2ScConvert();
      comp2D.setColor(Color.black);
      comp2D.setFont(new Font("Arial",Font.PLAIN,12));
      comp2D.drawLine(Sc[0],Sc[1],Sc_X0[0],Sc_X0[1]);
      comp2D.drawLine(Sc[0],Sc[1],Sc_Y0[0],Sc_Y0[1]);

      comp2D.drawString(new String("(" + grid.troncatedStringX(Pl[0],2) + "," + grid.troncatedStringY(Pl[1],2) +")"),Sc[0],Sc[1]);
      //comp2D.drawString(new String(""+Pl[0]),Sc_Y0[0],Sc_Y0[1]);
      //comp2D.drawString(new String(""+Pl[1]),Sc_X0[0],Sc_X0[1]);

    }
  }

}

