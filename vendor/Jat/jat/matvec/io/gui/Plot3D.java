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

//import jat.matvec.data.Matrix;

import jat.matvec.io.gui.plotTools.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class Plot3D extends JPanel implements MouseListener, MouseMotionListener {

  protected Dimension defaultSize = new Dimension(400,400);

  public static int PIXEL = DataPlot3D.PIXEL;
  public static int DOT = DataPlot3D.DOT;
  public static int LINE = DataPlot3D.LINE;
  public static int DOTLINE = DataPlot3D.DOTLINE;
  public static int BAR = DataPlot3D.BAR;
  public static int DOTBAR = DataPlot3D.DOTBAR;
  public static int HIST = DataPlot3D.HIST;
  public static int GRID = DataPlot3D.GRID;

  protected DataPlot3D[] plots;
  protected Axe3D axe;
  protected Grid3D grid;
  protected PlotAttributes PA;
  protected NotedPoint3D np;

  protected double[][] X ;
  protected double[][] Y ;
  protected double[][] Z ;
  protected double X0 = 0;
  protected double Y0 = 0;
  protected double Z0 = 0;
  protected double[][] widthX;
  protected double[][] widthY;
  protected double[][] widthZ;

  protected int prevx;
  protected int prevy;

  protected void setAppearence() {
    setPreferredSize(defaultSize);
    setSize(defaultSize);
  }

  public void setPlotType(int i,int type) {
    PA.typeList[i] = type;
  }

  public void setPlotColor(int i,Color color) {
    PA.colorList[i] = color;
  }

  public void setPlotLegend(String[] legend) {
    PA.legend = legend;
  }

  public void update() {
    axe = new Axe3D(X0,Y0,Z0,X[0],Y[0],Z[0],this,new Eye3D(),PA.legend);
    grid = new Grid3D(axe);
    plots = new DataPlot3D[X.length];
    np = new NotedPoint3D();
    for (int i=0;i<X.length;i++) {
      plots[i] = new DataPlot3D(X[i],Y[i],Z[i],widthX[i],widthY[i],widthZ[i],axe,PA.typeList[i],PA.colorList[i]);
    }

    addMouseListener(this);
    addMouseMotionListener(this);

    repaint();
  }

  public void paint(Graphics comp) {
    Graphics2D comp2D = (Graphics2D)comp;
    comp2D.setColor(getBackground());
    comp2D.fillRect(0,0,getSize().width,getSize().height);
    grid.draw(comp2D);
    axe.draw(comp2D);
    for (int i=0;i<plots.length;i++) {
      plots[i].draw(comp2D);
    }
    np.draw(comp2D);

    setBackground(Color.white);
  }

  public void mouseDragged(MouseEvent e) {

    int x = e.getX();
    int y = e.getY();

    double newTheta = axe.getEye3D().getTheta() - (double)(Math.PI*0.01*(x-prevx));
    double newPhi = axe.getEye3D().getPhi() - (double)(Math.PI*0.01*(prevy-y));

    axe.setEye3D(new Eye3D(newTheta,newPhi));

    for (int i=0;i<plots.length;i++) {
      plots[i].setAxe3D(axe);
    }

    repaint();

    prevx = x;
    prevy = y;
    e.consume();
  }

  public void mousePressed(MouseEvent e) {
    prevx = e.getX();
    prevy = e.getY();
    e.consume();
  }

  public void mouseClicked(MouseEvent e) {
  }
  public void mouseReleased(MouseEvent e) {
  }
  public void mouseEntered(MouseEvent e) {
  }
  public void mouseExited(MouseEvent e) {
  }
  public void mouseMoved(MouseEvent e) {

    np.setVisible(false);

    int Sc_x = e.getX();
    int Sc_y = e.getY();

    int[] Sc_XY;
    int Sc_X;
    int Sc_Y;

    double[] Pl_XYZ;
    double Pl_X;
    double Pl_Y;
    double Pl_Z;

    for (int i=0;i<plots.length;i++) {
      Coordinates3D[] all = plots[i].getCoord();
      for (int j = 0;j<all.length;j++) {
        Sc_XY = all[j].getSc();
        Sc_X = Sc_XY[0];
        Sc_Y = Sc_XY[1];
        if ((Math.abs(Sc_x-Sc_X)<5)&&(Math.abs(Sc_y-Sc_Y)<5)) {
          Pl_XYZ = all[j].getPl();
          Pl_X = Pl_XYZ[0];
          Pl_Y = Pl_XYZ[1];
          Pl_Z = Pl_XYZ[2];
          np = new NotedPoint3D(Pl_X,Pl_Y,Pl_Z,axe,grid);
        }
      }
    }
    repaint();
    e.consume();
  }
}
