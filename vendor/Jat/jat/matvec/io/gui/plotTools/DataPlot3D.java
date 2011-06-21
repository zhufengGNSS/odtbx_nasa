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

public class DataPlot3D {

  public static int PIXEL = PlotAttributes.PIXEL;
  public static int DOT = PlotAttributes.DOT;
  public static int LINE = PlotAttributes.LINE;
  public static int DOTLINE = PlotAttributes.DOTLINE;
  public static int BAR = PlotAttributes.BAR;
  public static int DOTBAR = PlotAttributes.DOTBAR;
  public static int HIST = PlotAttributes.HIST;
  public static int GRID = PlotAttributes.GRID;

  private Coordinates3D[] points;
  private double[] widthX;
  private double[] widthY;
  private double[] widthZ;
  private int type;
  private Color color;
  private Axe3D axe;

  public DataPlot3D(double[] x,double[] y,double[] z,double[] wX,double[] wY,double[] wZ,Axe3D ax,int typ,Color col) {
    axe = ax;
    points = new Coordinates3D[x.length];
    for (int i=0;i<x.length;i++) {
      points[i] = new Coordinates3D(x[i],y[i],z[i],axe);
    }
    color = col;
    type = typ;
    widthX = wX;
    widthY = wY;
    widthZ = wZ;
  }

  public Coordinates3D[] getCoord() {
    return points;
  }

  public void setAxe3D(Axe3D ax) {
    axe = ax;
    for (int i=0;i<points.length;i++) {
      points[i].setAxe3D(axe);
    }
  }

  private void drawDots(Graphics2D comp2D) {
    int d = PlotAttributes.dotSize;
    for (int i=0;i<points.length;i++) {
      comp2D.setColor(color);
      comp2D.fillOval(points[i].getSc()[0]-(int)(d/2),points[i].getSc()[1]-(int)(d/2),d,d);
      comp2D.drawOval(points[i].getSc()[0]-(int)(d/2),points[i].getSc()[1]-(int)(d/2),d,d);
      comp2D.setColor(Color.white);
      comp2D.fillOval(points[i].getSc()[0]-(int)(d/4),points[i].getSc()[1]-(int)(d/4),(int)(d/4),(int)(d/4));
    }
  }

  public void draw(Graphics2D comp2D) {
    switch (type) {
      case 0:
        comp2D.setColor(color);
        for (int i=0;i<points.length-1;i++) {
          comp2D.drawLine(points[i].getSc()[0],points[i].getSc()[1],points[i].getSc()[0],points[i].getSc()[1]);
        }
      break;
      case 1:
        drawDots(comp2D);
      break;
      case 2:
        comp2D.setColor(color);
        for (int i=0;i<points.length-1;i++) {
          comp2D.drawLine(points[i].getSc()[0],points[i].getSc()[1],points[i+1].getSc()[0],points[i+1].getSc()[1]);
        }
      break;
      case 3:
        drawDots(comp2D);
        comp2D.setColor(color);
        for (int i=0;i<points.length-1;i++) {
          comp2D.drawLine(points[i].getSc()[0],points[i].getSc()[1],points[i+1].getSc()[0],points[i+1].getSc()[1]);
        }
      break;
      case 4:
        comp2D.setColor(color);
        for (int i=0;i<points.length;i++) {
          Coordinates3D pointbase = points[i].copy();
          double[] coord = pointbase.getPl();
          coord[2] = 0;
          pointbase.setPl(coord);

          comp2D.drawLine(points[i].getSc()[0],points[i].getSc()[1],pointbase.getSc()[0],pointbase.getSc()[1]);
        }
      break;
      case 5:
        drawDots(comp2D);
        comp2D.setColor(color);
        for (int i=0;i<points.length;i++) {
          Coordinates3D pointbase = points[i].copy();
          double[] coord = pointbase.getPl();
          coord[2] = 0;
          pointbase.setPl(coord);

          comp2D.drawLine(points[i].getSc()[0],points[i].getSc()[1],pointbase.getSc()[0],pointbase.getSc()[1]);
        }
      break;
      case 6:
        comp2D.setColor(color);
        for (int i=0;i<points.length;i++) {
          Coordinates3D pointbase = points[i].copy();
          double[] coord = pointbase.getPl();
          coord[2] = 0;
          pointbase.setPl(coord);

          Coordinates3D pt1 = pointbase.addVector(-widthX[i]/2,-widthY[i]/2,0);
          Coordinates3D pt2 = pointbase.addVector(widthX[i]/2,-widthY[i]/2,0);
          Coordinates3D pt3 = pointbase.addVector(-widthX[i]/2,widthY[i]/2,0);
          Coordinates3D pt4 = pointbase.addVector(widthX[i]/2,widthY[i]/2,0);
          Coordinates3D pt5 = points[i].addVector(-widthX[i]/2,-widthY[i]/2,0);
          Coordinates3D pt6 = points[i].addVector(widthX[i]/2,-widthY[i]/2,0);
          Coordinates3D pt7 = points[i].addVector(-widthX[i]/2,widthY[i]/2,0);
          Coordinates3D pt8 = points[i].addVector(widthX[i]/2,widthY[i]/2,0);

          comp2D.drawLine(pt1.getSc()[0],pt1.getSc()[1],pt2.getSc()[0],pt2.getSc()[1]);
          comp2D.drawLine(pt2.getSc()[0],pt2.getSc()[1],pt4.getSc()[0],pt4.getSc()[1]);
          comp2D.drawLine(pt4.getSc()[0],pt4.getSc()[1],pt3.getSc()[0],pt3.getSc()[1]);
          comp2D.drawLine(pt3.getSc()[0],pt3.getSc()[1],pt1.getSc()[0],pt1.getSc()[1]);

          comp2D.drawLine(pt5.getSc()[0],pt5.getSc()[1],pt6.getSc()[0],pt6.getSc()[1]);
          comp2D.drawLine(pt6.getSc()[0],pt6.getSc()[1],pt8.getSc()[0],pt8.getSc()[1]);
          comp2D.drawLine(pt8.getSc()[0],pt8.getSc()[1],pt7.getSc()[0],pt7.getSc()[1]);
          comp2D.drawLine(pt7.getSc()[0],pt7.getSc()[1],pt5.getSc()[0],pt5.getSc()[1]);

          comp2D.drawLine(pt1.getSc()[0],pt1.getSc()[1],pt5.getSc()[0],pt5.getSc()[1]);
          comp2D.drawLine(pt2.getSc()[0],pt2.getSc()[1],pt6.getSc()[0],pt6.getSc()[1]);
          comp2D.drawLine(pt3.getSc()[0],pt3.getSc()[1],pt7.getSc()[0],pt7.getSc()[1]);
          comp2D.drawLine(pt4.getSc()[0],pt4.getSc()[1],pt8.getSc()[0],pt8.getSc()[1]);
        }
      break;
      case 7:
        comp2D.setColor(color);
        for (int i=0;i<points.length-1;i++) {
        }
      break;
    }
  }
}
