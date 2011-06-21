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

import javax.swing.*;
import java.awt.*;

public class Axe3D {

  private double[] Pl_XMin = new double[3];
  private double[] Pl_XMax = new double[3];
  private double[] Pl_YMin = new double[3];
  private double[] Pl_YMax = new double[3];
  private double[] Pl_ZMin = new double[3];
  private double[] Pl_ZMax = new double[3];

  private int[] Sc_XMin = new int[2];
  private int[] Sc_XMax = new int[2];
  private int[] Sc_YMin = new int[2];
  private int[] Sc_YMax = new int[2];
  private int[] Sc_ZMin = new int[2];
  private int[] Sc_ZMax = new int[2];

  private Eye3D eye;
  private String[] legend;

  private double[] Pl_0 = {0,0,0};
  private int[] Sc_0 = new int[2];

  private Dimension panelDimension;

  private static double delta = 0.2;

  public Axe3D(double[] x,double[] y,double[] z,JPanel panel,Eye3D ey,String[] leg) {
    legend = leg;
    eye = ey;
    Pl_XMin[0] = Min(Min(x),Pl_0[0])-(Max(x)-Min(x))*delta;
    Pl_XMin[1] = Pl_0[1];
    Pl_XMin[2] = Pl_0[2];
    Pl_XMax[0] = Max(Max(x),Pl_0[0])+(Max(x)-Min(x))*delta;
    Pl_XMax[1] = Pl_0[1];
    Pl_XMax[2] = Pl_0[2];
    Pl_YMin[0] = Pl_0[0];
    Pl_YMin[1] = Min(Min(y),Pl_0[1])-(Max(y)-Min(y))*delta;
    Pl_YMin[2] = Pl_0[2];
    Pl_YMax[0] = Pl_0[0];
    Pl_YMax[1] = Max(Max(y),Pl_0[1])+(Max(y)-Min(y))*delta;
    Pl_YMax[2] = Pl_0[2];
    Pl_ZMin[0] = Pl_0[0];
    Pl_ZMin[1] = Pl_0[1];
    Pl_ZMin[2] = Min(Min(z),Pl_0[2])-(Max(z)-Min(z))*delta;
    Pl_ZMax[0] = Pl_0[0];
    Pl_ZMax[1] = Pl_0[1];
    Pl_ZMax[2] = Max(Max(z),Pl_0[2])+(Max(z)-Min(z))*delta;

    panelDimension = panel.getSize();

    PlScConvert();

  }

  public Axe3D(double x0,double y0,double z0,double[] x,double[] y,double[] z,JPanel panel,Eye3D ey,String[] leg) {

    Pl_0[0] = x0;
    Pl_0[1] = y0;
    Pl_0[2] = z0;

    legend = leg;
    eye = ey;
    Pl_XMin[0] = Min(Min(x),Pl_0[0])-(Max(x)-Min(x))*delta;
    Pl_XMin[1] = Pl_0[1];
    Pl_XMin[2] = Pl_0[2];
    Pl_XMax[0] = Max(Max(x),Pl_0[0])+(Max(x)-Min(x))*delta;
    Pl_XMax[1] = Pl_0[1];
    Pl_XMax[2] = Pl_0[2];
    Pl_YMin[0] = Pl_0[0];
    Pl_YMin[1] = Min(Min(y),Pl_0[1])-(Max(y)-Min(y))*delta;
    Pl_YMin[2] = Pl_0[2];
    Pl_YMax[0] = Pl_0[0];
    Pl_YMax[1] = Max(Max(y),Pl_0[1])+(Max(y)-Min(y))*delta;
    Pl_YMax[2] = Pl_0[2];
    Pl_ZMin[0] = Pl_0[0];
    Pl_ZMin[1] = Pl_0[1];
    Pl_ZMin[2] = Min(Min(z),Pl_0[2])-(Max(z)-Min(z))*delta;
    Pl_ZMax[0] = Pl_0[0];
    Pl_ZMax[1] = Pl_0[1];
    Pl_ZMax[2] = Max(Max(z),Pl_0[2])+(Max(z)-Min(z))*delta;

    panelDimension = panel.getSize();

    PlScConvert();

  }

  private void PlScConvert() {

    Sc_0 = Pl2Sc(Pl_0);

    Sc_XMin = Pl2Sc(Pl_XMin);
    Sc_XMax = Pl2Sc(Pl_XMax);
    Sc_YMin = Pl2Sc(Pl_YMin);
    Sc_YMax = Pl2Sc(Pl_YMax);
    Sc_ZMin = Pl2Sc(Pl_ZMin);
    Sc_ZMax = Pl2Sc(Pl_ZMax);
  }

  private double[] D2D(double[] xyz) {
    double[] xy = new double[2];
    xy[0] = Math.cos(eye.getTheta())*((xyz[1]-(Pl_YMax[1]+Pl_YMin[1])/2)/(Pl_YMax[1]-Pl_YMin[1])) - Math.sin(eye.getTheta())*((xyz[0]-(Pl_XMax[0]+Pl_XMin[0])/2)/(Pl_XMax[0]-Pl_XMin[0]));
    xy[1] = Math.cos(eye.getPhi())*((xyz[2]-(Pl_ZMax[2]+Pl_ZMin[2])/2)/(Pl_ZMax[2]-Pl_ZMin[2])) - Math.sin(eye.getPhi())*Math.cos(eye.getTheta())*((xyz[0]-(Pl_XMax[0]+Pl_XMin[0])/2)/(Pl_XMax[0]-Pl_XMin[0])) - Math.sin(eye.getPhi())*Math.sin(eye.getTheta())*((xyz[1]-(Pl_YMax[1]+Pl_YMin[1])/2)/(Pl_YMax[1]-Pl_YMin[1]));
    return xy;
  }

  public int[][] Pl2Sc(double[][] XY) {
    int[][] ret = new int[XY.length][2];
    for (int i = 0;i<XY.length;i++) {
      double[] XY_temp = new double[3];
      XY_temp[0] = XY[i][0];
      XY_temp[1] = XY[i][1];
      XY_temp[2] = XY[i][2];
      int[] temp = Pl2Sc(XY_temp);
      ret[i][0] = temp[0];
      ret[i][1] = temp[1];
    }
    return ret;
  }

   public int[][][] Pl2Sc(double[][][] XY) {
    int[][][] ret = new int[XY.length][2][5];
    for (int i = 0;i<XY.length;i++) {
      for (int j=0;j<5;j++) {
      double[] XY_temp = new double[3];
      XY_temp[0] = XY[i][0][j];
      XY_temp[1] = XY[i][1][j];
      XY_temp[2] = XY[i][2][j];
      int[] temp = Pl2Sc(XY_temp);
      ret[i][0][j] = temp[0];
      ret[i][1][j] = temp[1];
      }
    }
    return ret;
  }

  public int[] Pl2Sc(double[] xyz) {

    int h = (int)panelDimension.getHeight();
    int w = (int)panelDimension.getWidth();

    double[] xy = D2D(xyz);

    int[] temp = new int[2];
    temp[0] = (int)(w/2)+(int)(((double)w)*(xy[0]/1.7/*/(Pl_XMax[0]-Pl_XMin[0])*/));
    temp[1] = (int)(h/2)-(int)(((double)h)*(xy[1]/1.7/*/(Pl_YMax[1]-Pl_YMin[1])*/));
    return temp;
  }

  /*public double[] Sc2Pl(int[] xy) {
  }*/

  public void setEye3D(Eye3D ey) {
    eye = ey;
  }

  public Eye3D getEye3D() {
    return (Eye3D)(eye.clone());
  }

  public void draw(Graphics2D comp2D) {
    PlScConvert();
    comp2D.setColor(Color.black);
    comp2D.drawLine(Sc_XMin[0],Sc_XMin[1],Sc_XMax[0],Sc_XMax[1]);
    comp2D.drawLine(Sc_YMin[0],Sc_YMin[1],Sc_YMax[0],Sc_YMax[1]);
    comp2D.drawLine(Sc_ZMin[0],Sc_ZMin[1],Sc_ZMax[0],Sc_ZMax[1]);
    comp2D.setFont(new Font("Arial",Font.BOLD,14));
    comp2D.drawString(legend[0],(int)(0.9*Sc_XMax[0]+0.1*Sc_XMin[0]),(int)(0.9*Sc_XMax[1]+0.1*Sc_XMin[1]));
    comp2D.drawString(legend[1],(int)(0.9*Sc_YMax[0]+0.1*Sc_YMin[0]),(int)(0.9*Sc_YMax[1]+0.1*Sc_YMin[1]));
    comp2D.drawString(legend[2],(int)(0.9*Sc_ZMax[0]+0.1*Sc_ZMin[0]),(int)(0.9*Sc_ZMax[1]+0.1*Sc_ZMin[1]));
  }

  public int[] getSc0() {
    return Sc_0;
  }

  public int[] getScXMin() {
    return Sc_XMin;
  }

  public int[] getScXMax() {
    return Sc_XMax;
  }

  public int[] getScYMin() {
    return Sc_YMin;
  }

  public int[] getScYMax() {
    return Sc_YMax;
  }

  public int[] getScZMin() {
    return Sc_ZMin;
  }

  public int[] getScZMax() {
    return Sc_ZMax;
  }

  public double[] getPlXMin() {
    return Pl_XMin;
  }

  public double[] getPl0() {
    return Pl_0;
  }

  public double[] getPlXMax() {
    return Pl_XMax;
  }

  public double[] getPlYMin() {
    return Pl_YMin;
  }

  public double[] getPlYMax() {
    return Pl_YMax;
  }

  public double[] getPlZMin() {
    return Pl_ZMin;
  }

  public double[] getPlZMax() {
    return Pl_ZMax;
  }

 private double Min(double[] list) {
    double temp = list[0];
    for (int i=0;i<list.length;i++) {
      temp = Math.min(temp,list[i]);
    }
    return temp;
  }

  private double Min(double a, double b) {
    return Math.min(a,b);
  }

  private double Max(double[] list) {
    double temp = list[0];
    for (int i=0;i<list.length;i++) {
      temp = Math.max(temp,list[i]);
    }
    return temp;
  }

  private double Max(double a,double b) {
    return Math.max(a,b);
  }

}
