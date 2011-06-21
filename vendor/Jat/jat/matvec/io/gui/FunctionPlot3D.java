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

import jat.matvec.function.DoubleFunction;

import jat.matvec.io.gui.plotTools.*;

public class FunctionPlot3D extends Plot3D {

  private DoubleFunction[] F;

  private double Xmin;
  private double Xmax;
  private double Ymin;
  private double Ymax;

  private static int nbPointsX = 20;
  private static int nbPointsY = 20;

  public FunctionPlot3D(DoubleFunction f, double xmin, double xmax, double ymin, double ymax) {
    setAppearence();
    setPlotAttributes();
    Xmin = xmin;
    Xmax = xmax;
    Ymin = ymin;
    Ymax = ymax;
    update(f);
  }

  public FunctionPlot3D(DoubleFunction[] f, double xmin, double xmax, double ymin, double ymax) {
    setAppearence();
    setPlotAttributes();
    Xmin = xmin;
    Xmax = xmax;
    Ymin = ymin;
    Ymax = ymax;
    update(f);
  }

  protected void setPlotAttributes() {
    PA = new PlotAttributes();
    PA.setTypeList(PIXEL);
    String[] leg = {"X","Y","Z"};
    PA.setLegend(leg);
  }

  public void update(DoubleFunction f) {
    checkArgNumber(f);
    F = new DoubleFunction[1];
    F[0] = f;

    setXYZ();
    update();
  }

  public void update(DoubleFunction[] f) {
    checkArgNumber(f);
    F = new DoubleFunction[f.length];
    for (int i = 0; i < f.length; i++) {
      F[i] = f[i];
    }

    setXYZ();
    update();
  }

  public void add(DoubleFunction f) {
    checkArgNumber(f);
    DoubleFunction[] F_tmp = new DoubleFunction[F.length + 1];
    for (int i = 0; i < F.length; i++) {
      F_tmp[i] = F[i];
    }
    F_tmp[F.length] = f;
    F = F_tmp;

    setXYZ();
    update();
  }

  public void setMinMax(double xmin, double xmax, double ymin, double ymax) {
    Xmin = xmin;
    Xmax = xmax;
    Ymin = ymin;
    Ymax = ymax;

    setXYZ();
    update();
  }

  private void setXYZ() {
    X = new double[F.length][nbPointsX*nbPointsY];
    Y = new double[F.length][nbPointsX*nbPointsY];
    Z = new double[F.length][nbPointsX*nbPointsY];
    widthX = new double[F.length][];
    widthY = new double[F.length][];
    widthZ = new double[F.length][];
    for (int i = 0; i < F.length; i++) {
      for (int j = 0; j < nbPointsX; j++) {
        for (int k = 0; k < nbPointsY; k++) {
          double[] xy = {Xmin + (Xmax - Xmin)*j/(nbPointsX-1),Ymin + (Ymax - Ymin)*k/(nbPointsY-1)};
          X[i][j+k*nbPointsX] = xy[0];
          Y[i][j+k*nbPointsX] = xy[1];
          Z[i][j+k*nbPointsX] = F[i].eval(xy);
        }
      }
      widthX[i] = new double[nbPointsX*nbPointsY];
      widthY[i] = new double[nbPointsX*nbPointsY];
      widthZ[i] = new double[nbPointsX*nbPointsY];
    }
  }

  /** Check if argNumber == 2.
  @param f   DoubleFunction.
   */

   private void checkArgNumber (DoubleFunction f) {
      f.checkArgNumber(2);
   }

  /** Check if argNumber == 2.
  @param F   DoubleFunction array.
   */

   private void checkArgNumber (DoubleFunction[] F) {
    for (int i = 0; i < F.length; i++)
      F[i].checkArgNumber(2);
   }


}
