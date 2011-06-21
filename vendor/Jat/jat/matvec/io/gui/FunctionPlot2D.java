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

public class FunctionPlot2D extends Plot2D {

  private DoubleFunction[] F;

  private double Xmin;
  private double Xmax;

  private static int nbPointsX = 100;

  public FunctionPlot2D(DoubleFunction f, double xmin, double xmax) {
    setAppearence();
    setPlotAttributes();
    Xmin = xmin;
    Xmax = xmax;
    update(f);
  }

  public FunctionPlot2D(DoubleFunction[] f, double xmin, double xmax) {
    setAppearence();
    setPlotAttributes();
    Xmin = xmin;
    Xmax = xmax;

    update(f);
  }

  protected void setPlotAttributes() {
    PA = new PlotAttributes();
    PA.setTypeList(LINE);
    String[] leg = {"X","Y"};
    PA.setLegend(leg);
  }

  public void update(DoubleFunction f) {
    checkArgNumber(f);
    F = new DoubleFunction[1];
    F[0] = f;

    setXY();
    update();
  }

  public void update(DoubleFunction[] f) {
    checkArgNumber(f);
    F = new DoubleFunction[f.length];
    for (int i = 0; i < f.length; i++) {
      F[i] = f[i];
    }

    setXY();
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

    setXY();
    update();
  }

  public void setMinMax(double xmin, double xmax) {
    Xmin = xmin;
    Xmax = xmax;

    setXY();
    update();
  }

  private void setXY() {
    X = new double[F.length][nbPointsX];
    Y = new double[F.length][nbPointsX];
    widthX = new double[F.length][];
    widthY = new double[F.length][];
    for (int i = 0; i < F.length; i++) {
      for (int j = 0; j < nbPointsX; j++) {
        double[] x = {Xmin + (Xmax - Xmin)*j/(nbPointsX-1)};
        X[i][j] = x[0];
        Y[i][j] = F[i].eval(x);
      }
      widthX[i] = new double[nbPointsX];
      widthY[i] = new double[nbPointsX];
    }
  }

  /** Check if argNumber == 1.
  @param f   DoubleFunction.
   */

   private void checkArgNumber (DoubleFunction f) {
      f.checkArgNumber(1);
   }

  /** Check if argNumber == 1.
  @param F   DoubleFunction array.
   */

   private void checkArgNumber (DoubleFunction[] F) {
    for (int i = 0; i < F.length; i++)
      F[i].checkArgNumber(1);
   }


}
