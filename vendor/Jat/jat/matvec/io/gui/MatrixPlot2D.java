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

import jat.matvec.data.Matrix;

import jat.matvec.io.gui.plotTools.*;

public class MatrixPlot2D extends Plot2D {

  private Matrix[] XY;

  public MatrixPlot2D(Matrix xy) {
    setAppearence();
    setPlotAttributes();
    update(xy);
  }

  public MatrixPlot2D(Matrix[] xy) {
    setAppearence();
    setPlotAttributes();
    update(xy);
  }

  public MatrixPlot2D(Matrix x,Matrix y) {
    setAppearence();
    setPlotAttributes();
    update(x,y);
  }

  protected void setPlotAttributes() {
    PA = new PlotAttributes();
    PA.setTypeList(DOT);
    String[] leg = {"X","Y"};
    PA.setLegend(leg);
  }

  public void update(Matrix xy) {
    checkColumnDimension(xy);
    XY = new Matrix[1];
    XY[0] = xy.copy();

    setXY();
    update();
  }

  public void update(Matrix[] xy) {
    checkColumnDimension(xy);
    XY = new Matrix[xy.length];
    for (int i = 0; i <xy.length; i++) {
      XY[i] = xy[i].copy();
    }

    setXY();
    update();
  }

  public void add(Matrix xy) {
    checkColumnDimension(xy);
    Matrix[] XY_tmp = new Matrix[XY.length + 1];
    for (int i = 0; i < XY.length; i++) {
      XY_tmp[i] = XY[i];
    }
    XY_tmp[XY.length] = xy.copy();
    XY = XY_tmp;

    setXY();
    update();
  }

  public void update(Matrix x,Matrix y) {
    checkDimensions(x,y);
    XY = new Matrix[1];
    Matrix xy = new Matrix(x.getRowDimension(),2);
    xy.setMatrix(0,0,x.copy());
    xy.setMatrix(0,1,y.copy());
    XY[0] = xy;

    setXY();
    update();
  }

  private void TransposeIfNecessary() {
    for (int i = 0; i < XY.length; i++) {
      if (XY[i].getRowDimension()<XY[i].getColumnDimension()) {
        XY[i] = XY[i].transpose();
      }
    }
  }

  private void setXY() {
    TransposeIfNecessary();
    X = new double[XY.length][];
    Y = new double[XY.length][];
    widthX = new double[XY.length][];
    widthY = new double[XY.length][];
    for (int i = 0; i < XY.length; i++) {
      X[i] = XY[i].getColumnArrayCopy(0);
      Y[i] = XY[i].getColumnArrayCopy(1);
      widthX[i] = new double[XY[i].getColumnDimension()];
      widthY[i] = new double[XY[i].getColumnDimension()];
    }
  }


  /** Check if ColumnDimension(xy) == 2
  @param xy   Matrix
   */

   private void checkColumnDimension (Matrix xy) {
      xy.checkColumnDimension(2);
   }

  /** Check if ColumnDimension(xy) == 2
  @param xy   Matrix
   */

   private void checkColumnDimension (Matrix[] xy) {
    for (int i = 0; i < xy.length; i++)
      xy[i].checkColumnDimension(2);
   }

  /** Check if size(x) == size(y)
  @param x   Matrix
  @param y   Matrix
   */

   private void checkDimensions(Matrix x,Matrix y) {
      x.checkColumnDimension(1);
      x.checkMatrixDimensions(y);
   }

}
