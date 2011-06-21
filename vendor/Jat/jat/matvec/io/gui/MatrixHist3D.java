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
import jat.matvec.data.RandomMatrix;

import jat.matvec.io.gui.plotTools.*;

public class MatrixHist3D extends Plot3D {

  private Matrix[] XYZ;
  private Matrix[] width;

  private int nbSlicesX;
  private int nbSlicesY;

  public MatrixHist3D(RandomMatrix xy,int n, int l) {
    setAppearence();
    setPlotAttributes();
    nbSlicesX = n;
    nbSlicesY = l;
    update(xy);
  }

  public MatrixHist3D(RandomMatrix[] xy,int n, int l) {
    setAppearence();
    setPlotAttributes();
    nbSlicesX = n;
    nbSlicesY = l;
    update(xy);
  }

  protected void setPlotAttributes() {
    PA = new PlotAttributes();
    PA.setTypeList(HIST);
    String[] leg = {"x","y",""};
    PA.setLegend(leg);
  }

  private void setSlicing(int nbx, int nby) {

    X0 = XYZ[0].sumRows().get(0,0)/XYZ[0].getRowDimension();
    Y0 = XYZ[0].sumRows().get(0,1)/XYZ[0].getRowDimension();

    nbSlicesX = nbx;
    nbSlicesY = nby;
    int n = XYZ.length;
    int m;

    Matrix[] mat = new Matrix[n];
    Matrix temp;

    width = new Matrix[n];

    double colminX;
    double colmaxX;
    double colpasX;
    double colminY;
    double colmaxY;
    double colpasY;

    for (int i = 0; i < n; i++) {
      temp = new Matrix(nbSlicesX*nbSlicesY,3);
      m = XYZ[i].getRowDimension();
      colminX = XYZ[i].minRows().get(0,0);
      colmaxX = XYZ[i].maxRows().get(0,0);
      colpasX = (colmaxX-colminX)/(nbSlicesX);
      colminY = XYZ[i].minRows().get(0,1);
      colmaxY = XYZ[i].maxRows().get(0,1);
      colpasY = (colmaxY-colminY)/(nbSlicesY);

      width[i] = new Matrix(nbSlicesX*nbSlicesY,2);

      for (int j = 0; j < nbSlicesX; j++) {
        for (int l = 0; l < nbSlicesY; l++) {
          temp.set(j*nbSlicesY+l,0,(j+.5)*colpasX+colminX);
          temp.set(j*nbSlicesY+l,1,(l+.5)*colpasY+colminY);
          width[i].set(j*nbSlicesY+l,0,colpasX);
          width[i].set(j*nbSlicesY+l,1,colpasY);
        }
      }

      for (int k = 0; k < m; k++) {
        int[] ind = {0,1};
        Matrix d = temp.getColumns(ind).dist(XYZ[i].getMatrix(k,k,0,1));
        int s = d.find(d.min().get(0,0))[0][0];
        temp.set(s,2,temp.get(s,2)+1);
      }
      mat[i] = temp.copy();
    }
    XYZ = mat;
  }


  public void update(RandomMatrix xy) {
    checkColumnDimension(xy);
    XYZ = new Matrix[1];
    XYZ[0] = xy.copy();

    setXYZ();
    update();
  }

  public void update(RandomMatrix[] xy) {
    checkColumnDimension(xy);
    XYZ = new Matrix[xy.length];
    for (int i = 0; i < xy.length; i++) {
      XYZ[i] = xy[i].copy();
    }

    setXYZ();
    update();
  }

  public void add(Matrix xy) {
    checkColumnDimension(xy);
    Matrix[] XYZ_tmp = new Matrix[XYZ.length + 1];
    for (int i = 0; i < XYZ.length; i++) {
      XYZ_tmp[i] = XYZ[i];
    }
    XYZ_tmp[XYZ.length] = xy.copy();
    XYZ = XYZ_tmp;

    setXYZ();
    update();
  }

  public void setNumberSlices(int nx,int ny) {
    nbSlicesX = nx;
    nbSlicesY = ny;

    setXYZ();
    update();
  }

  private void TransposeIfNecessary() {
    for (int i = 0; i < XYZ.length; i++) {
      if (XYZ[i].getRowDimension()<XYZ[i].getColumnDimension()) {
        XYZ[i] = XYZ[i].transpose();
      }
    }
  }

  private void setXYZ() {
    TransposeIfNecessary();
    setSlicing(nbSlicesX,nbSlicesY);
    X = new double[XYZ.length][];
    Y = new double[XYZ.length][];
    Z = new double[XYZ.length][];
    widthX = new double[XYZ.length][];
    widthY = new double[XYZ.length][];
    widthZ = new double[XYZ.length][];
    for (int i = 0; i < XYZ.length; i++) {
      X[i] = XYZ[i].getColumnArrayCopy(0);
      Y[i] = XYZ[i].getColumnArrayCopy(1);
      Z[i] = XYZ[i].getColumnArrayCopy(2);
      widthX[i] = width[i].getColumnArrayCopy(0);
      widthY[i] = width[i].getColumnArrayCopy(1);
      widthZ[i] = new double[XYZ[i].getColumnDimension()];
    }
  }

  /** Check if ColumnDimension(xy) == 2
  @param xy   Matrix
   */

   private void checkColumnDimension (Matrix xy) {
      xy.checkColumnDimension(2);
   }

  /** Check if ColumnDimension(xxy) == 2
  @param xy   Matrix
   */

   private void checkColumnDimension (Matrix[] xy) {
    for (int i = 0; i < xy.length; i++)
      xy[i].checkColumnDimension(2);
   }
}
