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

package jat.matvec.function;

import jat.matvec.data.Matrix;

public class DoubleFunctionInterpolation extends DoubleFunction {

  //private int argNumber;

  private double[][] X;
  private double[] Y;

  public DoubleFunctionInterpolation(double[][] in, double[] out) {
    if (in.length != out.length) {
      throw new IllegalArgumentException("p1 = double[points][coordinates] and p2 = double[points]");
    }
    argNumber = in[0].length;
    X = in;
    Y = out;
  }

  public double eval(double[] values) {
    checkArgNumber(values.length);
    return interpolate(values);
  }

  public double eval(double value) {
    checkArgNumber(1);
    double[] values = new double[1];
    values[0] = value;
    return interpolate(values);
  }

  private double interpolate(double[] x) {
    int[] close = closest(x,argNumber+1);
    Matrix M_X = new Matrix(X).getRows(close).transpose();
    Matrix M_Y = new Matrix(Y,Y.length).getRows(close);
    Matrix M_x = new Matrix(x,x.length);

    Matrix A = new Matrix(1,argNumber+1,1).merge(M_X);
    Matrix B = new Matrix(1,1,1).merge(M_x);
    Matrix W = A.solve(B);
    Matrix y = W.transpose().times(M_Y);

    return y.get(0,0);
  }

  private int[] closest(double[] x, int num) {
    double[] d = new double[X.length];
    for (int i = 0; i < d.length; i++) {
      d[i] = dist(x,X[i]);
    }
    int[] c = new int[num];
    double dMin;
    for (int j = 0; j < num; j++) {
    dMin = Double.MAX_VALUE;
      for (int i = 0; i < X.length; i++) {
        if ((d[i] < dMin)) {
          dMin = d[i];
          c[j] = i;
        }
      }
      d[c[j]] = Double.MAX_VALUE;
    }
    return c;
  }

  private double dist(double[] x, double[] y) {
    double d = 0;
    for (int i = 0; i < x.length; i++) {
       d = d + (x[i] - y[i])*(x[i] - y[i]);
    }
    return /*Math.sqrt(*/d/*)*/;
  }

}
