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

package jat.matvec.data.arrayTools;

   /** Find a value or a verified condition in a 2D-Array of double.
   */

public class Find {

/* ------------------------
   Class variables
 * ------------------------ */

   /** 2D-Array to test.
   */
  private double[][] A;

  /** Indices of elements found.
   */
  private int[][] indices;

   /** Value to find or to compare.
   */
  private double value;

   /** String of the test.
   */
  private String test;


/* ------------------------
   Constructors
 * ------------------------ */

   /** Find a value.
   @param a    Double array to test.
   @param v    Value to compare.
   */

  public Find(double[][] a, double v) {
    A = a;
    value = v;
    test = "==";
    indices =  find(A,test,value);
  }

   /** Find elements verifying a test.
   @param a    Double array to test.
   @param t    Test : "==", "<=", ">=", "<", ">", "!=".
   @param v    Value to compare.
   */


  public Find(double[][] a, String t, double v) {
    A = a;
    value = v;
    test = t;
    indices =  find(A,test,value);
  }


/* ------------------------
   Public Methods
 * ------------------------ */

  /** Get the indices verifying the test.
   @return  2D-indices.
   */

  public int[][] getIndices() {
    return indices;
  }

   /** Get the boolean array of the test (true / false).
   @return  2D-boolean array.
   */

  public boolean[][] getBooleanArray() {
    boolean[][] b = setBooleanArray(indices,A);
    return b;
  }

   /** Get the double array of the test (1D / 0D).
   @return  2D-double array of 1 or 0.
   */

  public double[][] getDoubleArray() {
    double[][] d = setDoubleArray(indices,A);
    return d;
  }

/* ------------------------
   Private Methods
 * ------------------------ */

   /** Set the boolean array.
   @param ind    Indices verifying the test.
   @param a    Array to test
   @return  Boolean array.
   */

  private boolean[][] setBooleanArray(int[][] ind, double[][] a) {
    boolean[][] b = new boolean[a.length][a[0].length];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        b[i][j] = false;
      }
    }
    for (int i = 0; i < ind.length; i++) {
      b[ind[i][0]][ind[i][1]] = true;
    }
    return b;
  }

     /** Set the double array.
   @param ind    Indices verifying the test.
   @param a    Array to test.
   @return  Double array.
   */

  private double[][] setDoubleArray(int[][] ind, double[][] a) {
    double[][] d = new double[a.length][a[0].length];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        d[i][j] = 0;
      }
    }
    for (int i = 0; i < ind.length; i++) {
      d[ind[i][0]][ind[i][1]] = 1;
    }
    return d;
  }

    /** Find elements verifying the test.
   @param a    Double array to test.
   @param t    String of the test.
   @param v    Value to test
   @return  Double array.
   */

  private int[][] find(double[][] a, String t, double v) {
    if (t.equals("==")) {
      return findEqual(a,v);
    } else if (t.equals("<=")) {
      return findInfEqual(a,v);
    } else if (t.equals(">=")) {
      return findSupEqual(a,v);
    } else if (t.equals("<")) {
      return findInf(a,v);
    } else if (t.equals(">")) {
      return findSup(a,v);
    } else if (t.equals("!=")) {
      return findDiff(a,v);
    } else {
      throw new IllegalArgumentException("Test String " + t + " is unknown.");
    }
  }


  private int[][] findEqual(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] == v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] findInfEqual(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] <= v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] findSupEqual(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] >= v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] findInf(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] < v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] findSup(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] > v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] findDiff(double[][] a, double v) {
    int[][] ind = new int[0][2];
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a[0].length; j++) {
        if (a[i][j] != v)
          ind = put(ind,i,j);
      }
    }
    return ind;
  }

  private int[][] put(int[][] ind, int i0, int j0) {
    int[][] new_ind = new int[ind.length+1][2];
    for (int i = 0; i < ind.length; i++) {
      for (int j = 0; j < 2; j++) {
        new_ind[i][j] = ind[i][j];
      }
    }
    new_ind[ind.length][0] = i0;
    new_ind[ind.length][1] = j0;
    return new_ind;
  }

}