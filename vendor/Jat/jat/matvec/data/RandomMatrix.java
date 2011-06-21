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

package jat.matvec.data;

import jat.matvec.data.Matrix;
import jat.matvec.data.RandomVariable;

import jat.matvec.io.gui.FrameView;
import jat.matvec.io.gui.MatrixHist2D;
import jat.matvec.io.gui.MatrixHist3D;

import jat.matvec.function.DoubleFunction;
import jat.matvec.data.arrayTools.Shuffle;

/**
<P>
   The RandomMatrix Class provides tools for statistical simulations,it extends the Matrix Class and adds many methods.

@author Yann RICHET.
@version 2.0
*/


public class RandomMatrix extends Matrix {

	private static final long serialVersionUID = 1L;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Is the RandomMatrix a sample or an overal population?
   */

    private boolean isSample = true;


/* ------------------------
   Constructors
 * ------------------------ */

   /** Construct an m-by-n matrix of 0.
   @param m    Number of rows.
   @param n    Number of columns.
   */

  public RandomMatrix(int m,int n) {
    super(m,n);
  }

   /** Construct a RandomMatrix from a Matrix (conversion).
   @param M    Matrix.
   */

  public RandomMatrix(Matrix M) {
    super(M.getArrayCopy());
  }

/* ----------------------
   Public Methods
 * ---------------------- */

////////////////////////
//Static constructors.//
////////////////////////

  /** Construct an m-by-n matrix of random numbers from a uniform random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix uniform(int m,int n,double min,double max) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.uniform(min, max);
        C[i][j] = x;
      }
    }
    return X;
  }

  /** Construct an m-by-n matrix of random numbers from a discrete random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param val_prob    Matrix of the discrete value and their probabilities.
   @return      A RandomMatrix.
   */

  public static RandomMatrix dirac(int m,int n,Matrix val_prob) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    double y;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.dirac(val_prob.getColumnArrayCopy(0),val_prob.getColumnArrayCopy(1));
        C[i][j] = x;
      }
    }
    return X;
  }

   /** Construct an m-by-n matrix of random numbers from a Gaussian (Normal) random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param mu    Mean of the random variable.
   @param sigma    Standard deviation of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix normal(int m,int n,double mu,double sigma) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.normal(mu,sigma);
        C[i][j] = x;
      }
    }
    return X;
  }

   /** Construct an m-by-n matrix of random numbers from a LogNormal random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param mu    Mean of the Normal random variable.
   @param sigma    Standard deviation of the Normal random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix logNormal(int m,int n,double mu,double sigma) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.logNormal(mu,sigma);
        C[i][j] = Math.exp(x);
      }
    }
    return X;
  }

   /** Construct an m-by-n matrix of random numbers from an exponantial random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param lambda    Parmaeter of the exponential random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix exponential(int m,int n,double lambda) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.exponential(lambda);
        C[i][j] = x;
      }
    }
    return X;
  }


   /** Construct an m-by-n matrix of random numbers from a symetric triangular random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix triangular(int m,int n,double min,double max) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.triangular(min,max);
        C[i][j] = x;
      }
    }
    return X;
  }

     /** Construct an m-by-n matrix of random numbers from a non-symetric triangular random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param min    Min of the random variable.
   @param med    Value of the random variable with max density.
   @param max    Max of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix triangular(int m,int n,double min,double med,double max) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    double y;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.triangular(min,med,max);
        C[i][j] = x;
      }
    }
    return X;
  }

     /** Construct an m-by-n matrix of random numbers from a Beta random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param a    First parameter of the Beta random variable.
   @param b    Second parameter of the Beta random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix beta(int m,int n,double a,double b) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    double y;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.beta(a,b);
        C[i][j] = x;
      }
    }
    return X;
  }

     /** Construct an m-by-n matrix of random numbers from a Cauchy random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param mu    Median of the Weibull random variable
   @param sigma    Second parameter of the Cauchy random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix cauchy(int m,int n,double mu,double sigma) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    double y;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.cauchy(mu,sigma);
        C[i][j] = x;
      }
    }
    return X;
  }

     /** Construct an m-by-n matrix of random numbers from a Weibull random variable.
   @param m    Number of rows.
   @param n    Number of columns.
   @param lambda    First parameter of the Weibull random variable.
   @param c    Second parameter of the Weibull random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix weibull(int m,int n,double lambda,double c) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    double y;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        x = RandomVariable.weibull(lambda,c);
        C[i][j] = x;
      }
    }
    return X;
  }

  /** Construct an m-by-n matrix of random numbers from a sample of this random variable, using the BOOTSTRAP technic.
   @param m    Number of rows.
   @param n    Number of columns.
   @param sample    RandomMatrix of a sample of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix bootstrap(int m,int n,RandomMatrix sample) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        int num = RandomVariable.randInt(0,sample.n*sample.m-1);
        x = sample.getColumnPackedCopy()[num];
        C[i][j] = x;

      }
    }
    return X;
  }

  /** Construct an m-by-n matrix of random numbers from a random variable definied by its density function, using the rejection technic.
   *  ! WARNING : this simulation technic can take a very long time !
   @param m    Number of rows.
   @param n    Number of columns.
   @param fun    Density function of the random variable.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A RandomMatrix.
   */

  public static RandomMatrix rejection(int m,int n,DoubleFunction fun, double min, double max) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();

    double maxFun = 0;
    for (int i = 0; i < (10*m*n); i++) {
      double[] val = {min + i*(max-min)/(10*m*n-1)};
      maxFun = Math.max(maxFun,fun.eval(val));
    }

    double try_x;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        try_x = RandomVariable.rejection(fun,maxFun,min,max);
        C[i][j] = try_x;
      }
    }
    return X;
  }

   /** Construct shuffled copy of a matrix.
   @param B    Matrix to shuffle.
   @return      A RandomMatrix.
   */

  public static RandomMatrix shuffle(Matrix B) {
    Shuffle order = new Shuffle(B.m*B.n);
    double[] b = B.getColumnPackedCopy();
    double[] x = new double[b.length];
    for (int i = 0; i < x.length; i++) {
      x[i] = b[order.getOrder(i)];
    }
    RandomMatrix X = new RandomMatrix(new Matrix(x,B.m));
    return X;
  }

   /** Construct a sample with replacement of a matrix.
   @param m    Number of rows.
   @param n    Number of columns.
   @param B    Matrix to sample.
   @return      A RandomMatrix.
   */

  public static RandomMatrix sampleWithReplacement(int m, int n,Matrix B) {
    RandomMatrix X = new RandomMatrix(m,n);
    double[][] C = X.getArray();
    double[] b = B.getColumnPackedCopy();
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        int num = RandomVariable.randInt(0,B.n*B.m-1);
        C[i][j] = b[num];
      }
    }
    return X;
  }

   /** Construct a sample without replacement of a matrix.
   @param m    Number of rows.
   @param n    Number of columns.
   @param B    Matrix to sample.
   @return      A RandomMatrix.
   */

  public static RandomMatrix sampleWithoutReplacement(int m, int n,Matrix B) {
    if ((m*n) > ((B.m)*(B.n))) {
      throw new IllegalArgumentException("Matrix number of elements must be < " + ((B.m)*(B.n)));
    }
    Shuffle order = new Shuffle(B.m*B.n);
    double[] b = B.getColumnPackedCopy();
    double[] x = new double[m*n];
    for (int i = 0; i < x.length; i++) {
      x[i] = b[order.getOrder(i)];
    }
    RandomMatrix X = new RandomMatrix(new Matrix(x,m));
    return X;
  }


////////////////////////////////////////
//Norms and characteristics of Matrix.//
////////////////////////////////////////

   /** Specify if the RandomMatrix is a sample of an overall population, or if it's an overall population.
    *  This information is needed to calculate unbiaised estimtors of variance for instance.
   @param is    Is sample?.
   */

    public void setIsSample(boolean is) {
      isSample = is;
    }

   /** Generate a row matrix, each column contents the mean value of the columns.
   @return     An 1-by-n matrix.
   */

   public Matrix mean() {
      return meanRows();
   }

   /** Generate a row matrix, each column contents the mean value of the columns.
   @return     An 1-by-n matrix.
   */

   public Matrix meanRows() {
      Matrix X = new Matrix(1,n);
      double[][] C = X.getArray();
      double s ;
      for (int j = 0; j < n; j++) {
          s = 0;
          for (int k = 0; k < m; k++) {
            s = s + A[k][j];
          }
          C[0][j] = s/m;
      }
      return X;
   }

   /** Generate a column matrix, each line contents the mean value of the lines.
   @return     An m-by-1 matrix.
   */

   public Matrix meanColumns() {
      Matrix X = new Matrix(m,1);
      double[][] C = X.getArray();
      double s ;
      for (int i = 0; i < m; i++) {
          s = 0;
          for (int k = 0; k < n; k++) {
            s = s + A[i][k];
          }
          C[i][0] = s/n;
      }
      return X;
   }

   /** Generate a covariance matrix, each column contains values of a pulling.
   @return     An n-by-n matrix.
   */

   public Matrix cov() {
      return covRows();
   }

   /** Generate a covariance matrix, each column contains values of a pulling.
   @return     An n-by-n matrix.
   */

   public Matrix covRows() {
      Matrix X = new Matrix(n,n);
      int degrees = (isSample) ? (m-1) : (m);
      double[][] C = X.getArray();
      double c ;
      double s1 ;
      double s2 ;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          c = 0;
          s1 = 0;
          s2 = 0;
          for (int k = 0; k < m; k++) {
            s1 = s1 + A[k][i];
            s2 = s2 + A[k][j];
          }
          s1 = s1/m;
          s2 = s2/m;
          for (int k = 0; k < m; k++) {
            c = c + (A[k][i] - s1)*(A[k][j] - s2);
          }
          C[i][j] = c/degrees;
        }
      }
      return X;
   }

   /** Generate a covariance matrix, each row contains values of a pulling.
   @return     An m-by-m matrix.
   */

   public Matrix covColumns() {
      Matrix X = new Matrix(m,m);
      int degrees = (isSample) ? (n-1) : (n);
      double[][] C = X.getArray();
      double c ;
      double s1 ;
      double s2 ;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          c = 0;
          s1 = 0;
          s2 = 0;
          for (int k = 0; k < n; k++) {
            s1 = s1 + A[i][k];
            s2 = s2 + A[j][k];
          }
          s1 = s1/n;
          s2 = s2/n;
          for (int k = 0; k < n; k++) {
            c = c + (A[i][k] - s1)*(A[j][k] - s2);
          }
          C[i][j] = c/degrees;
        }
      }
      return X;
   }

   /** Generate a correlation matrix, each column contains values of a pulling.
   @return     An n-by-n matrix.
   */

   public Matrix cor() {
      return corRows();
   }

   /** Generate a correlation matrix, each column contains values of a pulling.
   @return     An n-by-n matrix.
   */

   public Matrix corRows() {
      Matrix X = new Matrix(n,n);
      int degrees = (isSample) ? (m-1) : (m);
      double[][] C = X.getArray();
      double[][] V = new double[n][n];
      double c ;
      double s1 ;
      double s2 ;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          c = 0;
          s1 = 0;
          s2 = 0;
          for (int k = 0; k < m; k++) {
            s1 = s1 + A[k][i];
            s2 = s2 + A[k][j];
          }
          s1 = s1/m;
          s2 = s2/m;
          for (int k = 0; k < m; k++) {
            c = c + (A[k][i] - s1)*(A[k][j] - s2);
          }
          V[i][j] = c/degrees;
        }
      }
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          C[i][j] = V[i][j]/Math.sqrt(V[i][i]*V[j][j]);//c/degrees;
        }
      }
      return X;
   }

   /** Generate a correlation matrix, each row contains values of a pulling.
   @return     An m-by-m matrix.
   */

   public Matrix corColumns() {
      Matrix X = new Matrix(m,m);
      int degrees = (isSample) ? (n-1) : (n);
      double[][] C = X.getArray();
      double[][] V = new double[m][m];
      double c ;
      double s1 ;
      double s2 ;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          c = 0;
          s1 = 0;
          s2 = 0;
          for (int k = 0; k < n; k++) {
            s1 = s1 + A[i][k];
            s2 = s2 + A[j][k];
          }
          s1 = s1/n;
          s2 = s2/n;
          for (int k = 0; k < n; k++) {
            c = c + (A[i][k] - s1)*(A[j][k] - s2);
          }
          V[i][j] = c/degrees;
        }
      }
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          C[i][j] = V[i][j]/Math.sqrt(V[i][i]*V[j][j]);
        }
      }
      return X;
   }

   /** Generate a variance matrix, each column contains values of a pulling.
   @return     An 1-by-n matrix.
   */

   public Matrix var() {
      return varRows();
   }

   /** Generate a variance matrix, each column contains values of a pulling.
   @return     An 1-by-n matrix.
   */

   public Matrix varRows() {
      Matrix X = new Matrix(1,n);
      int degrees = (isSample) ? (m-1) : (m);
      double[][] C = X.getArray();
      double c ;
      double s ;
      for (int j = 0; j < n; j++) {
          c = 0;
          s = 0;
          for (int k = 0; k < m; k++) {
            s = s + A[k][j];
          }
          s = s/m;
          for (int k = 0; k < m; k++) {
            c = c + (A[k][j] - s)*(A[k][j] - s);
          }
          C[0][j] = c/degrees;
      }
      return X;
   }

   /** Generate a variance matrix, each row contains values of a pulling.
   @return     An m-by-1 matrix.
   */

   public Matrix varColumns() {
      Matrix X = new Matrix(m,1);
      int degrees = (isSample) ? (n-1) : (n);
      double[][] C = X.getArray();
      double c ;
      double s ;
      for (int i = 0; i < m; i++) {
          c = 0;
          s = 0;
          for (int k = 0; k < n; k++) {
            s = s + A[i][k];
          }
          s = s/n;
          for (int k = 0; k < n; k++) {
            c = c + (A[i][k] - s)*(A[i][k] - s);
          }
          C[i][0] = c/degrees;
      }
      return X;
   }

/////////////////////////////////////////////////////////
//Matrix io methods, in panels, frames or command line.//
/////////////////////////////////////////////////////////

   /** Plot the histogram of the Matrix Columns in a JPanel
   @param slices      Number of slices
   @return      A MatrixHist2D (extends a JPanel)
   */

  public MatrixHist2D toPanelHist2D(int slices) {
    MatrixHist2D mp2d = new MatrixHist2D(this,slices);
    return mp2d;
  }

   /** Plot the histogram of the Matrix Columns in a JFrame
   @param title Title of the JFrame.
   @param slices      Number of slices
   @return      A MatrixHist2D (extends a JPanel)
   */

  public MatrixHist2D toFrameHist2D(String title,int slices) {
    MatrixHist2D mp2d = toPanelHist2D(slices);
    FrameView fv = new FrameView(title,mp2d);
    return mp2d;
  }

   /** Plot the histogram of the Matrix Columns in a JPanel
   @param slicesX      Number of slices in X
   @param slicesY      Number of slices in Y
   @return      A MatrixHist3D (extends a JPanel)
   */

  public MatrixHist3D toPanelHist3D(int slicesX,int slicesY) {
    MatrixHist3D mp3d = new MatrixHist3D(this,slicesX,slicesY);
    return mp3d;
  }

   /** Plot the histogram of the Matrix Columns in a JFrame
   @param title Title of the JFrame.
   @param slicesX      Number of slices in X
   @param slicesY      Number of slices in Y
   @return      A MatrixHist3D (extends a JPanel)
   */

  public MatrixHist3D toFrameHist3D(String title,int slicesX,int slicesY) {
    MatrixHist3D mp3d = toPanelHist3D(slicesX,slicesY);
    FrameView fv = new FrameView(title,mp3d);
    return mp3d;
  }
}