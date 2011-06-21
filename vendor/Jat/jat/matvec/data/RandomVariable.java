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

import jat.matvec.function.DoubleFunction;

/**
<P>
   The RandomVaraibale Class provides static methods for generating random numbers.

@author Yann RICHET.
@version 2.0
*/


public class RandomVariable {

   /** Generate a random number between 0 and 1.
   @return      A double between 0 and 1.
   */

  protected static double rand() {
    double x = Math.random();
    return x;
  }

   /** Generate a random integer.
   @param i0    Min of the random variable.
   @param i1    Max of the random variable.
   @return      An int between i0 and i1.
   */

  protected static int randInt(int i0, int i1) {
    double x = rand();
    int i = i0 + new Double(Math.floor((i1-i0+1)*x)).intValue();
    return i;
  }

  /** Generate a random number from a uniform random variable.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A double.
   */

  public static double uniform(double min,double max) {
    double x = min+(max-min)*rand();
    return x;
  }

  /** Generate a random number from a discrete random variable.
   @param values    Discrete values.
   @param prob    Probability of each value.
   @return      A double.
   */

  public static double dirac(double[] values,double[] prob) {
    double[] prob_cumul = new double[values.length];
    prob_cumul[0] = prob[0];
    for (int i = 1; i < values.length; i++) {
      prob_cumul[i] = prob_cumul[i-1] + prob[i];
    }
    double y = rand();
    double x = 0;
    for (int i = 0; i < values.length-1; i++) {
      if ((y > prob_cumul[i])&(y < prob_cumul[i+1])) {
        x = values[i];
      }
    }
    return x;
  }

   /** Generate a random number from a Gaussian (Normal) random variable.
   @param mu    Mean of the random variable.
   @param sigma    Standard deviation of the random variable.
   @return      A double.
   */

  public static double normal(double mu,double sigma) {
    double x = mu + sigma*Math.cos(2*Math.PI*rand())*Math.sqrt(-2*Math.log(rand()));
    return x;
  }

   /** Generate a random number from a LogNormal random variable.
   @param mu    Mean of the Normal random variable.
   @param sigma    Standard deviation of the Normal random variable.
   @return      A double.
   */

  public static double logNormal(double mu,double sigma) {
    double x = mu + sigma*Math.cos(2*Math.PI*rand())*Math.sqrt(-2*Math.log(rand()));
    return x;
  }

   /** Generate a random number from an exponantial random variable (Mean = 1/lambda, variance = 1/lambda^2).
   @param lambda    Parmaeter of the exponential random variable.
   @return      A double.
   */

  public static double exponential(double lambda) {
    double x = -1/lambda*Math.log(rand());
    return x;
  }


   /** Generate a random number from a symetric triangular random variable.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A double.
   */

  public static double triangular(double min,double max) {
    double x = min/2+(max-min)*rand()/2 + min/2+(max-min)*rand()/2;
    return x;
  }

     /** Generate a random number from a non-symetric triangular random variable.
   @param min    Min of the random variable.
   @param med    Value of the random variable with max density.
   @param max    Max of the random variable.
   @return      A double.
   */

  public static double triangular(double min,double med,double max) {
    double y = rand();
    //if min < x < med, y = (x-min)?/(max-min)(med-min), else, med < x < max, and y = 1-(max-x)?/(max-min)(max-med)
    double x = (y < ((med-min)/(max-min))) ? (min+Math.sqrt(y*(max-min)*(med-min))) : (max-Math.sqrt((1-y)*(max-min)*(max-med)));
    return x;
  }

     /** Generate a random number from a beta random variable.
   @param a    First parameter of the Beta random variable.
   @param b    Second parameter of the Beta random variable.
   @return      A double.
   */

  public static double beta(double a,double b) {
    double try_x;
    double try_y;
    do {
      try_x = Math.pow(rand(),1/a);
      try_y = Math.pow(rand(),1/b);
    } while ((try_x + try_y) > 1);
    return try_x/(try_x + try_y);
  }

     /** Generate a random number from a Cauchy random variable (Mean = Inf, and Variance = Inf).
   @param mu    Median of the Weibull random variable
   @param sigma    Second parameter of the Cauchy random variable.
   @return      A double.
   */

  public static double cauchy(double mu,double sigma) {
    double x = sigma*Math.tan(Math.PI*(rand()-0.5)) + mu;
    return x;
  }

     /** Generate a random number from a Weibull random variable.
   @param lambda    First parameter of the Weibull random variable.
   @param c    Second parameter of the Weibull random variable.
   @return      A double.
   */

  public static double weibull(double lambda,double c) {
    double x = Math.pow(-Math.log(1-rand()),1/c)/lambda;
    return x;
  }

  /** Generate a random number from a random variable definied by its density function, using the rejection technic.
   *  !!! WARNING : this simulation technic can take a very long time !!!
   @param fun    Density function (may be not normalized) of the random variable.
   @param maxFun    Max of the function.
   @param min    Min of the random variable.
   @param max    Max of the random variable.
   @return      A double.
   */

  public static double rejection(DoubleFunction fun, double maxFun, double min, double max) {
    double[] try_x = new double[1];
    double try_y;
    do {
      try_x[0] = min + rand()*(max-min);
      try_y = rand()*maxFun;
    } while (fun.eval(try_x) < try_y);
    return try_x[0];
  }

}