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

package jat.matvec.test;

import jat.matvec.data.*;
import jat.matvec.function.*;
import jat.matvec.io.data.*;
import jat.matvec.io.gui.*;

import java.io.*;

//import javax.swing.*;
//import java.awt.*;


public class Testjmat {

  public Testjmat() {
    testText();
    testMatrix();
    testRandomMatrix();
    testFunction();
  }

  private void title(String s) {
    System.out.println("--------------- "+s+" ---------------");
  }

  private void print(String s) {
    System.out.print(s);
  }
  private void println(String s) {
    System.out.println(s);
  }

  private void printMatrixResult(double x,String s) {
    System.out.println(s + " = " + x);
  }

  private void printMatrixResult(Matrix X,String s) {
    X.toCommandLine(s);
  }


  private void testText() {
    title("test de la classe Text");

    Text t = new Text("texte");

    t.setString("123456...");
    t.setString(Matrix.random(5,5));

    t.merge("abcde...");
    t.merge(t);
    t.merge(Matrix.random(3,3));

    testTextFile();
    testTextWindow();
  }

  private  void testTextFile() {
    Text t1 = new Text("abcd");

    t1.toFile("text.dat");
    t1.toFile(new File("text.dat"));

    Text t2 = Text.fromFile("text.dat");
    t2.toCommandLine();
    TextFile f2 = new TextFile("text2.dat",t2);
    f2.append("suite...");
  }

  private  void testTextWindow() {
    Text t1 = new Text("abcd");
    t1.toFrame("texte");
  }


  private  void testMatrix() {
    title("test de la classe Matrix");

    Matrix A = Matrix.random(5,5);
    Matrix B = Matrix.random(3,3);
    Matrix C = Matrix.random(10,4);
    int[] I = {0,2,1};
    int[] J = {3,4,1};
    int[] I2 = {0,2,1};
    int[] J2 = {3,4,1};

    int[] i = {0,2,4};

    printMatrixResult(A,"A");
    printMatrixResult(B,"B");
    printMatrixResult(C,"C");

    printMatrixResult(A.get(0,0),"A.get(0,0)");
    printMatrixResult(A.get(I,J),"A.get(I,J)");
    printMatrixResult(A.getMatrix(0,2,0,3),"A.getMatrix(0,2,0,3)");
    printMatrixResult(A.getRow(2),"A.getRow(2)");
    printMatrixResult(A.getRows(i),"A.getRows(i)");
    printMatrixResult(A.getColumn(1),"A.getColumn(1)");
    printMatrixResult(C.diag(),"C.diag()");
    printMatrixResult(C.diag(-2),"C.diag(-2)");
    printMatrixResult(C.diag(-7),"C.diag(-7)");
    printMatrixResult(C.diag(2),"C.diag(2)");

    A.set(0,1,2.2);
    printMatrixResult(A,"A.set(0,1,2.2)");
    A.set(I2,J2,5);
    printMatrixResult(A,"A.set(I2,J2,5)");
    A.setMatrix(0,0,B);
    printMatrixResult(A,"A.set(0,0,B)");
    A.setMatrix(1,0,Matrix.random(3,3));
    printMatrixResult(A,"A.setMatrix(1,0,Matrix.random(3,3))");
    A.setMatrix(0,2,0,2,Matrix.random(3,3));
    printMatrixResult(A,"A.setMatrix(0,2,0,2,Matrix.random(3,3))");
    A.setMatrix(0,2,1,1,0.2);
    printMatrixResult(A,"A.setMatrix(0,2,1,1,0.2)");
    A.setRow(1,Matrix.random(1,5));
    printMatrixResult(A,"A.setRow(1,Matrix.random(1,5))");
    A.setRows(i,Matrix.random(3,5));
    printMatrixResult(A,"A.setRows(i,Matrix.random(3,5))");
    A.setColumn(1,Matrix.random(5,1));
    printMatrixResult(A,"A.setColumn(1,Matrix.random(5,1))");
    A.setColumns(i,Matrix.random(5,3));
    printMatrixResult(A,"A.setColumns(i,Matrix.random(5,3))");

    printMatrixResult(C.reshapeRows(20,2),"C.reshapeRows(20,2)");
    printMatrixResult(C.reshapeColumns(20,2),"C.reshapeColumns(20,2)");
    printMatrixResult(C.transpose(),"C.transpose()");
    printMatrixResult(C.merge(Matrix.random(2,4)),"C.merge(Matrix.random(2,4))");
    printMatrixResult(A.norm1(),"A.norm1()");
    printMatrixResult(A.norm2(),"A.norm2()");
    printMatrixResult(A.normInf(),"A.normInf()");
    printMatrixResult(A.normF(),"A.normF()");
    printMatrixResult(A.det(),"A.det()");
    printMatrixResult(A.rank(),"A.rank()");
    printMatrixResult(A.cond(),"A.cond()");
    printMatrixResult(A.trace(),"A.trace()");
    printMatrixResult(A.min(),"A.min()");
    printMatrixResult(A.max(),"A.max()");
    printMatrixResult(A.sum(),"A.sum()");
    printMatrixResult(A.prod(),"A.prod()");

    printMatrixResult(A.plus(A),"A.plus(A)");
    printMatrixResult(A.minus(A),"A.minus(A)");
    printMatrixResult(A.uminus(),"A.uminus()");
    printMatrixResult(A.times(Matrix.random(5,10)),"A.times(Matrix.random(5,10))");
    printMatrixResult(A.times(3.5),"A.times(3.5)");
    printMatrixResult(A.divide(2.5),"A.divide(2.5)");
    printMatrixResult(A.divide(Matrix.random(5,5)),"A.divide(Matrix.random(5,5))");
    printMatrixResult(A.solve(Matrix.random(5,5)),"A.solve(Matrix.random(5,5))");
    printMatrixResult(A.inverse(),"A.inverse()");
    printMatrixResult(A.ebeTimes(Matrix.random(5,5)),"A.ebeTimes(Matrix.random(5,5))");
    printMatrixResult(A.ebeDivide(Matrix.random(5,5)),"A.ebeDivide(Matrix.random(5,5))");
    printMatrixResult(A.ebeCos(),"A.ebeCos()");
    printMatrixResult(A.ebeSin(),"A.ebeSin()");
    printMatrixResult(A.ebeExp(),"A.ebeExp()");
    printMatrixResult(A.ebePow(Matrix.random(5,5)),"A.ebePow(Matrix.random(5,5))");
    printMatrixResult(A.ebeLog(),"A.ebeLog()");
    printMatrixResult(A.ebeInv(),"A.ebeInv()");
    printMatrixResult(A.ebeSqrt(),"A.ebeSqrt()");
    printMatrixResult(A.ebeFun(new DoubleFunctionExpression("3*x+5-x^2","x")),"A.ebeFun(new DoubleFunctionExpression('3*x+5-x^2','x'))");
    printMatrixResult(A.sortedMatrix(0),"A.sortMatrix(0)");
    printMatrixResult(A.findMatrix(0.5),"A.findMatrix(0.5)");
    printMatrixResult(A.findMatrix(">",0.5),"A.findMatrix('!=',0.5)");

    title("test des decompositions matricielles");

    title("test decomposition LU");
    printMatrixResult(A.lu().getL(),"A.lu().getL()");
    printMatrixResult(A.lu().getU(),"A.lu().getU()");
    printMatrixResult(A.lu().getP().times(A),"P*A");
    printMatrixResult(A.lu().getL().times(A.lu().getU()),"A.lu().getL().times(A.lu().getU())");

    title("test decomposition QR");
    printMatrixResult(A.qr().getQ(),"A.qr().getQ()");
    printMatrixResult(A.qr().getR(),"A.qr().getR()");
    printMatrixResult(A,"A");
    printMatrixResult(A.qr().getQ().times(A.qr().getR()),"Q*R");

    title("test decomposition LL'");
    printMatrixResult(A.chol().getL(),"A.chol().getL()");
    printMatrixResult(A,"A");
    printMatrixResult(A.chol().getL().times(A.chol().getL().transpose()),"L*L'");

    title("test decomposition Singular Value");
    printMatrixResult(A.svd().getU(),"A.svd().getU()");
    printMatrixResult(A.svd().getS(),"A.svd().getS()");
    printMatrixResult(A.svd().getV(),"A.svd().getV()");
    printMatrixResult(A,"A");
    printMatrixResult(A.svd().getU().times(A.svd().getS().times(A.svd().getV().transpose())),"U*S*V'");

    title("test decomposition Diagonale");
    printMatrixResult(A.eig().getD(),"A.eig().getD()");
    printMatrixResult(A.eig().getV(),"A.eig().getV()");
    printMatrixResult(A.times(A.eig().getV()),"A*V");
    printMatrixResult(A.eig().getV().times(A.eig().getD()),"V*D");


    testMatrixFile();
    testMatrixPlot();
    testMatrixTable();
  }

  private  void testRandomMatrix() {
    title("test de la classe RandomMatrix");

    RandomMatrix A = RandomMatrix.normal(400,1,0,1);
    A.toFrameHist2D("normal",10);
    A = RandomMatrix.exponential(400,1,1);
    A.toFrameHist2D("exponential",10);
    A = RandomMatrix.logNormal(400,1,0,1);
    A.toFrameHist2D("logNormal",10);
    A = RandomMatrix.triangular(400,1,0,1);
    A.toFrameHist2D("triangular",10);
    A = RandomMatrix.triangular(400,1,0,0.1,1);
    A.toFrameHist2D("triangular",10);
    A = RandomMatrix.uniform(400,1,0,1);
    A.toFrameHist2D("uniform",10);
    A = RandomMatrix.beta(400,1,0.5,0.5);
    A.toFrameHist2D("beta",10);
    A = RandomMatrix.cauchy(400,1,10,10);
    A.toFrameHist2D("cauchy",10);
    A = RandomMatrix.weibull(400,1,0.5,0.5);
    A.toFrameHist2D("weibull",10);
    A = RandomMatrix.rejection(400,1,new DoubleFunctionExpression("sin(x*3.14)","x"),0,1);
    A.toFrameHist2D("rejection",10);


    A = RandomMatrix.uniform(5,3,0,1);
    printMatrixResult(A,"A");
    printMatrixResult(A.mean(),"A.mean()");
    printMatrixResult(A.cov(),"A.cov()");
    printMatrixResult(A.cor(),"A.cor()");
    printMatrixResult(A.var(),"A.vari()");

    printMatrixResult(RandomMatrix.shuffle(A),"RandomMatrix.shuffle(A)");
    printMatrixResult(RandomMatrix.bootstrap(5,5,A),"RandomMatrix.bootstrap(5,5,A)");
    printMatrixResult(RandomMatrix.sampleWithoutReplacement(2,2,A),"RandomMatrix.sampleWithoutReplacement(2,2,A)");
    printMatrixResult(RandomMatrix.sampleWithReplacement(6,6,A),"RandomMatrix.sampleWithReplacement(6,6,A)");

    testMatrixHist();

  }

  private  void testMatrixFile() {
    Matrix M1 = Matrix.random(5,5);

    M1.toFile(new File("essai.mat"));
    M1.toFile("essai.mat");

    Matrix M2 = Matrix.fromFile(new File("essai.mat"));
  }

  private  void testMatrixPlot() {
    Matrix M1 = RandomMatrix.normal(10,2,5,5);
    M1.toFramePlot2D("normal");

    Matrix M2 = RandomMatrix.normal(10,2,0,1);
    Matrix M3 = RandomMatrix.normal(10,2,0,1);
    Matrix[] M = {M2,M3};
    FrameView fv = new FrameView("normal", new MatrixPlot2D(M));

    Matrix MM1 = RandomMatrix.normal(10,3,5,5);
    MM1.toFramePlot3D("normal");

    Matrix MM2 = RandomMatrix.normal(10,3,0,1);
    Matrix MM3 = RandomMatrix.normal(10,3,0,1);
    Matrix[] MM = {MM2,MM3};
    FrameView fvv = new FrameView("normal", new MatrixPlot3D(MM));
  }

  private  void testMatrixTable() {
    Matrix M1 = Matrix.random(50,15);

    M1.toFrameTable("random");
  }

  private  void testMatrixHist() {
    RandomMatrix rand1 = RandomMatrix.normal(40,1,0,1);
    rand1.toFrameHist2D("normal",5);

    RandomMatrix rand2 = RandomMatrix.normal(40,1,0,1);
    RandomMatrix rand3 = RandomMatrix.normal(40,1,0,1);
    RandomMatrix[] rand = {rand2,rand3};
    FrameView fv = new FrameView("normal", new MatrixHist2D(rand,5));

    RandomMatrix randd1 = RandomMatrix.normal(100,2,0,1);
    randd1.toFrameHist3D("normal",5,5);

    RandomMatrix randd2 = RandomMatrix.normal(100,2,0,1);
    RandomMatrix randd3 = RandomMatrix.normal(100,2,0,1);
    RandomMatrix[] randd = {randd2,randd3};
    FrameView fvv = new FrameView("normal", new MatrixHist3D(randd,5,5));

  }

  private void testFunction() {
    title("test de la classe Function");

    testDoubleFunctionExpression();
    testMatrixFunctionExpression();
    testTestDoubleFunctionExpression();
    testInvokeDoubleFunction();
    testInvokeMatrixFunction();
    testDoubleFunctionInterpolation();

    testFunctionPlot();
  }

  private  void testFunctionPlot() {
    DoubleFunction[] df1 = {new DoubleFunctionExpression("cos(x)","x"),new DoubleFunctionExpression("sin(x)","x")};
    FrameView fv1 = new FrameView("cos(x)", new FunctionPlot2D(df1,-6,6));

    String[] s = {"x","y"};
    DoubleFunction df2 = new DoubleFunctionExpression("3*x/y",s);
    FrameView fv2 = new FrameView("3*x/y", new FunctionPlot3D(df2,0.1,1,0.1,1));
  }

  private  void testDoubleFunctionExpression() {
    title("test de la fonction de double");


    String[] vars = new String[2];
    vars[0] = "x";
    vars[1] = "y";
    String expression = new String("3*x+1/y");
    DoubleFunctionExpression f = new DoubleFunctionExpression(expression,vars);

    Matrix vals = Matrix.random(2,1);
    vals.toCommandLine("[x,y]");
    double result1 = f.eval(vals);
    println("resultat :");
    println(expression+" = "+result1);
  }

  private  void testMatrixFunctionExpression() {
    title("test de la fonction matricielle");

    String[] vars = new String[3];
    vars[0] = "x";
    vars[1] = "y";
    vars[2] = "z";
    String expression = new String("x*y+diag(z)");
    MatrixFunctionExpression f = new MatrixFunctionExpression(expression,vars);

    Matrix[] vals = new Matrix[3];
    vals[0] = Matrix.random(3,1);
    vals[1] = Matrix.random(1,1);
    vals[2] = Matrix.random(3,3);
    Matrix result1 = f.eval(vals);
    vals[0].toCommandLine("x");
    vals[1].toCommandLine("y");
    vals[2].toCommandLine("z");
    result1.toCommandLine("resultat : " + expression);
  }

  private  void testTestDoubleFunctionExpression() {
    title("test de la fonction bool�enne");

    String[] vars = new String[3];
    vars[0] = "x";
    vars[1] = "y";
    vars[2] = "z";
    String expression = new String("3*x+1/y<z");
    TestDoubleFunctionExpression f = new TestDoubleFunctionExpression(expression,vars);

    Matrix vals = Matrix.random(3,1);
    vals.toCommandLine("[x,y,z]");
    boolean result1 = f.eval(vals);
    println("resultat :");
    println(expression+" = "+result1);

    vals.set(2,0,10);
    vals.set(1,0,1);

    vals.toCommandLine("[x,y,z]");
    boolean result2 = f.eval(vals);
    println("resultat :");
    println(expression+" = "+result2);

  }

  private  void testInvokeDoubleFunction() {
    title("test de la fonction d'invocation - double");
    InvokeDoubleFunction i = new InvokeDoubleFunction("essai.bat","essai.mat");
    println("resultat :" + i.eval());
  }

  private  void testInvokeMatrixFunction() {
    title("test de la fonction d'invocation - matrix");
    InvokeMatrixFunction i = new InvokeMatrixFunction("essai.bat","essai.mat");
    i.eval().toCommandLine("resultat");
  }

  private  void testDoubleFunctionInterpolation() {
    title("test de la fonction d'interpolation");

    double[][] X = {{0,0},{2,0},{0,2},{2,2}};
    double[] Y = {1,2,3,4};
    DoubleFunctionInterpolation dfi = new DoubleFunctionInterpolation(X,Y);
    double[] x = {1.5,0.5};
    double d = dfi.eval(x);

    println("points � interpoler :");
    for (int i = 0; i < X.length; i++) {
      print("(");
      for (int j = 0; j < X[i].length; j++) {
        print(X[i][j]  + ",");
      }
      println(") -> " + Y[i]);
    }
    println("resultat :");
    print("(");
    for (int j = 0; j < x.length; j++) {
      print(x[j]  + ",");
    }
    print(") -> " + d);

  }

  public static void main(String[] args) {
    Testjmat t = new Testjmat();
  }

}
