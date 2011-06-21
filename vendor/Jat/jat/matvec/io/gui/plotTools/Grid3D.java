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

//import javax.swing.*;
import java.awt.*;

public class Grid3D {

  private Axe3D axe;

  private int numGrid = 5;

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

  private double[] Pl_0 = {0,0,0};
  private int[] Sc_0 = new int[2];

  private double[][][] Pl_Xgrid ;
  private double[][][] Pl_Ygrid ;
  private double[][][] Pl_Zgrid ;

  private int[][][] Sc_Xgrid ;
  private int[][][] Sc_Ygrid ;
  private int[][][] Sc_Zgrid ;

  private double Xpow;
  private double Ypow;
  private double Zpow;

  public Grid3D(Axe3D ax) {
    axe = ax;

    Pl_XMin = axe.getPlXMin();
    Pl_XMax = axe.getPlXMax();
    Pl_YMin = axe.getPlYMin();
    Pl_YMax = axe.getPlYMax();
    Pl_ZMin = axe.getPlZMin();
    Pl_ZMax = axe.getPlZMax();
    Pl_0 = axe.getPl0();
    Sc_XMin = axe.getScXMin();
    Sc_XMax = axe.getScXMax();
    Sc_YMin = axe.getScYMin();
    Sc_YMax = axe.getScYMax();
    Sc_ZMin = axe.getScZMin();
    Sc_ZMax = axe.getScZMax();
    Sc_0 = axe.getSc0();

    setGridStep();

    PlScConvert();
  }

  private void PlScConvert() {

    Sc_Xgrid = axe.Pl2Sc(Pl_Xgrid);
    Sc_Ygrid = axe.Pl2Sc(Pl_Ygrid);
    Sc_Zgrid = axe.Pl2Sc(Pl_Zgrid);
  }



  private void setGridStep() {

    Xpow = Math.floor(Math.log((Pl_XMax[0] - Pl_XMin[0])/numGrid)/Math.log(10));
    double Xpas = approx((Pl_XMax[0] - Pl_XMin[0])/numGrid,Xpow);
    int Xnb_pos = (int)(Math.floor((Pl_XMax[0] - Pl_0[0])/Xpas));
    int Xnb_neg = (int)(Math.floor((Pl_0[0] - Pl_XMin[0])/Xpas));

    Pl_Xgrid = new double[Xnb_pos+Xnb_neg+2][3][5];
    Sc_Xgrid = new int[Xnb_pos+Xnb_neg+2][2][5];

    Ypow = Math.floor(Math.log((Pl_YMax[1] - Pl_YMin[1])/numGrid)/Math.log(10));
    double Ypas = approx((Pl_YMax[1] - Pl_YMin[1])/numGrid,Ypow);
    int Ynb_pos = (int)(Math.floor((Pl_YMax[1] - Pl_0[1])/Ypas));
    int Ynb_neg = (int)(Math.floor((Pl_0[1] - Pl_YMin[1])/Ypas));

    Pl_Ygrid = new double[Ynb_pos+Ynb_neg+2][3][5];
    Sc_Ygrid = new int[Ynb_pos+Ynb_neg+2][2][5];

    Zpow = Math.floor(Math.log((Pl_ZMax[2] - Pl_ZMin[2])/numGrid)/Math.log(10));
    double Zpas = approx((Pl_ZMax[2] - Pl_ZMin[2])/numGrid,Zpow);
    int Znb_pos = (int)(Math.floor((Pl_ZMax[2] - Pl_0[2])/Zpas));
    int Znb_neg = (int)(Math.floor((Pl_0[2] - Pl_ZMin[2])/Zpas));

    Pl_Zgrid = new double[Znb_pos+Znb_neg+2][3][5];
    Sc_Zgrid = new int[Znb_pos+Znb_neg+2][2][5];


    Pl_Xgrid[0][0][0] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][0] = Pl_YMin[1];
    Pl_Xgrid[0][2][0] = Pl_YMin[2];
    Pl_Xgrid[0][0][1] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][1] = Pl_YMax[1];
    Pl_Xgrid[0][2][1] = Pl_YMax[2];
    Pl_Xgrid[0][0][2] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][2] = Pl_0[1];
    Pl_Xgrid[0][2][2] = Pl_0[2];
    Pl_Xgrid[0][0][3] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][3] = Pl_ZMin[1];
    Pl_Xgrid[0][2][3] = Pl_ZMin[2];
    Pl_Xgrid[0][0][4] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][4] = Pl_ZMax[1];
    Pl_Xgrid[0][2][4] = Pl_ZMax[2];

    Pl_Ygrid[0][0][0] = Pl_XMin[0];
    Pl_Ygrid[0][1][0] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][2][0] = Pl_XMin[2];
    Pl_Ygrid[0][0][1] = Pl_XMax[0];
    Pl_Ygrid[0][1][1] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][2][1] = Pl_XMax[2];
    Pl_Ygrid[0][0][2] = Pl_0[0];
    Pl_Ygrid[0][1][2] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][2][2] = Pl_0[2];
    Pl_Ygrid[0][0][3] = Pl_ZMin[0];
    Pl_Ygrid[0][1][3] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][2][3] = Pl_ZMin[2];
    Pl_Ygrid[0][0][4] = Pl_ZMax[0];
    Pl_Ygrid[0][1][4] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][2][4] = Pl_ZMax[2];

    Pl_Zgrid[0][0][0] = Pl_XMin[0];
    Pl_Zgrid[0][1][0] = Pl_XMin[1];
    Pl_Zgrid[0][2][0] = approx(Pl_0[2] - Znb_neg*Zpas,Zpow);
    Pl_Zgrid[0][0][1] = Pl_XMax[0];
    Pl_Zgrid[0][1][1] = Pl_XMax[1];
    Pl_Zgrid[0][2][1] = approx(Pl_0[2] - Znb_neg*Zpas,Zpow);
    Pl_Zgrid[0][0][2] = Pl_0[0];
    Pl_Zgrid[0][1][2] = Pl_0[1];
    Pl_Zgrid[0][2][2] = approx(Pl_0[2] - Znb_neg*Zpas,Zpow);
    Pl_Zgrid[0][0][3] = Pl_YMin[0];
    Pl_Zgrid[0][1][3] = Pl_YMin[1];
    Pl_Zgrid[0][2][3] = approx(Pl_0[2] - Znb_neg*Zpas,Zpow);
    Pl_Zgrid[0][0][4] = Pl_YMax[0];
    Pl_Zgrid[0][1][4] = Pl_YMax[1];
    Pl_Zgrid[0][2][4] = approx(Pl_0[2] - Znb_neg*Zpas,Zpow);


    for (int i =1;i<Pl_Xgrid.length-1;i++) {
      Pl_Xgrid[i][0][0] = Pl_Xgrid[i-1][0][0] + Xpas;
      Pl_Xgrid[i][1][0] = Pl_Xgrid[i-1][1][0];
      Pl_Xgrid[i][2][0] = Pl_Xgrid[i-1][2][0];
      Pl_Xgrid[i][0][1] = Pl_Xgrid[i-1][0][1] + Xpas;
      Pl_Xgrid[i][1][1] = Pl_Xgrid[i-1][1][1];
      Pl_Xgrid[i][2][1] = Pl_Xgrid[i-1][2][1];
      Pl_Xgrid[i][0][2] = Pl_Xgrid[i-1][0][2] + Xpas;
      Pl_Xgrid[i][1][2] = Pl_Xgrid[i-1][1][2];
      Pl_Xgrid[i][2][2] = Pl_Xgrid[i-1][2][2];
      Pl_Xgrid[i][0][3] = Pl_Xgrid[i-1][0][3] + Xpas;
      Pl_Xgrid[i][1][3] = Pl_Xgrid[i-1][1][3];
      Pl_Xgrid[i][2][3] = Pl_Xgrid[i-1][2][3];
      Pl_Xgrid[i][0][4] = Pl_Xgrid[i-1][0][4] + Xpas;
      Pl_Xgrid[i][1][4] = Pl_Xgrid[i-1][1][4];
      Pl_Xgrid[i][2][4] = Pl_Xgrid[i-1][2][4];
    }
    for (int j =1;j<Pl_Ygrid.length-1;j++) {
      Pl_Ygrid[j][0][0] = Pl_Ygrid[j-1][0][0];
      Pl_Ygrid[j][1][0] = Pl_Ygrid[j-1][1][0] + Ypas;
      Pl_Ygrid[j][2][0] = Pl_Ygrid[j-1][2][0];
      Pl_Ygrid[j][0][1] = Pl_Ygrid[j-1][0][1];
      Pl_Ygrid[j][1][1] = Pl_Ygrid[j-1][1][1] + Ypas;
      Pl_Ygrid[j][2][1] = Pl_Ygrid[j-1][2][1];
      Pl_Ygrid[j][0][2] = Pl_Ygrid[j-1][0][2];
      Pl_Ygrid[j][1][2] = Pl_Ygrid[j-1][1][2] + Ypas;
      Pl_Ygrid[j][2][2] = Pl_Ygrid[j-1][2][2];
      Pl_Ygrid[j][0][3] = Pl_Ygrid[j-1][0][3];
      Pl_Ygrid[j][1][3] = Pl_Ygrid[j-1][1][3] + Ypas;
      Pl_Ygrid[j][2][3] = Pl_Ygrid[j-1][2][3];
      Pl_Ygrid[j][0][4] = Pl_Ygrid[j-1][0][4];
      Pl_Ygrid[j][1][4] = Pl_Ygrid[j-1][1][4] + Ypas;
      Pl_Ygrid[j][2][4] = Pl_Ygrid[j-1][2][4];
    }
    for (int k =1;k<Pl_Zgrid.length-1;k++) {
      Pl_Zgrid[k][0][0] = Pl_Zgrid[k-1][0][0];
      Pl_Zgrid[k][1][0] = Pl_Zgrid[k-1][1][0];
      Pl_Zgrid[k][2][0] = Pl_Zgrid[k-1][2][0] + Zpas;
      Pl_Zgrid[k][0][1] = Pl_Zgrid[k-1][0][1];
      Pl_Zgrid[k][1][1] = Pl_Zgrid[k-1][1][1];
      Pl_Zgrid[k][2][1] = Pl_Zgrid[k-1][2][1] + Zpas;
      Pl_Zgrid[k][0][2] = Pl_Zgrid[k-1][0][2];
      Pl_Zgrid[k][1][2] = Pl_Zgrid[k-1][1][2];
      Pl_Zgrid[k][2][2] = Pl_Zgrid[k-1][2][2] + Zpas;
      Pl_Zgrid[k][0][3] = Pl_Zgrid[k-1][0][3];
      Pl_Zgrid[k][1][3] = Pl_Zgrid[k-1][1][3];
      Pl_Zgrid[k][2][3] = Pl_Zgrid[k-1][2][3] + Zpas;
      Pl_Zgrid[k][0][4] = Pl_Zgrid[k-1][0][4];
      Pl_Zgrid[k][1][4] = Pl_Zgrid[k-1][1][4];
      Pl_Zgrid[k][2][4] = Pl_Zgrid[k-1][2][4] + Zpas;
    }
  }

  private String troncatedString(double d, double power,int prec) {
    return new String("" + approx(d,(double)(power-prec)));
  }

  public String troncatedStringX(double d,int prec) {
    return troncatedString(d,Xpow,prec);
  }
  public String troncatedStringY(double d,int prec) {
    return troncatedString(d,Ypow,prec);
  }
  public String troncatedStringZ(double d,int prec) {
    return troncatedString(d,Zpow,prec);
  }

  private double approx(double X, double power) {
    int i = (int)(Math.rint(X*Math.pow(10,-Math.rint(power))));
    double Y = i/Math.pow(10,Math.rint(-power));
    return Y;
  }

  public void draw(Graphics2D comp2D) {
    PlScConvert();
    comp2D.setColor(Color.darkGray);
    comp2D.setFont(new Font("Arial",Font.PLAIN,12));
    comp2D.setColor(Color.lightGray);
    for (int i=0;i<Pl_Xgrid.length;i++) {
        comp2D.drawLine(Sc_Xgrid[i][0][0],Sc_Xgrid[i][1][0],Sc_Xgrid[i][0][1],Sc_Xgrid[i][1][1]);
        comp2D.drawLine(Sc_Xgrid[i][0][3],Sc_Xgrid[i][1][3],Sc_Xgrid[i][0][4],Sc_Xgrid[i][1][4]);
        comp2D.drawString(troncatedStringX(Pl_Xgrid[i][0][2],0),Sc_Xgrid[i][0][2],Sc_Xgrid[i][1][2]);
    }
    for (int j=0;j<Pl_Ygrid.length;j++) {
        comp2D.drawLine(Sc_Ygrid[j][0][0],Sc_Ygrid[j][1][0],Sc_Ygrid[j][0][1],Sc_Ygrid[j][1][1]);
        comp2D.drawLine(Sc_Ygrid[j][0][3],Sc_Ygrid[j][1][3],Sc_Ygrid[j][0][4],Sc_Ygrid[j][1][4]);
        comp2D.drawString(troncatedStringY(Pl_Ygrid[j][1][2],0),Sc_Ygrid[j][0][2],Sc_Ygrid[j][1][2]);
    }
    for (int k=0;k<Pl_Zgrid.length;k++) {
        comp2D.drawLine(Sc_Zgrid[k][0][0],Sc_Zgrid[k][1][0],Sc_Zgrid[k][0][1],Sc_Zgrid[k][1][1]);
        comp2D.drawLine(Sc_Zgrid[k][0][3],Sc_Zgrid[k][1][3],Sc_Zgrid[k][0][4],Sc_Zgrid[k][1][4]);
        comp2D.drawString(troncatedStringZ(Pl_Zgrid[k][2][2],0),Sc_Zgrid[k][0][2],Sc_Zgrid[k][1][2]);
    }
  }

}
