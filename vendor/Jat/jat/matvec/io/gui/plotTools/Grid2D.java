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

public class Grid2D {

  private Axe2D axe;

  private int numGrid = 5;

  private double[] Pl_XMin = new double[2];
  private double[] Pl_XMax = new double[2];
  private double[] Pl_YMin = new double[2];
  private double[] Pl_YMax = new double[2];

  private int[] Sc_XMin = new int[2];
  private int[] Sc_XMax = new int[2];
  private int[] Sc_YMin = new int[2];
  private int[] Sc_YMax = new int[2];

  private double[] Pl_0 = {0,0};
  private int[] Sc_0 = new int[2];

  private double[][][] Pl_Xgrid ;
  private double[][][] Pl_Ygrid ;

  private int[][][] Sc_Xgrid ;
  private int[][][] Sc_Ygrid ;

  private int Xpow;
  private int Ypow;

  public Grid2D(Axe2D ax) {
    axe = ax;

    Pl_XMin = axe.getPlXMin();
    Pl_XMax = axe.getPlXMax();
    Pl_YMin = axe.getPlYMin();
    Pl_YMax = axe.getPlYMax();
    Pl_0 = axe.getPl0();
    Sc_XMin = axe.getScXMin();
    Sc_XMax = axe.getScXMax();
    Sc_YMin = axe.getScYMin();
    Sc_YMax = axe.getScYMax();
    Sc_0 = axe.getSc0();

    setGridStep();

    PlScConvert();
  }

  private void PlScConvert() {

    Sc_Xgrid = axe.Pl2Sc(Pl_Xgrid);
    Sc_Ygrid = axe.Pl2Sc(Pl_Ygrid);
  }



  private void setGridStep() {

    Xpow = (int)Math.floor(Math.log((Pl_XMax[0] - Pl_XMin[0])/numGrid)/Math.log(10));
    double Xpas = approx((Pl_XMax[0] - Pl_XMin[0])/numGrid,Xpow);
    int Xnb_pos = (int)(Math.floor((Pl_XMax[0] - Pl_0[0])/Xpas));
    int Xnb_neg = (int)(Math.floor((Pl_0[0] - Pl_XMin[0])/Xpas));

    Pl_Xgrid = new double[Xnb_pos+Xnb_neg+2][2][3];
    Sc_Xgrid = new int[Xnb_pos+Xnb_neg+2][2][3];


    Ypow = (int)Math.floor(Math.log((Pl_YMax[1] - Pl_YMin[1])/numGrid)/Math.log(10));
    double Ypas = approx((Pl_YMax[1] - Pl_YMin[1])/numGrid,Ypow);
    int Ynb_pos = (int)(Math.floor((Pl_YMax[1] - Pl_0[1])/Ypas));
    int Ynb_neg = (int)(Math.floor((Pl_0[1] - Pl_YMin[1])/Ypas));

    Pl_Ygrid = new double[Ynb_pos+Ynb_neg+2][2][3];
    Sc_Ygrid = new int[Ynb_pos+Ynb_neg+2][2][3];


    Pl_Xgrid[0][0][0] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][0] = Pl_YMin[1];
    Pl_Xgrid[0][0][1] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][1] = Pl_YMax[1];
    Pl_Xgrid[0][0][2] = approx(Pl_0[0] - Xnb_neg*Xpas,Xpow);
    Pl_Xgrid[0][1][2] = Pl_0[1];

    Pl_Ygrid[0][0][0] = Pl_XMin[0];
    Pl_Ygrid[0][1][0] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][0][1] = Pl_XMax[0];
    Pl_Ygrid[0][1][1] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);
    Pl_Ygrid[0][0][2] = Pl_0[0];
    Pl_Ygrid[0][1][2] = approx(Pl_0[1] - Ynb_neg*Ypas,Ypow);

    for (int i =1;i<Pl_Xgrid.length-1;i++) {
      Pl_Xgrid[i][0][0] = Pl_Xgrid[i-1][0][0] + Xpas;
      Pl_Xgrid[i][1][0] = Pl_Xgrid[i-1][1][0];
      Pl_Xgrid[i][0][1] = Pl_Xgrid[i-1][0][1] + Xpas;
      Pl_Xgrid[i][1][1] = Pl_Xgrid[i-1][1][1];
      Pl_Xgrid[i][0][2] = Pl_Xgrid[i-1][0][2] + Xpas;
      Pl_Xgrid[i][1][2] = Pl_Xgrid[i-1][1][2];
    }
    for (int j =1;j<Pl_Ygrid.length-1;j++) {
      Pl_Ygrid[j][0][0] = Pl_Ygrid[j-1][0][0];
      Pl_Ygrid[j][1][0] = Pl_Ygrid[j-1][1][0] + Ypas;
      Pl_Ygrid[j][0][1] = Pl_Ygrid[j-1][0][1];
      Pl_Ygrid[j][1][1] = Pl_Ygrid[j-1][1][1] + Ypas;
      Pl_Ygrid[j][0][2] = Pl_Ygrid[j-1][0][2];
      Pl_Ygrid[j][1][2] = Pl_Ygrid[j-1][1][2] + Ypas;
    }
  }

  private String troncatedString(double d, int power, int prec) {
    return new String("" + approx(d,(double)(power-prec)));
  }

  public String troncatedStringX(double d, int prec) {
    return troncatedString(d,Xpow,prec);
  }
  public String troncatedStringY(double d, int prec) {
    return troncatedString(d,Ypow,prec);
  }

  private double approx(double X, double power) {
    int i = (int)(Math.rint(X*Math.pow(10,-Math.rint(power))));
    double Y = i/Math.pow(10,Math.rint(-power));
    return Y;
  }

  public void draw(Graphics2D comp2D) {
    //PlScConvert();
    comp2D.setColor(Color.darkGray);
    comp2D.setFont(new Font("Arial",Font.PLAIN,12));
    comp2D.setColor(Color.lightGray);
    for (int i=0;i<Pl_Xgrid.length;i++) {
        comp2D.drawLine(Sc_Xgrid[i][0][0],Sc_Xgrid[i][1][0],Sc_Xgrid[i][0][1],Sc_Xgrid[i][1][1]);
        comp2D.drawString(troncatedStringX(Pl_Xgrid[i][0][2],0),Sc_Xgrid[i][0][2],Sc_Xgrid[i][1][2]);
    }
    for (int j=0;j<Pl_Ygrid.length;j++) {
        comp2D.drawLine(Sc_Ygrid[j][0][0],Sc_Ygrid[j][1][0],Sc_Ygrid[j][0][1],Sc_Ygrid[j][1][1]);
        comp2D.drawString(troncatedStringY(Pl_Ygrid[j][1][2],0),Sc_Ygrid[j][0][2],Sc_Ygrid[j][1][2]);
    }
  }

}
