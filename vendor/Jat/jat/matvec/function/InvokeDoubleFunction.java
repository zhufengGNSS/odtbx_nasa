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
import jat.matvec.io.data.MatrixFile;
import java.io.File;
//import java.io.IOException;

public class InvokeDoubleFunction {

  File functionFile;
  File resultFile;

  public InvokeDoubleFunction(String fn,String rf) {
      functionFile = new File(fn);
      resultFile = new File(rf);
  }

  public InvokeDoubleFunction(File fn,File rf) {
      functionFile = fn;
      resultFile = rf;
  }

  public double eval() {
    try {
      Process p = Runtime.getRuntime().exec(functionFile.getName());
      p.waitFor();
      MatrixFile mf = new MatrixFile(resultFile);
      Matrix X = mf.getMatrix();
      return new Double(X.get(0,0)).doubleValue();
    } catch (Exception e) {
      System.out.println("Error : File " + resultFile +" unreadable : "+e);
      return Double.NaN;
    }
  }

}