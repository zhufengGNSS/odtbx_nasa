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

package jat.matvec.io.data;

import jat.matvec.data.Matrix;
import jat.matvec.io.data.fileTools.MatrixString;
import jat.matvec.io.data.fileTools.CharFile;

import java.io.File;
//import java.io.IOException;

public class MatrixFile {

  private Matrix M;
  private File file;

  public MatrixFile(File f,Matrix m) {
    M = m;
    file = f;
    CharFile.toFile(file,MatrixString.printMatrix(M));
  }

  public MatrixFile(String fn,Matrix m) {
    M = m;
    file = new File(fn);
    CharFile.toFile(file,MatrixString.printMatrix(M));
  }

  public MatrixFile(File f) {
    file = f;
    if (file.exists()) {
      M = MatrixString.readMatrix(CharFile.fromFile(file));
    } else {
      M = new Matrix(0,0);
      throw new IllegalArgumentException("File does not exist.");
    }
  }

  public MatrixFile(String fn) {
    file = new File(fn);
    if (file.exists()) {
      M = MatrixString.readMatrix(CharFile.fromFile(file));
    } else {
      M = new Matrix(0,0);
      throw new IllegalArgumentException("File does not exist.");
    }
  }

  public Matrix getMatrix() {
    return M;
  }

  public File getFile() {
    return file;
  }

  public String getFileName() {
    return file.getName();
  }
}