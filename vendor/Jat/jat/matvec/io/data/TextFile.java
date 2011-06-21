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

import jat.matvec.data.Text;
import jat.matvec.data.Matrix;
import jat.matvec.io.data.fileTools.CharFile;

import java.io.File;
//import java.io.IOException;

public class TextFile {

  private Text text = new Text("");
  private File file;

  public TextFile(File f) {
    file = f;
    if (file.exists())
      text = new Text(CharFile.fromFile(file));
    else
      text = new Text("");
  }

  public TextFile(File f, Text t) {
    text = t;
    file = f;
    CharFile.toFile(file,text.getString());
  }

  public TextFile(File f, String s) {
    text = new Text(s);
    file = f;
    CharFile.toFile(file,text.getString());
  }

  public TextFile(File f, Matrix X) {
    text = new Text(X);
    file = f;
    CharFile.toFile(file,text.getString());
  }

  public TextFile(String fileName) {
    file = new File(fileName);
    if (file.exists())
      text = new Text(CharFile.fromFile(file));
    else
      text = new Text("");
  }

  public TextFile(String fileName, Text t) {
    text = t;
    file = new File(fileName);
    CharFile.toFile(file,text.getString());
  }

  public TextFile(String fileName, String s) {
    text = new Text(s);
    file = new File(fileName);
    CharFile.toFile(file,text.getString());
  }

  public TextFile(String fileName, Matrix X) {
    text = new Text(X);
    file = new File(fileName);
    CharFile.toFile(file,text.getString());
  }

  public void append(Text t) {
    text = new Text(text.getString() + "\n" + t.getString());
    CharFile.toFile(file,text.getString());
  }

  public void append(String s) {
    text = new Text(text.getString() + "\n" + s);
    CharFile.toFile(file,text.getString());
  }

  public void append(Matrix X) {
    text = new Text(text.getString() + "\n" + new Text(X).getString());
    CharFile.toFile(file,text.getString());
  }

  public Text getText() {
    return text;
  }

  public File getFile() {
    return file;
  }

  public String getFileName() {
    return file.getName();
  }
}