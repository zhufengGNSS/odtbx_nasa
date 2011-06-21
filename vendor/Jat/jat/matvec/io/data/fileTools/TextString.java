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

package jat.matvec.io.data.fileTools;

import jat.matvec.data.Text;

public class TextString {

  private static int decimalSize = 10;

  private Text T;
  private String S;

  public TextString(Text t) {
    T = t;
    S = TextString.printText(T);
  }

  public TextString(String s) {
    S = s;
    T = readText(S);
  }

  public Text getText() {
    return T;
  }

  public String getString() {
    return S;
  }

  public static Text readText(String s) {
    return new Text(s);
  }

  public static String printText(Text t) {
    return t.getString();
  }

}