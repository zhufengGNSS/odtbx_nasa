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

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;

/**
 * <p>Titre : JAva MAtrix TOols</p>
 * <p>Description : </p>
 * <p>Copyright : Copyright (c) 2002</p>
 * <p>Socit : IRSN</p>
 * @author Yann RICHET
 * @version 1.0
 */

public class CharFile {

  public static String fromFile(File file)  {
    String string = new String("");
    try {
      FileReader fr = new FileReader(file);
      BufferedReader b = new BufferedReader(fr);
      boolean eof = false;
      while (!eof) {
        String line = b.readLine();
        if (line == null) {
          eof = true;
          string = string.substring(0,string.length()-1);
        } else
          string = string + line + "\n";
      }
      b.close();
    } catch (IOException e) {
      System.out.println("File " + file.getName() + " is unreadable.");
    }
    return string;
  }

  public static void toFile(File file,String s) {
    try {
      FileWriter fw = new FileWriter(file);
      BufferedWriter bw = new BufferedWriter(fw);
      bw.write(s);
      bw.close();
    } catch (IOException e) {
      System.out.println("File " + file.getName() + " is unwritable.");
    }
  }

}