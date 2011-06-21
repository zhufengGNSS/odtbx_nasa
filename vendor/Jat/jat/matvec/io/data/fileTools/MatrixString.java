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

import jat.matvec.data.Matrix;

//import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
//import java.text.FieldPosition;

import java.io.StringReader;
//import java.io.PrintWriter;
//import java.io.BufferedReader;
import java.io.StreamTokenizer;

import java.io.IOException;

public class MatrixString {

  private static int decimalSize = 10;

  private Matrix M;
  private String S;

  public MatrixString(Matrix m) {
    M = m;
    S = printMatrix(M);
  }

  public MatrixString(String s) {
    S = s;
    M = readMatrix(S);
  }

  public Matrix getMatrix() {
    return M;
  }

  public String getString() {
    return S;
  }

  public static String printMatrix(Matrix m) {

      String str = new String("");

      DecimalFormat format = new DecimalFormat();
      format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
      format.setMinimumIntegerDigits(1);
      format.setMaximumFractionDigits(decimalSize);
      format.setMinimumFractionDigits(decimalSize);
      format.setGroupingUsed(false);


      for (int i = 0; i < m.getRowDimension(); i++) {
         for (int j = 0; j < m.getColumnDimension(); j++) {
            String s = format.format(m.get(i,j)); // format the number
            int padding = Math.max(1,decimalSize+2-s.length()); // At _least_ 1 space
            for (int k = 0; k < padding; k++)
               str = str + " ";
            str = str + s;
         }
         str = str + "\n";
      }
      str = str + "\n";   // end with blank line.

    return str;
  }



   /** Read a matrix from a stream.  The format is the same the print method,
     * so printed matrices can be read back in (provided they were printed using
     * US Locale).  Elements are separated by
     * whitespace, all the elements for each row appear on a single line,
     * the last row is followed by a blank line.
   @param s  the input String.
   @return A Matrix.
   */

   public static Matrix readMatrix (String s) {
      try {
        StreamTokenizer tokenizer= new StreamTokenizer(new StringReader(s));

      // Although StreamTokenizer will parse numbers, it doesn't recognize
      // scientific notation (E or D); however, Double.valueOf does.
      // The strategy here is to disable StreamTokenizer's number parsing.
      // We'll only get whitespace delimited words, EOL's and EOF's.
      // These words should all be numbers, for Double.valueOf to parse.

        tokenizer.resetSyntax();
        tokenizer.wordChars(0,255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);
        java.util.Vector v = new java.util.Vector();

        // Ignore initial empty lines
        while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
        if (tokenizer.ttype == StreamTokenizer.TT_EOF)
	  throw new IOException("Unexpected EOF on matrix read.");
        do {
          v.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
        } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

        int n = v.size();  // Now we've got the number of columns!
        double row[] = new double[n];
        for (int j=0; j<n; j++)  // extract the elements of the 1st row.
          row[j]=((Double)v.elementAt(j)).doubleValue();
        v.removeAllElements();
        v.addElement(row);  // Start storing rows instead of columns.
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
          // While non-empty lines
          v.addElement(row = new double[n]);
          int j = 0;
          do {
              if (j >= n) throw new IOException("Row " + v.size() + " is too long.");
              row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
           } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
           if (j < n) throw new IOException("Row " + v.size() + " is too short.");
        }
        int m = v.size();  // Now we've got the number of rows.
        double[][] A = new double[m][];
        v.copyInto(A);  // copy the rows out of the vector
        return new Matrix(A);
      } catch (IOException e) {
        System.out.println("Error while reading a Matrix : " + e);
        return new Matrix(0,0);
      }
   }
}