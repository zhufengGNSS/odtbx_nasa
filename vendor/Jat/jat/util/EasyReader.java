/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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
 *
 */

package jat.util;

import java.io.*;

/**
 *
 *  EasyReader provides simple methods for reading the console and
 *  for opening and reading text files.  All exceptions are handled
 *  inside the class and are hidden from the user.
 *
 *  <xmp>
 *  Example:
 *  =======
 *
 *  EasyReader console = new EasyReader();
 *  System.out.print("Enter input file name: ");
 *  String fileName = console.readLine();
 *
 *  EasyReader inFile = new EasyReader(fileName);
 *  if (inFile.bad())
 *  {
 *    System.err.println("Can't open " + fileName);
 *    System.exit(1);
 *  }
 *
 *  String firstLine = inFile.readLine();
 *  if (!inFile.eof())   // or:  if (firstLine != null)
 *    System.out.println("The first line is : " + firstLine);
 *
 *  System.out.print("Enter the maximum number of integers to read: ");
 *  int maxCount = console.readInt();
 *  int k, count = 0;
 *
 *  while (count < maxCount && !inFile.eof())
 *  {
 *    k = inFile.readInt();
 *    if (!inFile.eof())
 *    {
 *      // process or store this number
 *      count++;
 *    }
 *  }
 *
 *  inFile.close();    // optional
 *  System.out.println(count + " numbers read");
 *  </xmp>
 *
 */

public class EasyReader
{
  protected String myFileName;
  protected BufferedReader myInFile;
  protected int myErrorFlags = 0;
  protected static final int OPENERROR = 0x0001;
  protected static final int CLOSEERROR = 0x0002;
  protected static final int READERROR = 0x0004;
  protected static final int EOF = 0x0100;

  /**
   *  Constructor.  Prepares console (System.in) for reading
   */
  public EasyReader()
  {
    myFileName = null;
    myErrorFlags = 0;
    myInFile = new BufferedReader(
                            new InputStreamReader(System.in), 128);
  }

  /**
   *  Constructor.  opens a file for reading
   *  @param fileName the name or pathname of the file
   */
  public EasyReader(String fileName)
  {
    myFileName = fileName;
    myErrorFlags = 0;
    try
    {
      myInFile = new BufferedReader(new FileReader(fileName), 1024);
    }
    catch (FileNotFoundException e)
    {
      myErrorFlags |= OPENERROR;
      myFileName = null;
    }
  }

  /**
   *  Closes the file
   */
  public void close()
  {
    if (myFileName == null)
      return;
    try
    {
      myInFile.close();
    }
    catch (IOException e)
    {
      System.err.println("Error closing " + myFileName + "\n");
      myErrorFlags |= CLOSEERROR;
    }
  }

  /**
   *  Checks the status of the file
   *  @return true if en error occurred opening or reading the file,
   *  false otherwise
   */
  public boolean bad()
  {
    return myErrorFlags != 0;
  }

  /**
   *  Checks the EOF status of the file
   *  @return true if EOF was encountered in the previous read
   *  operation, false otherwise
   */
  public boolean eof()
  {
    return (myErrorFlags & EOF) != 0;
  }

  private boolean ready() throws IOException
  {
    return myFileName == null || myInFile.ready();
  }

  /**
   *  Reads the next character from a file (any character including
   *  a space or a newline character).
   *  @return character read or <code>null</code> character
   *  (Unicode 0) if trying to read beyond the EOF
   */
  public char readChar()
  {
    char ch = '\u0000';

    try
    {
      if (ready())
      {
         ch = (char)myInFile.read();
      }
    }
    catch (IOException e)
    {
      if (myFileName != null)
        System.err.println("Error reading " + myFileName + "\n");
      myErrorFlags |= READERROR;
    }

    if (ch == '\u0000')
      myErrorFlags |= EOF;

    return ch;
  }

  /**
   *  Reads from the current position in the file up to and including
   *  the next newline character.  The newline character is thrown away
   *  @return the read string (excluding the newline character) or
   *  null if trying to read beyond the EOF
   */
  public String readLine()
  {
    String s = null;

    try
    {
      s = myInFile.readLine();
    }
    catch (IOException e)
    {
      if (myFileName != null)
        System.err.println("Error reading " + myFileName + "\n");
      myErrorFlags |= READERROR;
    }

    if (s == null)
      myErrorFlags |= EOF;
    return s;
  }

  /**
   *  Skips whitespace and reads the next word (a string of consecutive
   *  non-whitespace characters (up to but excluding the next space,
   *  newline, etc.)
   *  @return the read string or null if trying to read beyond the EOF
   */
  public String readWord()
  {
    StringBuffer buffer = new StringBuffer(128);
    char ch = ' ';
    int count = 0;
    String s = null;

    try
    {
      while (ready() && Character.isWhitespace(ch))
        ch = (char)myInFile.read();
      while (ready() && !Character.isWhitespace(ch))
      {
        count++;
        buffer.append(ch);
        myInFile.mark(1);
        ch = (char)myInFile.read();
      };

      if (count > 0)
      {
        myInFile.reset();
        s = buffer.toString();
      }
      else
      {
        myErrorFlags |= EOF;
      }
    }

    catch (IOException e)
    {
      if (myFileName != null)
        System.err.println("Error reading " + myFileName + "\n");
      myErrorFlags |= READERROR;
    }

    return s;
  }

  /**
   *  Reads the next integer (without validating its format)
   *  @return the integer read or 0 if trying to read beyond the EOF
   */
  public int readInt()
  {
    String s = readWord();
    if (s != null)
      return Integer.parseInt(s);
    else
      return 0;
  }

  /**
   *  Reads the next double (without validating its format)
   *  @return the number read or 0 if trying to read beyond the EOF
   */
  public double readDouble()
  {
	String s = readWord();
	if (s != null)
	  return Double.parseDouble(s);
	  // in Java 1, use: return Double.valueOf(s).doubleValue();
	else
	  return 0.0;
  }

  /**
   *  Reads the next double (without validating its format)
   *  FORTRAN exponential format has D instead of E
   *  @return the number read or 0 if trying to read beyond the EOF
   */
  public double readFORTRANDouble()
  {
	String s = readWord();
	s=s.replace('D','E');
	if (s != null)
	  return Double.parseDouble(s);
	  // in Java 1, use: return Double.valueOf(s).doubleValue();
	else
	  return 0.0;
  }
}

