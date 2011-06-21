/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2010 United States Government as represented by the
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

package jat.spacetime.unittest;

import jat.matvec.data.VectorN;
import jat.spacetime.EarthRef;
import jat.spacetime.Time;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import junit.framework.TestCase;

public class EarthFixedRefTest extends TestCase {

  /** A simple structure holding the contents of one line of */
  public static class Entry
  {
    public Time time;
    public VectorN eciPos;
    public VectorN eciVel;
    public VectorN ecfPos;
    public VectorN ecfVel;
  }
  
  private int numEntries = 0;
  private double totalPosErr = 0.0;
  private double totalVelErr = 0.0;
    
  public static void main(String[] args) throws IOException {
	  if ((args.length == 1) && args[0].equals("stats")) {
		  EarthFixedRefTest test = new EarthFixedRefTest();
		  test.testECItoECEF();
		  test.reportStats();
	  }
	  else {
		  junit.textui.TestRunner.run(EarthFixedRefTest.class);
	  }
}

  public void reportStats() {
	  DecimalFormat fmt = new DecimalFormat("##0.00%");
	  System.out.println("Average position error = " + fmt.format(totalPosErr/numEntries));
	  System.out.println("Average velocity error = " + fmt.format(totalVelErr/numEntries));
  }
  /*
   * Test method for 'jat.spacetime.EarthFixedRef.getTranslater(ReferenceFrame, Time)'
   */
  public void testECItoECEF() throws IOException {
    
    String testDataFile = "jat/spacetime/unittest/elo_ascii.txt";
    InputStream strm = getClass().getClassLoader().getResourceAsStream(
    testDataFile);
    BufferedReader rdr = new BufferedReader(new InputStreamReader(strm));
    
    EarthRef inertial = null;
    int numEntriesRead = 0;
    int lineCtr = 0;
    boolean done = false;
    boolean parsedSomething = false;
    while (!done) {
      Entry nextEntry = new Entry();
      ++lineCtr;
      String nextLine = rdr.readLine();      
      done = (nextLine == null);
      if ((nextLine != null) && !nextLine.startsWith("#") &&
    		  !nextLine.trim().equals(""))
      {
        parseNext(nextLine, nextEntry, lineCtr);
        parsedSomething = true;
        inertial = new EarthRef(nextEntry.time);
        // Compute the translations and errors
        VectorN xlated = 
        	inertial.eci2ecf(nextEntry.eciPos, nextEntry.eciVel, nextEntry.time);
        VectorN r_ecef = xlated.get(0, 3);
        VectorN v_ecef = xlated.get(3, 3);
        
        double rErr = r_ecef.minus(nextEntry.ecfPos).mag()/nextEntry.ecfPos.mag();
        double vErr = v_ecef.minus(nextEntry.ecfVel).mag()/nextEntry.ecfVel.mag();
        handle(nextEntry, r_ecef, v_ecef, rErr, vErr, lineCtr);        
        ++numEntriesRead;
      }
    }
    assertTrue("Could not parse any data from " + testDataFile, parsedSomething);
  }
  
  private void handle(Entry entry, VectorN r_ecef, VectorN v_ecef, double rErr,
		  double vErr, int line) {
	  final double MARGIN = 0.00001;
	  ++numEntries;
	  totalPosErr += rErr;
	  totalVelErr += vErr;
	  if ((rErr > MARGIN) || (vErr > MARGIN)) {
		  DecimalFormat fmt = new DecimalFormat("##0.0000000%");
		  assertTrue("Line " + line + " position: " + r_ecef + " is off from " + 
	              entry.ecfPos + " by " + fmt.format(rErr),  rErr <= MARGIN);
		  assertTrue("Line " + line + " velocity: " + v_ecef + " is off from " + 
	              entry.ecfVel + " by " + fmt.format(vErr),  vErr <= MARGIN);
	  }
  }
  
  private void parseNext(String input, Entry entry, int lineNum)
  {
    try {
      // The line should have 14 doubles separated by spaces.
      String[] doubles = input.trim().split("\\s+");
      assertEquals("Line " + String.valueOf(lineNum) + " needs 14 number separated " +
          "by spaces.", 14, doubles.length);
      
      // Time comes in two columns, days and seconds
      entry.time = new Time(Integer.parseInt(doubles[0]));
      entry.time.update(Double.parseDouble(doubles[1]));
      
      entry.eciPos = new VectorN(3);
      entry.eciPos.set(0, Double.parseDouble(doubles[2]));
      entry.eciPos.set(1, Double.parseDouble(doubles[3]));
      entry.eciPos.set(2, Double.parseDouble(doubles[4]));
      entry.eciVel = new VectorN(3);
      entry.eciVel.set(0, Double.parseDouble(doubles[5]));
      entry.eciVel.set(1, Double.parseDouble(doubles[6]));
      entry.eciVel.set(2, Double.parseDouble(doubles[7]));
      entry.ecfPos = new VectorN(3);
      entry.ecfPos.set(0, Double.parseDouble(doubles[8]));
      entry.ecfPos.set(1, Double.parseDouble(doubles[9]));
      entry.ecfPos.set(2, Double.parseDouble(doubles[10]));
      entry.ecfVel = new VectorN(3);
      entry.ecfVel.set(0, Double.parseDouble(doubles[11]));
      entry.ecfVel.set(1, Double.parseDouble(doubles[12]));
      entry.ecfVel.set(2, Double.parseDouble(doubles[13]));
    }
    catch (NullPointerException e) {
      fail("Failure to parse doubles on line " + lineNum + ". " +
          e.getMessage());
    }
  }
  
}
