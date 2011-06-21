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
 *
 * Emergent Space Technologies
 * File created by Rob Antonucci 
 **/
package jat.spacetime.unittest;

import jat.eph.*;
import jat.matvec.data.VectorN;
import jat.spacetime.BodyCenteredInertialRef;
import jat.spacetime.EarthRef;
import jat.spacetime.LunaFixedRef;
import jat.spacetime.LunaRef;
import jat.spacetime.ReferenceFrame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import junit.framework.TestCase;

public class LunaFixedRefTest extends TestCase {

  /** A simple structure holding the contents of one line of */
  public static class Entry
  {
    public double epochSecs;
    public VectorN lciPos;
    public VectorN lciVel;
    public VectorN lcfPos;
    public VectorN lcfVel;
  }
    
  public static void main(String[] args) {
    junit.textui.TestRunner.run(LunaFixedRefTest.class);
  }

  /*
   * Test method for 'jat.spacetime.LunaFixedRef.getTranslater(ReferenceFrame, Time)'
   */
  public void testLCItoLCF() throws IOException {
    final String COMMENT = "#";
    final String DECLARATION = "%";
    final String NUM_LINES_VAR = "num";
    final String EPOCH_VAR = "jd0";
    
    InputStream strm = getClass().getClassLoader().getResourceAsStream(
    "jat/spacetime/unittest/llo_ascii.txt");
    BufferedReader rdr = new BufferedReader(new InputStreamReader(strm));
    
    // The first lines will specify the number of data lines
    // and the start of the epoch.
    int numLines = 0;
    double epoch = 0;
    int lineCtr=0;
    while ((numLines == 0) || (epoch == 0)) {
      ++lineCtr;
      String nextLine = rdr.readLine();
      if (nextLine.startsWith(COMMENT)) {
        // It's a comment.  Ignore.
      }
      else if (nextLine.startsWith(DECLARATION)) {
        // Variable definition.  Parse.
        int equalIndex = nextLine.indexOf("=");
        assertTrue("Declaraion on line " + lineCtr + "not of the form: " +
            "% var = value", equalIndex > 2);
        String varName = nextLine.substring(1, equalIndex).trim();
        String value = nextLine.substring(equalIndex+1).trim();
        if (varName.equals(NUM_LINES_VAR)) {
          numLines = Integer.parseInt(value);
        }
        else if (varName.equals(EPOCH_VAR)) {
          epoch = Double.parseDouble(value);
        }
      }
      else if (!nextLine.trim().equals("")) {
        fail("Found data line before " + NUM_LINES_VAR +
            " and " + EPOCH_VAR + " were defined.  Error at line " + lineCtr);
      }
    }
    
    compareData(rdr, numLines, epoch, lineCtr);
  }
  
  private void compareData(BufferedReader rdr, int numLines, double epoch,
      int startLine) throws IOException
  {
    final double MARGIN = 0.01;
    BodyCenteredInertialRef inertial = new BodyCenteredInertialRef(DE405_Body.GEOCENTRIC_MOON);
    LunaFixedRef fixed = new LunaFixedRef();
    int numEntriesRead = 0;
    int lineCtr = startLine;
    boolean done = false;
    Time t = new Time(TimeUtils.JDtoMJD(epoch));
    while ((numEntriesRead < numLines) && !done) {
      Entry nextEntry = new Entry();
      ++lineCtr;
      String nextLine = rdr.readLine();      
      done = (nextLine == null);
      if ((nextLine != null) && !nextLine.trim().equals(""))
      {
        parseNext(nextLine, nextEntry, lineCtr);
        t.update(nextEntry.epochSecs);
        
        // Verify the point translations
        ReferenceFrameTranslater xlater = 
          new ReferenceFrameTranslater(inertial, fixed, t);
        VectorN computedF = xlater.translatePoint(nextEntry.lciPos);
        double distance = computedF.minus(nextEntry.lcfPos).mag();
        double allowed = nextEntry.lcfPos.mag() * MARGIN;
        assertTrue("Line " + lineCtr + " position: " + computedF + " is off from " + 
            nextEntry.lcfPos + " by " + distance,  distance <= allowed);
        ReferenceFrameTranslater backXlater = 
          new ReferenceFrameTranslater(fixed, inertial, t);
        VectorN computedI = backXlater.translatePoint(nextEntry.lcfPos);
        distance = computedI.minus(nextEntry.lciPos).mag();
        allowed = nextEntry.lciPos.mag() * MARGIN;
        assertTrue("Line " + lineCtr + " position: " + computedI + " is off from " + 
            nextEntry.lciPos + " by " + distance,  distance <= allowed);
        VectorN recomputedI = xlater.translatePointBack(nextEntry.lcfPos);
        distance = recomputedI.minus(nextEntry.lciPos).mag();
        assertTrue("Line " + lineCtr + " position: " + recomputedI + " is off from " + 
            nextEntry.lciPos + " by " + distance,  distance <= allowed);
        
        // Verify the velocity translations
        computedF = xlater.translateVelocity(nextEntry.lciVel, nextEntry.lciPos);
        double difference = computedF.minus(nextEntry.lcfVel).mag();
        allowed = nextEntry.lcfVel.mag() * MARGIN;
        assertTrue("Line " + lineCtr + " velocity: " + computedF + " is off from " + 
            nextEntry.lcfVel + " by " + difference,  difference <= allowed);
        computedI = backXlater.translateVelocity(nextEntry.lcfVel, nextEntry.lcfPos);
        difference = computedI.minus(nextEntry.lciVel).mag();
        allowed = nextEntry.lciVel.mag() * MARGIN;
        assertTrue("Line " + lineCtr + " velocity: " + computedI + " is off from " + 
            nextEntry.lciVel + " by " + difference,  difference <= allowed);
        recomputedI = xlater.translateVelocityBack(nextEntry.lcfVel, nextEntry.lcfPos);
        difference = recomputedI.minus(nextEntry.lciVel).mag();
        assertTrue("Line " + lineCtr + " velocity: " + recomputedI + " is off from " + 
            nextEntry.lciVel + " by " + difference,  difference <= allowed);
        
        ++numEntriesRead;
      }
    }
    assertEquals("Encountered " + numEntriesRead + " lines of data when " +
        numLines + " were expected.", numLines, numEntriesRead);
  }
  
  private void parseNext(String input, Entry entry, int lineNum)
  {
    try {
      // The line should have 13 doubles separated by spaces.
      String[] doubles = input.trim().split("\\s+");
      assertEquals("Line " + String.valueOf(lineNum) + " needs 13 number separated " +
          "by spaces.", 13, doubles.length);
      NumberFormat parser = new DecimalFormat("0.0000000000000000E000");
      entry.epochSecs = Double.parseDouble(doubles[0]);
      entry.lcfPos = new VectorN(3);
      entry.lcfPos.set(0, Double.parseDouble(doubles[1]));
      entry.lcfPos.set(1, Double.parseDouble(doubles[2]));
      entry.lcfPos.set(2, Double.parseDouble(doubles[3]));
      entry.lcfVel = new VectorN(3);
      entry.lcfVel.set(0, Double.parseDouble(doubles[4]));
      entry.lcfVel.set(1, Double.parseDouble(doubles[5]));
      entry.lcfVel.set(2, Double.parseDouble(doubles[6]));
      entry.lciPos = new VectorN(3);
      entry.lciPos.set(0, Double.parseDouble(doubles[7]));
      entry.lciPos.set(1, Double.parseDouble(doubles[8]));
      entry.lciPos.set(2, Double.parseDouble(doubles[9]));
      entry.lciVel = new VectorN(3);
      entry.lciVel.set(0, Double.parseDouble(doubles[10]));
      entry.lciVel.set(1, Double.parseDouble(doubles[11]));
      entry.lciVel.set(2, Double.parseDouble(doubles[12]));
    }
    catch (NullPointerException e) {
      fail("Failure to parse doubles on line " + lineNum + ". " +
          e.getMessage());
    }
  }
  
  public void testECIVelocityToLCFVelocity() {
    //double[] startEciPos = {-2.7555712051348615E8, 2.278862267130929E8, 1.0302526370467307E8};
    //double[] startEciVel = {719.703479618275,  737.5548929186978,   220.44009064585677};
    double[] startEciPos = {1.7232607076253737e+002,  1.6379287856614317e+003, -8.1592709137349675e+002};
    double[] startEciVel = {1.3836671664810035e-001,  7.1393322649210633e-001,  1.4624051437086556e+000};
    double timestep = 0.5; // in seconds
    double ERROR = 0.01;
    
    Time t = new Time(TimeUtils.JDtoMJD(2.458232499999999e+006));
    ReferenceFrame eciRef = new LunaRef();
    ReferenceFrame lcfRef = new LunaFixedRef();
    
    VectorN eciPos1 = new VectorN(startEciPos);
    VectorN eciVel1 = new VectorN(startEciVel);
    System.out.println("ECI Position =\t\t" + eciPos1);
    System.out.println("ECI Velocity =\t\t" + eciVel1);
    ReferenceFrameTranslater xlater = new ReferenceFrameTranslater(eciRef, lcfRef, t);
    VectorN lcfPos1 = xlater.translatePoint(eciPos1);
    VectorN lcfVel1 = xlater.translateVelocity(eciVel1, eciPos1);
    System.out.println("LCF Position =\t\t" + lcfPos1);
    System.out.println("LCF Velocity =\t\t" + lcfVel1);
    
    t.update(timestep);
    VectorN eciPos2 = eciVel1.times(timestep).plus(eciPos1);
    VectorN lcfPos2 = lcfVel1.times(timestep).plus(lcfPos1);
    System.out.println("ECI New Position =\t" + eciPos2);
    xlater = new ReferenceFrameTranslater(eciRef, lcfRef, t);
    VectorN lcfPos2True = xlater.translatePoint(eciPos2);
    System.out.println("LCF New Position =\t" + lcfPos2True);
    System.out.println("LCF Estimated =\t\t" + lcfPos2);
    
    VectorN change = lcfPos2.minus(lcfPos1);
    VectorN changeTrue = lcfPos2True.minus(lcfPos1);
    VectorN error = lcfPos2.minus(lcfPos2True);
    System.out.println("LCF Position changed by\t" + changeTrue);
    System.out.println("LCF estimated to change\t" + change);
    System.out.println("Error was\t\t" + error);
    assertTrue("velocity produced change " + change + " instead of " +
        changeTrue, error.mag() < changeTrue.mag() * ERROR);
    
    
  }

}
