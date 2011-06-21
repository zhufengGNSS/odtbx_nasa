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

package jat.forces.unittest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import jat.forces.ForceModel;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.EarthRef;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import junit.framework.TestCase;

public class ForceModelTest extends TestCase {
  
  /** Used to indicate comments in the test data files */
  public static final String COMMENT_PREFIX = "#";
  
  /** Acceptable error margin.  Error should only result from
   * floating point arithmetic, so errors should be VERY small. */
  private static double MARGIN = 0.000001;

  public static void main(String[] args) {
    junit.textui.TestRunner.run(ForceModelTest.class);
  }

  /*
   * This is a common test method for testing force models.
   * Relies on base force model testers to provide force
   * and spacecraft.
   */
  protected void testForceModelAcceleration(Spacecraft sc, ForceModel force,
      String targetFilename, String targetName) throws IOException {
    // We assume the positions are ECI
    EarthRef ref = null;
    
    InputStream rstrm = getClass().getClassLoader().
      getResourceAsStream("jat/forces/unittest/leo_orbit.txt");
    BufferedReader rrdr = new BufferedReader(new InputStreamReader(rstrm));
    String rLine = rrdr.readLine();
    // First lines may be comments and can be skipped.
    while (rLine.trim().startsWith(COMMENT_PREFIX)) {
      rLine = rrdr.readLine();
    }
    
    String forceFile = "jat/forces/unittest/" + targetFilename;
    InputStream fstrm = 
      getClass().getClassLoader().getResourceAsStream(forceFile);
    BufferedReader frdr = null;
    BufferedWriter fwrtr = null;
    String fLine = null;
    if (fstrm == null) {
      // We cannot read the target file.  Create a target file instead.
      System.out.println("Could not open test data file " +
          forceFile);
      File outFile = new File(System.getProperty("java.io.tmpdir"),
          targetFilename);
      fwrtr = new BufferedWriter(new FileWriter(outFile));
      System.out.println("Creating " + targetName + " file " +
          outFile.getAbsolutePath());
      fwrtr.write(COMMENT_PREFIX + "\t" + targetName);
      fwrtr.newLine();
      fwrtr.write(COMMENT_PREFIX + "\tx\t\t\t\ty\t\t\t\tz");
      fwrtr.newLine();
    }
    else {
      frdr = new BufferedReader(new InputStreamReader(fstrm));
      fLine = frdr.readLine();
      // First lines may be comments and can be skipped.
      while (fLine.trim().startsWith(COMMENT_PREFIX)) {
        fLine = frdr.readLine();
      }
    }
  
    int lineNum = 1;
    Time time = null;
    while (rLine != null) {
      VectorN position = new VectorN(3);
      VectorN velocity = new VectorN(3);
      Time newTime = parse(rLine, position, velocity);
      if (time == null) {
        time = newTime;
        ref = new EarthRef(time.mjd_ut1(), time.mjd_tt());
        ref.initializeSunEphem(time.mjd_tt());
        ref.set_use_sun(true);
      }
      else {
        double seconds = (newTime.mjd_utc() - time.mjd_utc()) * TimeUtils.sec2days;
        time.update(seconds);
        ref.update(time);
      }
      sc.updateMotion(position, velocity);
      VectorN computed = force.acceleration(time, ref, sc);
      if (fwrtr == null) {
        // Doing compare
        VectorN target = readTarget(fLine);
        double distance = target.minus(computed).mag();
        double allowed = target.mag() * MARGIN;
        if (distance > allowed) {
        	handleError(targetName, lineNum, computed, target);
        }
        else {
        	handlePass(targetName, lineNum, computed, target);
        }
        assertTrue("Error computing " + targetName + " at " +
            "time/position specified on line " + lineNum + ".  " + computed + 
            " is off from " + target + " by " + distance, distance <= allowed);
        fLine = frdr.readLine();
      }
      else {
        // Doing generation of target file
        writeTarget(computed, fwrtr);
      }
      rLine = rrdr.readLine();
      ++lineNum;
    }
    
    if (fwrtr != null) {
      fwrtr.close();
      fail("Could not find target file to which to compare computed results.");
    }
  }
  
  private void handleError(String targetName, int lineNum, VectorN computed,
		  VectorN target) {
      double distance = target.minus(computed).mag();
      double allowed = target.mag() * MARGIN;
      assertTrue("Error computing " + targetName + " at " +
              "time/position specified on line " + lineNum + ".  " + computed + 
              " is off from " + target + " by " + distance, distance <= allowed);

  }

  private void handlePass(String targetName, int lineNum, VectorN computed,
		  VectorN target) {
      // Do nothing.

  }

  /**
   * Parses the input orbit file containing position and velocity in km & km/s,
   * with time in MJD UTC.
   * @param positionInput
   * @param position (m)
   * @param velocity (m/s)
   * @return Time 
   */
  private Time parse(String positionInput, VectorN position, VectorN velocity) {
    // First parse the line into four doubles.
    String[] numberStrs = positionInput.trim().split("\\s+");
    double t = Double.parseDouble(numberStrs[0]);
    Time time = new Time(t);
    position.set(0, Double.parseDouble(numberStrs[1])*1000);
    position.set(1, Double.parseDouble(numberStrs[2])*1000);
    position.set(2, Double.parseDouble(numberStrs[3])*1000);
    if (numberStrs.length >= 7) {
      velocity.set(0, Double.parseDouble(numberStrs[4])*1000);
      velocity.set(1, Double.parseDouble(numberStrs[5])*1000);
      velocity.set(2, Double.parseDouble(numberStrs[6])*1000);
    }
    return time;
  }

  private VectorN readTarget(String positionInput) {
    // First parse the line into three doubles.
    String[] numberStrs = positionInput.trim().split("\\s+");
    VectorN force = new VectorN(3);
    force.set(0, Double.parseDouble(numberStrs[0]));
    force.set(1, Double.parseDouble(numberStrs[1]));
    force.set(2, Double.parseDouble(numberStrs[2]));
    return force;
  }

  private void writeTarget(VectorN force, BufferedWriter wrtr)
    throws IOException {
    wrtr.write(force.get(0) + "\t" + force.get(1) + "\t" + force.get(2));
    wrtr.newLine();
  }
}
