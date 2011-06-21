/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
package jat.eph.unittest;

import jat.eph.DE405;
import jat.eph.DE405_Body;
import jat.matvec.data.VectorN;
import jat.util.FileUtil;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

/**
 * JUnit tests for jat.eph.DE405.
 */
public class DE405Test extends TestCase {

	public static void main(String[] args) {
		junit.textui.TestRunner.run(DE405Test.class);
	}

	/**
	 * A test that verifies the ability to detect the required data files.
	 */
	public void testDataFiles() {

		// a completely bogus directory
		try {
			new DE405("BogusPathname");
			fail("The bogus pathname should have failed.");
		} catch (java.util.MissingResourceException ex) {
			// expected
		}

		// an existing directory, but with no data files
		String filesep = FileUtil.file_separator();
		String directory;
		try {
			directory = FileUtil.getClassFilePath("jat.eph", "DE405");
		} catch (Exception e) {
			System.err.println("Error: Could not read default DE405 path.");
			directory = "C:/Code/Jat/jat/eph/";
		}
		try {
			final String p = directory + filesep;
			new DE405(p);
			fail("The empty pathname: "
					+ p
					+ " should have failed. You don't have DE405 files there, do you?");
		} catch (java.util.MissingResourceException ex) {
			// expected
		}

	}

	/**
	 * A test to check the allowable DE405 time ranges.
	 * Note that floating-point representation and calculations
	 * may make some inputs fail even though they are typed 
	 * exactly.
	 */
	public void testTimeRange() {

		DE405 ephem = new DE405();

		double t = 0;
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
			fail("Out-of-range time " + t + " should have failed.");
		}
		catch (IllegalArgumentException ex) {
			// expected
		}

		t = 14991;
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
			fail("Out-of-range time " + t + " should have failed.");
		}
		catch (IllegalArgumentException ex) {
			// expected
		}

		// NOTE: t = 14992 should pass, but floating point 
		// math made it fail!

		t = 14992;
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
			fail("Out-of-range time " + t + " should have failed.");
		}
		catch (IllegalArgumentException ex) {
			// expected
		}

		t = 14992.0001; // this should pass
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
		}
		catch (IllegalArgumentException ex) {
			fail("In-range time " + t + " should not have failed.");
		}

		t = 124623.99999; // this should pass
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
		}
		catch (IllegalArgumentException ex) {
			fail("In-range time " + t + " should not have failed.");
		}

		// NOTE, t = 124624.0 should pass but floating-point math 
		// made it fail.

		t = 124624.0001; // this should also fail
		try {
			ephem.get_planet_pos(DE405_Body.EARTH, t);
			fail("Out-of-range time " + t + " should have failed.");
		}
		catch (IllegalArgumentException ex) {
			// expected
		}

	}

	/** Tests DE405 positions and velocities against data from SPICE DE405.
	 * Reads in pre-computed SPICE DE405 data for each DE405_Body combination
	 * and compares it to a pre-calculated tolerance.  The tolerance is from
	 * a comparison between JAT and SPICE and may include a wider tolerance
	 * for other effects not related to DE405 (such as time handling).
	 */
	public void testPosVel() {
		DE405 ephem = new DE405();

		// the SPICE body file name portions, matched with the DE405_Body, below
		List<String> p = new ArrayList<String>();
		p.add("MERCURY BARYCENTER"); 
		p.add("VENUS BARYCENTER");
		p.add("EARTH");
		p.add("MARS BARYCENTER"); 
		p.add("JUPITER BARYCENTER"); 
		p.add("SATURN BARYCENTER");
		p.add("URANUS BARYCENTER");
		p.add("NEPTUNE BARYCENTER"); 
		p.add("PLUTO BARYCENTER");
		p.add("SUN");
		p.add("MOON"); 
		p.add("EARTH BARYCENTER"); 
		p.add("SUN");
		p.add("MOON"); 
		p.add("SOLAR SYSTEM BARYCENTER");

		DE405_Body[] bod = {
				DE405_Body.MERCURY,
				DE405_Body.VENUS,
				DE405_Body.EARTH,
				DE405_Body.MARS, 
				DE405_Body.JUPITER,
				DE405_Body.SATURN,
				DE405_Body.URANUS,
				DE405_Body.NEPTUNE,
				DE405_Body.PLUTO,
				DE405_Body.GEOCENTRIC_SUN,
				DE405_Body.GEOCENTRIC_MOON,
				DE405_Body.EM_BARY,
				DE405_Body.SUN,
				DE405_Body.MOON,
				DE405_Body.SOLAR_SYSTEM_BARY };

		// read the tolerance data
		List<double[]> tollist = readTargetPosVels("TOL");
		assertTrue(tollist.size() > 0);

		double TOL_SF = 1.1; // "bump" to the tolerance data due to reading from a file

		for (int planet = 0; planet < 6; ++planet) {
			// tolerance data, use indices 1-6 for position & velocity
			double tol[] = (double[]) tollist.get(planet); 

			List<double[]> l = readTargetPosVels(p.get(planet));
			assertTrue(l.size() > 0);

			for (int d = 0; d < l.size(); d++) {

				double spicedat[] = (double[]) l.get(d); 
				VectorN posVel = ephem.get_planet_posvel(bod[planet], spicedat[0]);

				// Position check
				for(int index=0; index<3; ++index) {
					if (Math.abs(posVel.get(index)-spicedat[index+1]) > (tol[index+1]*TOL_SF)) {
						String msg = "DE405Test.java failed on:"
							+ "planet " + p.get(planet) + ", index " + index
							+ ", JAT: " + posVel.get(index)
							+ ", SPICE: " + spicedat[index+1] 
							                         + ", diff: " + Math.abs(posVel.get(index)-spicedat[index+1]) 
							                         + ", (tol " + tol[index+1] + ")\n";
						fail(msg);
					}
				}
				// Velocity check
				for(int index=3; index<6; ++index) {
					if (Math.abs(posVel.get(index)-spicedat[index+1]) > (tol[index+1]*TOL_SF)) {
						String msg = "DE405Test.java failed on:"
							+ "planet " + p.get(planet) + ", index " + index
							+ ", JAT: " + posVel.get(index)
							+ ", SPICE: " + spicedat[index+1] 
							                         + ", diff: " + Math.abs(posVel.get(index)-spicedat[index+1]) 
							                         + ", (tol " + tol[index+1] + ")\n";
						fail(msg);
					}
				}
			}
		}
	}

	/**
	 * The test for DE405 moon libration angles (needs to be implemented).
	 */
	public void testLibration() {
		// TODO - needs to be implemented
	}


	/**
	 * Reads a space-delimited file that contains SPICE data with the following
	 * format, 7 items per line, time position and velocity: 
	 * time (TT, mjd) r_x r_y r_z (m) v_x v_y v_z (m/s) 
	 * @param planetname The planet name part to fill in the filename
	 * @return List of double[7]
	 */
	private List<double[]> readTargetPosVels(String planetname) {

		List<double[]> l = new ArrayList<double[]>();
		String filename = "jat/eph/unittest/de405_" + planetname + "_data.dat";
		try {
			InputStream rstrm = getClass().getClassLoader().
			getResourceAsStream(filename);
			if (rstrm == null) {
				fail("readTargetPosVels() Failed to open file: " + filename);
			}
			else {
				BufferedReader rrdr = new BufferedReader(new InputStreamReader(rstrm));
				String rLine = rrdr.readLine();

				while (true) {
					rLine = rrdr.readLine();
					if (rLine == null) {
						break;  // done reading
					}
					else if (rLine.trim().startsWith("#")) {
						continue; // skip comments
					}
					else {
						String[] numberStrs = rLine.trim().split("\\s+");
						double[] d = new double[7];
						for(int index=0; index<7; ++index) {
							d[index] = Double.parseDouble(numberStrs[index]);
						}
						l.add(d);
					}
				}
			}
		}
		catch (IOException ioe) {
			System.out.println("Error reading file: " + filename + ", error was: " + ioe.getMessage());
			l.clear();
		}

		return l;
	}
}
