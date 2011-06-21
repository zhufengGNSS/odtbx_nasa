/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2008 United States Government as represented by the
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

package jat.spacetime.unittest;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.io.IOException;

import jat.spacetime.IAU2000;
import jat.spacetime.Time;
import jat.spacetime.IAU2000.iau80data;
import jat.util.FileUtil;

import org.junit.Test;

/**
 * JUnit test for IAU2000 methods.
 */
public class IAU2000Test {

	@Test
	public void testIau76fk5_itrf_gcrf() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testPrecess() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testNutation() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testSidereal() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testPolarm() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testFundarg() {
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test of gstime().  This tests 0.0 (4713 BC), the J2000 epoch, and a date
	 * in 2008 against the MATLAB gstime.m implementation.
	 */
	@Test
	public void testGstime() {
		double TOL = 1e-8;
		double calc = IAU2000.gstime(0.0);
		double mltruth = 4.24714537401247e+000; // MATLAB format long eng
		double diff = Math.abs(calc - mltruth);
		assertTrue("gstime(0.0) failed: " + diff + ", calc: " + calc + ", MATLAB= " + mltruth,diff < TOL);
		
		Time j2kepoch = new Time(2000, 1, 0, 0, 0, 0.0); // 2.451543500000000e+006
		calc = IAU2000.gstime(j2kepoch.jd_ut1());
		mltruth = 1.727564371525316;
		diff = Math.abs(calc - mltruth);
		assertTrue("gstime(J2000 epoch) failed: " + diff + ", calc: " + calc + ", MATLAB= " + mltruth,diff < TOL);
		
		Time testnow = new Time(2008, 8, 13, 1, 2, 3.4); // 2.454691543094907e+006
		calc = IAU2000.gstime(testnow.jd_ut1());
		mltruth = 5.887985209292070;
		diff = Math.abs(calc - mltruth);
		assertTrue("gstime(2008/8/12 1:2:3.4s) failed: " + diff + ", calc: " + calc + ", MATLAB= " + mltruth,diff < TOL);
		
		
	}

	/**
	 * Test of iau80in(). This tests:
	 * <ol>
	 * <li>Unknown input file</li>
	 * <li>An existing but empty file</li>
	 * <li>A file with incorrect data</li>
	 * <li>A file with truncated data (line truncated)</li>
	 * <li>A file with full lines, but not enough coefficients</li>
	 * <li>A proper file with full coefficients (selected coefficients).</li>
	 * </ol>
	 */
	@Test
	public void testIau80in() {
		String filename;

		// 1: unknown input file
		try {
			IAU2000.iau80in("UNKNOWNFILENAME");
		} catch (FileNotFoundException e) {
			// success
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with unknown file name.");
		}

		// 2: existing, empty file
		filename = FileUtil.getClassFilePath(this.getClass())
				+ FileUtil.file_separator() + "iau2000test_empty.dat";
		try {
			assertNull(IAU2000.iau80in(filename));
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with empty file.");
		}

		// 3: incorrect data
		filename = FileUtil.getClassFilePath(this.getClass())
				+ FileUtil.file_separator() + "iau2000test_improper.dat";
		try {
			IAU2000.iau80in(filename);
		} catch (IOException e) {
			// success
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with improper data.");
		}

		// 4: line truncation
		filename = FileUtil.getClassFilePath(this.getClass())
				+ FileUtil.file_separator() + "iau2000test_line.dat";
		try {
			assertNull(IAU2000.iau80in(filename));
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with line-terminated data.");
		}

		// 5: limited lines
		filename = FileUtil.getClassFilePath(this.getClass())
				+ FileUtil.file_separator() + "iau2000test_limited.dat";
		try {
			assertNull(IAU2000.iau80in(filename));
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with limited data.");
		}

		// 5: limited lines
		filename = FileUtil.getClassFilePath(jat.spacetime.IAU2000.class)
				+ FileUtil.file_separator() + "nut80.dat";
		iau80data data = null;
		try {
			data = IAU2000.iau80in(filename);
			assertNotNull(data);
		} catch (Exception e) {
			fail("Caught unexpected " + e.getClass().getSimpleName()
					+ " exception with limited data.");
		}
		if (data == null) {
			return; // test over
		}
		
		// contents check against the file contents
		// Remember, numbering started with 1 for both indices!
		assertTrue(data.iar80[1][1] == 0);
		assertTrue(data.iar80[1][5] == 1);
		assertTrue(data.iar80[3][3] == 2);
		assertTrue(data.iar80[4][4] == 0);
		assertTrue(data.iar80[10][2] == -1);
		
		double tol = 1E-15;
		double convrt = 0.0001 / 3600.0; // reduction factor from IAU2000
		
		assertTrue((data.rar80[1][1] - -171996.0*convrt) < tol);
		assertTrue((data.rar80[1][4] - 8.9*convrt) < tol);
		assertTrue((data.rar80[3][3] - 977.0*convrt) < tol);
		assertTrue((data.rar80[4][4] - 0.5*convrt) < tol);
		assertTrue((data.rar80[10][2] - -0.5*convrt) < tol);

	}

}
