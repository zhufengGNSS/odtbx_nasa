/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2009 United States Government as represented by the
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

import jat.spacetime.LeapSecondTable;
import jat.spacetime.LeapSecondTable.LeapSecData;
import jat.util.FileUtil;

import java.io.IOException;

import junit.framework.TestCase;

import org.junit.Test;

public class LeapSecondTableTest extends TestCase {

	/**
	 * Test of constructor with filename.
	 * 
	 * <ol>The following tests are run:
	 * <li>empty string</li>
	 * <li>incorrect filename</li>
	 * <li>correct filename but no data</li>
	 * <li>correct filename with data</li>
	 * </ol>
	 * 
	 * Also tests hasLeapSecs().
	 */
	@Test
	public void testLeapSecondTableString() {
		
		// TEST: empty string
		try {
			LeapSecondTable leapsectbl = new LeapSecondTable("");
			fail();
		} catch (IOException e) {
			// expected
		}
		
		// TST: incorrect filename
		try {
			LeapSecondTable leapsectbl = new LeapSecondTable("BOGUS.DAT");
			fail();
		} catch (IOException e) {
			// expected
		}
		
		// TEST: correct filename without data 
		String fs = FileUtil.file_separator();
		String directory = FileUtil.getClassFilePath("jat.spacetime.unittest","LeapSecondTableTest");
		LeapSecondTable leapsectbl;
		try {
			leapsectbl = new LeapSecondTable(directory+"emptyleapsec.dat");
			assertFalse(leapsectbl.hasLeapSecs());
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
		
		// TEST: correct filename with data
		leapsectbl = null;
		fs = FileUtil.file_separator();
		directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
		try {
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
			assertTrue(leapsectbl.hasLeapSecs());
			int maxleapsecs = leapsectbl.getAccumLeapSecs(54832.5);
			assertTrue(maxleapsecs == 34);
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
		
	}

	/**
	 * Test the copy constructor
	 */
	@Test
	public void testLeapSecondTableLeapSecondTable() {
		String fs = FileUtil.file_separator();
		String directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
		LeapSecondTable leapsectbl;
		try {
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
			assertTrue(leapsectbl.hasLeapSecs());
			
			LeapSecondTable leapcopy = new LeapSecondTable(leapsectbl);
			assertTrue(leapcopy.hasLeapSecs());
			
			double mjd = 54832.5;
			assertTrue(leapsectbl.getAccumLeapSecs(mjd) == leapcopy.getAccumLeapSecs(mjd));
			
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
		
	}

	/**
	 * Tests that the filename of a properly-loaded class is retrievable.
	 */
	@Test
	public void testGetFilename() {
		String fs = FileUtil.file_separator();
		String directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
		LeapSecondTable leapsectbl;
		try {
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
			assertTrue(leapsectbl.hasLeapSecs());
			
			String filename = leapsectbl.getFilename();
			assertTrue(filename.lastIndexOf("leapsec.dat")> 0);
			
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
	}


	/**
	 * Tests latestLeapSec().
	 * NOTE: This test will have to be updated as the leapsec.dat file changes.
	 */
	@Test
	public void testLatestLeapSec() {
		String fs = FileUtil.file_separator();
		String directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
		LeapSecondTable leapsectbl;
		try {
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
			assertTrue(leapsectbl.hasLeapSecs());
			
			// Corresponds to leapsec.dat:
			// 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0
			assertTrue(leapsectbl.latestLeapSec() == 54832);
			
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	/**
	 * Tests the ability to access data in the table via
	 * getData() and getAccumLeapSecs().
	 * 
	 * <ol>The following tests are run:
	 * <li>Tests against an empty table instance.</li>
	 * <li>dates before the table</li>
	 * <li>dates after the table</li>
	 * <li>dates on leap sec boundaries</li>
	 * <li>dates before leap sec boundaries</li>
	 * </ol>
	 * NOTE: Parts of this test may have to be updated as 
	 * the leapsec.dat file changes.
	 */
	@Test
	public void testGetData() {
		
		// TEST: against empty instance
		String fs = FileUtil.file_separator();
		String directory = FileUtil.getClassFilePath("jat.spacetime.unittest","LeapSecondTableTest");
		LeapSecondTable leapsectbl;
		try {
			leapsectbl = new LeapSecondTable(directory+"emptyleapsec.dat");
			assertFalse(leapsectbl.hasLeapSecs());
			
			LeapSecData lsd = leapsectbl.getData(42778);
			assertTrue(lsd.accumLeapSecs == 0);
			assertTrue(lsd.mjd_utc_start == 0.0);
			assertTrue(lsd.mjd_utc_next == Double.MAX_VALUE);
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
		
		/*
		 * Create a full instance for all subsequent tests.
		 * (using leapsec.dat)
		 */

		fs = FileUtil.file_separator();
		directory = FileUtil.getClassFilePath("jat.spacetime","TimeUtils");
		leapsectbl = null;
		try {
			leapsectbl = new LeapSecondTable(directory+"leapsec.dat");
			assertTrue(leapsectbl.hasLeapSecs());
			
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
		
		LeapSecData lsd;
		int ls;
		
		// TEST: dates before the table
		// leapsec.dat:
		// 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0
		double mjd = Double.MIN_VALUE;
		ls = 0;
		lsd = leapsectbl.getData(mjd);
		assertTrue(lsd.accumLeapSecs == ls);
		assertTrue(lsd.mjd_utc_start == 0.0);
		assertTrue(lsd.mjd_utc_next == 41317.0);
		assertTrue(leapsectbl.getAccumLeapSecs(mjd) == ls);
		
		mjd = 0.0;
		ls = 0;
		lsd = leapsectbl.getData(mjd);
		assertTrue(lsd.accumLeapSecs == ls);
		assertTrue(lsd.mjd_utc_start == 0.0);
		assertTrue(lsd.mjd_utc_next == 41317.0);
		assertTrue(leapsectbl.getAccumLeapSecs(mjd) == ls);

		mjd = 41317.0 - 0.0001;
		ls = 0;
		lsd = leapsectbl.getData(mjd);
		assertTrue(lsd.accumLeapSecs == ls);
		assertTrue(lsd.mjd_utc_start == 0.0);
		assertTrue(lsd.mjd_utc_next == 41317.0);
		assertTrue(leapsectbl.getAccumLeapSecs(mjd) == ls);
		
		// TEST: dates after the table
		// leapsec.dat:
		// 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0
		mjd = 60000.0;
		ls = 34;
		lsd = leapsectbl.getData(mjd);
		assertTrue(lsd.accumLeapSecs == ls);
		assertTrue(lsd.mjd_utc_start == 54832.0);
		assertTrue(lsd.mjd_utc_next == Double.MAX_VALUE);
		assertTrue(leapsectbl.getAccumLeapSecs(mjd) == ls);
		
		mjd = 54832.001;
		ls = 34;
		lsd = leapsectbl.getData(mjd);
		assertTrue(lsd.accumLeapSecs == ls);
		assertTrue(lsd.mjd_utc_start == 54832.0);
		assertTrue(lsd.mjd_utc_next == Double.MAX_VALUE);
		assertTrue(leapsectbl.getAccumLeapSecs(mjd) == ls);
		
		// TEST: dates on and before leap sec boundaries
		// from leapsec.dat, the first few leap seconds and the last two:
		double[] startdates = {41317, 41499, 41683, 53736, 54832};
		double[] nextdates  = {41499, 41683, 42048, 54832, Double.MAX_VALUE};
		int[] leaps = {10, 11, 12, 33, 34};
		
		// the same file, for testing times before the leap seconds
		double[] earlystartdates = {0.0, 41317, 41499, 51179, 53736};
		double[] earlynextdates  = startdates;
		int[] earlyleaps = {0, 10, 11, 32, 33};
		
		// test on the dates:
		for (int i = 0; i < startdates.length; i++) {
			lsd = leapsectbl.getData(startdates[i]);
			//System.out.println(lsd.toString());
			assertTrue("Failed when i="+i,lsd.accumLeapSecs == leaps[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_start == startdates[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_next == nextdates[i]);
			assertTrue("Failed when i="+i,leapsectbl.getAccumLeapSecs(startdates[i]) == leaps[i]);
		}
		
		// test just after the dates:
		for (int i = 0; i < startdates.length; i++) {
			double mjddelta = startdates[i] + (1.0/86400.0)/1000;
			lsd = leapsectbl.getData(mjddelta);
			//System.out.println(lsd.toString());
			//System.out.println(mjddelta);
			assertTrue("Failed when i="+i,lsd.accumLeapSecs == leaps[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_start == startdates[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_next == nextdates[i]);
			assertTrue("Failed when i="+i,leapsectbl.getAccumLeapSecs(mjddelta) == leaps[i]);
		}
		
		// test just before the dates:
		for (int i = 0; i < startdates.length; i++) {
			double mjddelta = startdates[i] - (1.0/86400.0)/1000;
			lsd = leapsectbl.getData(mjddelta);
			//System.out.println(lsd.toString());
			//System.out.println(mjddelta);
			assertTrue("Failed when i="+i,lsd.accumLeapSecs == earlyleaps[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_start == earlystartdates[i]);
			assertTrue("Failed when i="+i,lsd.mjd_utc_next == earlynextdates[i]);
			assertTrue("Failed when i="+i,leapsectbl.getAccumLeapSecs(mjddelta) == earlyleaps[i]);
		}
		
	}

}
