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
 * 
 */
package jat.spacetime.unittest;

import jat.spacetime.FitIERS;
import jat.util.FileUtil;
import junit.framework.TestCase;

/**
 * JUnit test for {@link jat.spacetime.FitIERS}
 * <br>
 * This test tries to document some of the unexpected side-effects
 * of the current static cache design in FitIERS.
 * @author abrown
 *
 */
public class FitIERSTest extends TestCase {
	
	protected static final double filetol = 1e-15;
	
	protected static final boolean enableCacheTest = true;

	/**
	 * Test method for {@link jat.spacetime.FitIERS#FitIERS()}.
	 * Ensures the constructor works.
	 */
	public void testFitIERS() {
		FitIERS f = null;
		try {
			f = new FitIERS();
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(): " + t.toString());
		}
		assertNotNull(f);
		
		assertTrue(f.getCacheHits() == 0);
		assertTrue(f.getCacheMisses() == 0);
	}

	/**
	 * Test method for {@link jat.spacetime.FitIERS#FitIERS(java.lang.String)}.
	 * Tests against the same file that is built into the constructor and 
	 * against a bogus (non-existent) filename.
	 */
	public void testFitIERSString() {

		// Test against a vaild file
        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = null;
		try {
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(String): " + t.toString());
		}
		assertNotNull(f);
		
		assertTrue(f.getCacheHits() == 0);
		assertTrue(f.getCacheMisses() == 0);

		/*
		 * Test against an invalid filename: the prior cached data isn't changed.
		 * THIS MAY BE AN UNDESIRABLE DESIGN! see below
		 */
		filename = directory+"BOGUS.dat";
		
		f = null;
		boolean gotThrowable = false;
		try {
			System.err.println("THE FOLLOWING STACKDUMP OF java.io.FileNotFoundException ON STDERR IS EXPECTED: \n");
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			// expected
			gotThrowable = true;
		}
		// Unfortunately, this class design doesn't throw an 
		// exception with a bad constructor argument:
		assertFalse(gotThrowable);
		
		// And we still get a usable class...
		assertNotNull(f);
		
		assertTrue(f.getCacheHits() == 0);
		assertTrue(f.getCacheMisses() == 0);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		// That has the old data...
		// EOP.dat, before line 1: clamp to line 1 data
		//41683.00
		double[] l0 = f.search(41683.00);
		assertTrue(l0.length == 3);
		assertTrue(StrictMath.abs(l0[0] - 0.120724) < filetol);
		
		assertTrue(f.getCacheHits() == 0);
		assertTrue(f.getCacheMisses() == 0);
	}

	/**
	 * Test method for {@link jat.spacetime.FitIERS#search(double)}.
	 * Tests the search method against a known file.
	 * Tests against exact times and in-between times for proper linear 
	 * interpolation.
	 * Tests clamping of data for out-of-range numbers as well.
	 */
	public void testSearch() {
		
		// for tracking our expectations of the cache operation
		int myhitcount = 0;
		int mymisscount = 0;

        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = new FitIERS(filename);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, before line 1: clamp to line 1 data
		//41683.00
		double[] l0 = f.search(41683.00);
		assertTrue(l0.length == 3);
		assertTrue(StrictMath.abs(l0[0] - 0.120724) < filetol);
		assertTrue(StrictMath.abs(l0[1] - 0.137066) < filetol);
		assertTrue(StrictMath.abs(l0[2] - 0.8084319) < filetol);

		// cache check is bypassed by out-of-range checks
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat line 1:
		//41684.00   .120724   .009786   .137066   .015902   .8084319   .0002710
		double[] l1 = f.search(41684.00);
		assertTrue(l1.length == 3);
		assertTrue(StrictMath.abs(l1[0] - 0.120724) < filetol);
		assertTrue(StrictMath.abs(l1[1] - 0.137066) < filetol);
		assertTrue(StrictMath.abs(l1[2] - 0.8084319) < filetol);

		myhitcount++;  // cached 41684.0 from 41683.0 (before table)
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat line 1: again, to check caching
		//41684.00   .120724   .009786   .137066   .015902   .8084319   .0002710
		double[] l1b = f.search(41684.00);
		assertTrue(l1b.length == 3);
		assertTrue(StrictMath.abs(l1b[0] - 0.120724) < filetol);
		assertTrue(StrictMath.abs(l1b[1] - 0.137066) < filetol);
		assertTrue(StrictMath.abs(l1b[2] - 0.8084319) < filetol);

		myhitcount++;  // dupe time
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		//EOP.dat line 14:
		//41697.00   .097304   .020404   .120936   .014510   .7657376   .0003986
		double[] l14 = f.search(41697.00);
		assertTrue(l14.length == 3);
		assertTrue(StrictMath.abs(l14[0] - 0.097304) < filetol);
		assertTrue(StrictMath.abs(l14[1] - 0.120936) < filetol);
		assertTrue(StrictMath.abs(l14[2] - 0.7657376) < filetol);
		
		mymisscount++; // expected cache miss
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat line 1: yet again, to check caching
		//41684.00   .120724   .009786   .137066   .015902   .8084319   .0002710
		double[] l1c = f.search(41684.00);
		assertTrue(l1c.length == 3);
		assertTrue(StrictMath.abs(l1c[0] - 0.120724) < filetol);
		assertTrue(StrictMath.abs(l1c[1] - 0.137066) < filetol);
		assertTrue(StrictMath.abs(l1c[2] - 0.8084319) < filetol);
		
		mymisscount++; // expected cache miss, jumping around
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat line 2:
		//41685.00   .118971   .011039   .135756   .013616   .8056304   .0002710
		double[] l2 = f.search(41685.00);
		assertTrue(l2.length == 3);
		assertTrue(StrictMath.abs(l2[0] - 0.118971) < filetol);
		assertTrue(StrictMath.abs(l2[1] - 0.135756) < filetol);
		assertTrue(StrictMath.abs(l2[2] - 0.8056304) < filetol);
		
		mymisscount++; // expected cache miss, the next day
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, some line from the middle of the file...
		//48542.00   .260215   .000207   .396078   .000143   .0537336   .0000070
		double[] lm = f.search(48542.00);
		assertTrue(lm.length == 3);
		assertTrue(StrictMath.abs(lm[0] - 0.260215) < filetol);
		assertTrue(StrictMath.abs(lm[1] - 0.396078) < filetol);
		assertTrue(StrictMath.abs(lm[2] - 0.0537336) < filetol);
		
		mymisscount++; // expected cache miss
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);

		
		/* 
		 * Check linear interpolation between lines 1 & 2
		 * Double the filetol due to subtraction.
		 */
		double[] l1li = f.search( 41684.50 );
		assertTrue(l1li.length == 3);
		double m = l1li[0] - ((l1[0]+l2[0])/2.0);
		assertTrue(StrictMath.abs(m) < filetol*2.0);
		m = l1li[1] - ((l1[1]+l2[1])/2.0);
		assertTrue(StrictMath.abs(m) < filetol*2.0);
		m = l1li[2] - ((l1[2]+l2[2])/2.0);
		assertTrue(StrictMath.abs(m) < filetol*2.0);
		
		mymisscount++; // expected cache miss
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, last line of valid data
		//54054.00 -0.038151  0.021996  0.322786  0.021996  0.2564057  0.0251596
		double[] lend = f.search(54054.00);
		assertTrue(lend.length == 3);
		assertTrue(StrictMath.abs(lend[0] - -0.038151) < filetol);
		assertTrue(StrictMath.abs(lend[1] - 0.322786) < filetol);
		assertTrue(StrictMath.abs(lend[2] - 0.2564057) < filetol);
		
		mymisscount++; // expected cache miss
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, the last line of data in the file
		//54055.00, all zeros for data
		double[] lpast = f.search(54055.00);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);
		
		// no cache hit/miss change off the end of the table
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, off the table, check clamping to all zeros (54055)
		//54056.00
		lpast = f.search(54056.00);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);
		
		// no cache hit/miss change off the end of the table
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
		
		// EOP.dat, and just a little later with fractional parts
		// still off the table, clamping to zeros (54055)
		//54056.000754444445
		lpast = f.search(54056.000754444445);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);
		
		// no cache hit/miss change off the end of the table
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
	}


	/**
	 * Test method for {@link jat.spacetime.FitIERS#search(double)}.
	 * Tests caching behavior with a "cold-start" on an instance.
	 * Replicates ODTBX Mantis issue 149.
	 */
	public void testSearchM149_A() {
		
		// for tracking our expectations of the cache operation
		int myhitcount = 0;
		int mymisscount = 0;

        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = new FitIERS(filename);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		// EOP.dat, 
		//54056.00
		double[] lpast = f.search(54056.00);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);
		
		// no cache change
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
	}

	/**
	 * Test method for {@link jat.spacetime.FitIERS#search(double)}.
	 * Tests caching behavior with an instance that has already
	 * been called at least once.
	 * Replicates ODTBX Mantis issue 149.
	 */
	public void testSearchM149_B() {
		
		// for tracking our expectations of the cache operation
		int myhitcount = 0;
		int mymisscount = 0;

        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = new FitIERS(filename);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		// EOP.dat, a later line of data
		//54055.00
		double[] lpast = f.search(54055.00);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);	
		
		// no cache change
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);	
		
		// EOP.dat, 
		//54056.00
		lpast = f.search(54056.00);
		assertTrue(lpast.length == 3);
		assertTrue(lpast[0] == 0.0);
		assertTrue(lpast[1] == 0.0);
		assertTrue(lpast[2] == 0.0);
		
		// no cache change
		assertTrue(f.getCacheHits() == myhitcount);
		assertTrue(f.getCacheMisses() == mymisscount);
	}
	
	/**
	 * Test method for {@link jat.spacetime.FitIERS#getEarliestTime(double)}.
	 * Tests against a known input file with valid data.
	 * Also tests against a file with invalid data (all zeros).
	 */
	public void testgetEarliestTime() {
		
		// Test against a vaild file
        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = null;
		try {
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(String): " + t.toString());
		}
		assertNotNull(f);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		// EOP.dat line 1: yet again, to check caching
		//41684.00   .120724   .009786   .137066   .015902   .8084319   .0002710
		double et = FitIERS.getEarliestTime();
		assertTrue(StrictMath.abs(et - 41684.00) < filetol);
		
		
		// Test against a valid file with only zero data
		directory = FileUtil.getClassFilePath("jat.spacetime.unittest","FitIERSTest")+fs;
		filename = directory+"EOP_one_zeroset.dat";
		try {
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(String): " + t.toString());
		}
		assertNotNull(f);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		et = FitIERS.getEarliestTime();
		
		// EOP_one_zeroset.dat:
		// 2
		// 75001.00   .000000   .000000   .000000   .000000   .0000000   .0000000
		// 75002.00   .000000   .000000   .000000   .000000   .0000000   .0000000
		assertTrue(StrictMath.abs(et - 0.0) < filetol);
	}

	/**
	 * Test method for {@link jat.spacetime.FitIERS#getLatestTime(double)}.
	 * Tests against a known input file with valid data.
	 * Also tests against a file with invalid data (all zeros).
	 */
	public void testgetLatestTime() {
		
		// Test against a vaild file
        final String fs = FileUtil.file_separator();
		String directory;
		directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		String filename = directory+"EOP.dat";
		
		FitIERS f = null;
		try {
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(String): " + t.toString());
		}
		assertNotNull(f);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		// EOP.dat, last line of valid data
		//54054.00 -0.038151  0.021996  0.322786  0.021996  0.2564057  0.0251596
		double et = FitIERS.getLatestTime();
		assertTrue(StrictMath.abs(et - 54054.00) < filetol);

		// Test against a valid file with only zero data
		directory = FileUtil.getClassFilePath("jat.spacetime.unittest","FitIERSTest")+fs;
		filename = directory+"EOP_one_zeroset.dat";
		try {
			f = new FitIERS(filename);
		}
		catch (Throwable t) {
			fail("Caught unexpected Throwable in FitIERS(String): " + t.toString());
		}
		assertNotNull(f);
		
		if (enableCacheTest) {
			f.enableInstanceCache(); // apply the cache checks
		}
		else {
			f.disableInstanceCache(); // remove the cache from all checks
		}
		
		et = FitIERS.getLatestTime();
		
		// EOP_one_zeroset.dat:
		// 2
		// 75001.00   .000000   .000000   .000000   .000000   .0000000   .0000000
		// 75002.00   .000000   .000000   .000000   .000000   .0000000   .0000000
		assertTrue(StrictMath.abs(et - 0.0) < filetol);
	}
}
