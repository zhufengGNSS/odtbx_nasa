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

import jat.spacetime.TimeUtils;
import junit.framework.TestCase;

import org.junit.Test;

/**
 * This tests the UTC-related methods of TimeUtils.
 * @author abrown
 *
 */
public class TimeUtilsUTCTest extends TestCase {
	
	// from hardcoded TimeUtils tai_utc() before updating:
	private double[] leaps = { 41317.0, 41499.0, 41683.0, 42048.0, 
			42413.0, 42778.0, 43144.0, 43509.0, 43874.0, 44239.0,
			44786.0, 45151.0, 45516.0, 46247.0, 47161.0, 47892.0,
			48257.0, 48804.0, 49169.0, 49534.0, 50083.0, 50630.0, 
			51179.0, 53736.0};
	private int[] ls = {10, 11, 12, 13, 
			14, 15, 16, 17, 18, 19,
			20, 21, 22, 23, 24, 25, 
			26, 27, 28, 29, 30, 31,
			32, 33 };
	
	/**
	 * The number of seconds before or after a leap second to test.
	 */
	private double testrange = 33*2;

	/**
	 * Tests the tai_utc() call around the defined UTC leaps. 
	 */
	@Test
	public void testTai_utc() {
		
		// test way before the table
		assertTrue(TimeUtils.tai_utc(-1.0) == 0.0);
		
		// test after the table
		assertTrue(TimeUtils.tai_utc(leaps[leaps.length-1]) == ls[ls.length-1]);
		assertTrue(TimeUtils.tai_utc(leaps[leaps.length-1] + 10.0) == ls[ls.length-1]);
		
		// run through the table
		for (int i = 0; i < leaps.length; i++) {
			// test one second before the leap
			if (i == 0) {
				// before the table
				assertTrue(TimeUtils.tai_utc(leaps[i] - TimeUtils.sec2days) == 0.0);
			} 
			else {
				assertTrue(TimeUtils.tai_utc(leaps[i] - TimeUtils.sec2days) == ls[i-1]);
			}
			//test exactly at the leap
			assertTrue(TimeUtils.tai_utc(leaps[i]) == ls[i]);
			// test one second after the leap
			assertTrue(TimeUtils.tai_utc(leaps[i] + TimeUtils.sec2days) == ls[i]);
		}
		
	}

	/**
	 * Tests that the tai_utc_from_tai() matches tai_utc(). 
	 * Tests before the leap second table, after the leap second table, 
	 * and exercises each second before, through, and after each leap when
	 * driven from UTC.
	 * Also drives the test from TAI, and shows how tai_utc() misses the 
	 * leap seconds.
	 */
	@Test
	public void testTai_utc_from_tai() {
		
		// test way before the table
		double utc = -1.0;
		int dATutc = TimeUtils.tai_utc(utc);
		double tai = utc + dATutc*TimeUtils.sec2days;
		int dATtai = TimeUtils.tai_utc_from_tai(tai);
		assertTrue(dATutc == dATtai);
		
		// test after the table
		utc = leaps[leaps.length-1];
		dATutc = TimeUtils.tai_utc(utc);
		tai = utc + dATutc*TimeUtils.sec2days;
		dATtai = TimeUtils.tai_utc_from_tai(tai); 
		assertTrue(dATutc == dATtai);
		
		utc = leaps[leaps.length-1] + 10.0;
		dATutc = TimeUtils.tai_utc(utc);
		tai = utc + dATutc*TimeUtils.sec2days;
		dATtai = TimeUtils.tai_utc_from_tai(tai); 
		assertTrue(dATutc == dATtai);
		
		// run through the table, test at various seconds before the leap
		for (int i = 0; i < leaps.length; i++) {
			for (double j = -testrange; j < 0; j++) {
				utc = leaps[i] + j*TimeUtils.sec2days;
				dATutc = TimeUtils.tai_utc(utc);
				tai = utc + dATutc*TimeUtils.sec2days; 
				if (dATutc != 0) {
					assertTrue("i="+i+", j="+j+", "+tai+" TAI, "+utc+" UTC, dATutc="+dATutc,
							tai > utc);
				}
				dATtai = TimeUtils.tai_utc_from_tai(tai); 
				assertTrue("Mistmatch at "+utc+" UTC, dATutc="+dATutc+", dATtai="+dATtai,dATutc == dATtai);
			}
		}

		// run through the table, test exactly at the leap
		for (int i = 0; i < leaps.length; i++) {
			utc = leaps[i];
			dATutc = TimeUtils.tai_utc(utc);
			tai = utc + dATutc*TimeUtils.sec2days; 
			assertTrue(tai > utc);
			dATtai = TimeUtils.tai_utc_from_tai(tai); 
			assertTrue(dATutc == dATtai);
		}
		
		// run through the table, test at various seconds after the leap
		for (int i = 0; i < leaps.length; i++) {
			for (double j = 1; j <= testrange; j++) {
				utc = leaps[i] + j*TimeUtils.sec2days;
				dATutc = TimeUtils.tai_utc(utc);
				tai = utc + dATutc*TimeUtils.sec2days; 
				assertTrue(tai > utc);
				dATtai = TimeUtils.tai_utc_from_tai(tai); 
				assertTrue("Mistmatch at "+utc+" UTC, dATutc="+dATutc+", dATtai="+dATtai,dATutc == dATtai);
			}
		}
		
		// This time, drive linearly using TAI instead of UTC.  Show that
		// two TAI seconds are the same UTC leap second.  Step at 250msec steps.
		int i = 3;
		for (double j = -testrange; j <= testrange; j += 0.25) {
			tai = leaps[i] + j*TimeUtils.sec2days;
			dATtai = TimeUtils.tai_utc_from_tai(tai);
			utc = tai - dATtai*TimeUtils.sec2days; 
			if (dATutc != 0) {
				assertTrue("i="+i+", j="+j+", "+tai+" TAI, "+utc+" UTC, dATutc="+dATutc,
						tai > utc);
			}
			dATutc = TimeUtils.tai_utc(utc);
			//System.out.println(tai+" TAI, "+utc+" UTC, dATutc="+dATutc+", dATtai="+dATtai);
			
			// This shows how tai_utc() misses the leap second while tai_utc_from_tai()
			// accounts for the leap second:
			// 42048.000127314815 TAI, 42047.99998842592 UTC, dATutc=12, dATtai=12
			// 42048.00013888889 TAI, 42048.0 UTC, dATutc=13, dATtai=12  <-- leap second
			// 42048.00015046296 TAI, 42048.0 UTC, dATutc=13, dATtai=13  <-- leap second
			// 42048.00016203704 TAI, 42048.00001157408 UTC, dATutc=13, dATtai=13
			
			if (dATutc != dATtai) {
				// mismatch, check for a leap second
				double utcleapdiff = Math.abs(utc - leaps[i]);
				if (utcleapdiff > TimeUtils.sec2days) {
					fail("Mistmatch during leap second at "+utc+" UTC, dATutc="+dATutc+", dATtai="+dATtai);
				}
			}
			else {
				// only check the results outside of a leap second
				assertTrue("Mistmatch at "+utc+" UTC, dATutc="+dATutc+", dATtai="+dATtai,dATutc == dATtai);
			}
		}		
	}

	/**
	 * Tests gps2utc() by treating GPS as a constant bias from TAI and comparing
	 * it to a UTC result computed by tai_utc_from_tai(), tests to the microsecond.
	 * Also tests the GPS epoch of 6 Jan 1980.
	 */
	@Test
	public void testGps2utc() {
		
		// test way before the table
		double tai = -1.0;
		int dATtai = TimeUtils.tai_utc_from_tai(tai);
		double utc = tai - dATtai*TimeUtils.sec2days; 
		double gps = tai - TimeUtils.TAI_GPS*TimeUtils.sec2days;
		double testutc = TimeUtils.gps2utc(gps);
		String msg = "TAI="+tai+", UTC="+utc+", GPS="+gps+", UTCtest="+testutc+", delta="+(testutc-utc);
		//System.out.println(msg);
		// The UTC results from TAI or GPS must match.
		assertTrue(msg, Math.abs(testutc - utc) < TimeUtils.sec2days/1e6);
		
		// Run the leap second table and test a range on either side of the leap seconds.
		for (int i = 0; i < leaps.length; i++) {
			for (double j = -testrange; j <= testrange; j += 0.25) {
				tai = leaps[i] + j*TimeUtils.sec2days;
				dATtai = TimeUtils.tai_utc_from_tai(tai);
				utc = tai - dATtai*TimeUtils.sec2days; 
				gps = tai - TimeUtils.TAI_GPS*TimeUtils.sec2days;
				testutc = TimeUtils.gps2utc(gps);
				msg = "TAI="+tai+", UTC="+utc+", GPS="+gps+", UTCtest="+testutc+", delta="+(testutc-utc);
				//System.out.println(msg);
				// The UTC results from TAI or GPS must match.
				assertTrue(msg, Math.abs(testutc - utc) < TimeUtils.sec2days/1e6);
			}
		}
		
		// Finally, check the GPS epoch:
		double gpsmjd_epoch = 44244.0; // UTC MJD of GPS epoch, 6 Jan 1980
		testutc = TimeUtils.gps2utc(gpsmjd_epoch);
		msg = "GPS epoch="+gpsmjd_epoch+", testutc="+testutc+", delta="+(testutc-gpsmjd_epoch);
		// There must be no difference at the epoch date.
		assertTrue(msg, Math.abs(testutc-gpsmjd_epoch) < TimeUtils.sec2days/1e6);
		
	}
	
	/**
	 * Tests the conversion between TT to UTC and back, exercises
	 * UTCtoTT() and TTToUTC().
	 */
	@Test
	public void testTTtoUTC() {
		
		// test way before the table
		double utc = -1.0;
		double tt = TimeUtils.UTCtoTT(utc);
		double testutc = TimeUtils.TTtoUTC(tt); 
		assertTrue(Math.abs(testutc-utc) < TimeUtils.sec2days/1e6);
		
		// test after the table
		utc = leaps[leaps.length-1];
		tt = TimeUtils.UTCtoTT(utc);
		testutc = TimeUtils.TTtoUTC(tt); 
		assertTrue(Math.abs(testutc-utc) < TimeUtils.sec2days/1e6);
		
		utc = leaps[leaps.length-1] + 10.0;
		tt = TimeUtils.UTCtoTT(utc);
		testutc = TimeUtils.TTtoUTC(tt); 
		assertTrue(Math.abs(testutc-utc) < TimeUtils.sec2days/1e6);
		

		// Run the leap second table and test a range on either side of the leap seconds.
		// Drive this via TT to uncover the leap seconds.
		for (int i = 0; i < leaps.length; i++) {
			for (double j = -testrange; j <= testrange; j += 0.25) {
				//System.out.println("i="+i+", j="+j);
				tt = leaps[i] + j*TimeUtils.sec2days;
				utc = TimeUtils.TTtoUTC(tt);
				double testtt = TimeUtils.UTCtoTT(utc);
				double diff = testtt - tt;
				String msg = "UTC="+utc+", TT="+tt+", TTtest="+testtt+", delta="+diff*TimeUtils.days2sec+"s";
				//System.out.println(msg);
				if (Math.abs(diff) >= TimeUtils.sec2days/1e6) {
					// either this is a problem or during a leap second,
					// check for a leap second, and since the "first" leap second in 1972 started at 10s,
					// always compare to the number of seconds leapt at one time (compare to the prior day)
					double utcleapdiff = Math.abs(utc - leaps[i]);
					double deltaleapsecs = TimeUtils.tai_utc(utc) - TimeUtils.tai_utc(utc-1.0); 
					if (utcleapdiff > deltaleapsecs) {
						msg = "Mistmatch at: "+msg+" ("+utcleapdiff*TimeUtils.days2sec+"s from leap second)";
						System.out.println(msg);
						fail(msg);
					}
					else {
						//msg = "Caught leap at: "+msg+" ("+utcleapdiff*TimeUtils.days2sec+"s from leap second)";
						//System.out.println(msg);
					}
				}
			}
		}
	}

}
