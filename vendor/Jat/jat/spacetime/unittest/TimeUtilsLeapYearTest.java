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

import org.junit.Assert;

import jat.spacetime.TimeUtils;
import junit.framework.TestCase;

/**
 * This tests the leap year handling of TimeUtils.
 * @author abrown
 *
 */
public class TimeUtilsLeapYearTest extends TestCase {

	/**
	 * Run through and test the leap year identification
	 * from 1601 to 2201.
	 */
	public void testIsLeapYear() {
		
		// The centuries that are not leap years (2000 was a leap year):
		int common[] = {1700 , 1800 , 1900 , 2100 , 2200 };
		
		for (int y = 1601; y < 2201; y++) {
			if ( (y % 4) == 0) {
				// even div by 4
				boolean ans = TimeUtils.isLeapYear(y);
				if (!ans) {
					/* The simple rule said it should have been a leap year
					 * but the returned value said it wasn't.  Check to
					 * see if this variation from the simple y%4 rule is
					 * correct.  This should only occur at certain century 
					 * years.
					 */
					boolean failedCommonCheck = true;
					for (int c = 0; c < common.length; c++) {
						if (y == common[c]) {
							failedCommonCheck = false;
							break;
						}
					}
					if (failedCommonCheck) {
						Assert.fail("Failed on year " + y);
					}
				}
			}
			else {
				// uneven div by 4
				assertFalse(TimeUtils.isLeapYear(y));
			}	
		}	
	}

	/**
	 * Tests the leap year handling of day2doy().
	 * Checks that Mar 1st of each year is either the
	 * 59th or 60th day of the year depending on if it is
	 * a leap year.  Tests  from 1601 to 2201.
	 */
	public void testLeapYear_day2doy() {
		for (int y = 1601; y < 2201; y++) {
			if (TimeUtils.isLeapYear(y)) {
				assertTrue(TimeUtils.day2doy(y, 3, 1) == 61);
			}
			else {
				assertTrue(TimeUtils.day2doy(y, 3, 1) == 60);
			}
		}	
	}

	/**
	 * Tests the leap year handling of doy2day().
	 * Checks that Mar 1st of each year is either the
	 * 60th or 61st day of the year depending on if it is
	 * a leap year.  Tests  from 1601 to 2201.
	 */
	public void testLeapYear_doy2day() {
		for (int y = 1601; y < 2201; y++) {
			if (TimeUtils.isLeapYear(y)) {
				assertTrue(TimeUtils.doy2day(y, 61) == 1);
			}
			else {
				assertTrue(TimeUtils.doy2day(y, 60) == 1);
			}
		}	
	}

	/**
	 * Tests the leap year handling of doy2month().
	 * Checks that Mar 1st of each year is either the
	 * 60th or 61st day of the year depending on if it is
	 * a leap year.  Tests  from 1601 to 2201.
	 */
	public void testLeapYear_doy2month() {
		for (int y = 1601; y < 2201; y++) {
			if (TimeUtils.isLeapYear(y)) {
				assertTrue(TimeUtils.doy2month(y, 61) == 3);
			}
			else {
				assertTrue(TimeUtils.doy2month(y, 60) == 3);
			}
		}	
	}
	
}
