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

package jat.demo.JulianDate;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

import jat.cm.*;

/**
 * @author Tobias Berthold
 *  Date        :   10-15-2003
 *  Description :   JAT Julian date demo
 * *
 */
public class CalendarDate
{
	public static void main(String argv[])
	{	Date dd;

		SimpleDateFormat sdf;
		sdf = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss aaa");
		
		// Fixed time
		//double JD = cm.juliandate(1985, 2, 17, 7, 38, 16);
		//double JD = cm.juliandate(2001, 1, 1, 12, 0, 0);
		double JD = cm.juliandate(2005, 3, 8, 14, 23, 56);
		Calendar cal = cm.Calendar(JD);
		System.out.println("JD: " + JD + "   Greg: " + sdf.format(cal.getTime()));

		// Current Time
		dd = new Date(java.lang.System.currentTimeMillis());
		//dd.setTime(java.lang.System.currentTimeMillis());
		cal.setTime(dd);
		JD = cm.juliandate(cal);
		cal = cm.Calendar(JD);
		System.out.println("JD: " + JD + "   Greg: " + sdf.format(cal.getTime()));

	}
}

// Practical Astronomy with your calculator by Peter Duffett-Smith 
// 2446113.75 is 2-17-1985 at 6am
