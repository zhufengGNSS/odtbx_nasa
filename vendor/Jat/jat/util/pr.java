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
 /*
 *  File        :   pr.java
 *  Author      :   Tobias Berthold
 *  Date        :   11-17-2002
 *  Change      :
 *  Description :   Print utilities
 */

package jat.util;

public class pr
{

	/* Format double with Fixed width double to string  Fw.d */
    public static String fwdts (double x, int w, int d)
    {
        java.text.DecimalFormat fmt = new java.text.DecimalFormat();
        fmt.setMaximumFractionDigits(d);
        fmt.setMinimumFractionDigits(d);
        fmt.setGroupingUsed(false);
        String s = fmt.format(x);
        while (s.length() < w)
        {
            s = " " + s;
        }
        return s;
    }

}   //  end class
