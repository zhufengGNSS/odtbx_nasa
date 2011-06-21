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
 *  File        :   ephemeris.java
 *  Author      :   Tobias Berthold
 *  Date        :   10-9-2002
 *  Change      :
 *  Description :   main program, JAT ephemeris demo
 */

package jat.demo.ephemeris;

import jat.cm.*;
import jat.util.*;
import jat.eph.*;
import jat.matvec.data.*;

/**
 * @author Tobias Berthold
 * @version 1.0
 */
public class ephemeris
{
    public static void main (String argv[])
    {
        double jd=cm.juliandate(2002, 2, 17, 12, 0, 0);
		String fs = FileUtil.file_separator();
		DE405 my_eph = new DE405(FileUtil.getClassFilePath("jat.eph","DE405")+fs+"DE405data"+fs);
        VectorN rv=my_eph.get_planet_posvel(DE405_Body.MARS, jd);

        System.out.println("The position of Mars on 10-17-2002 at 12:00pm is ");
        System.out.println("x= "+rv.get(0)+" km");
        System.out.println("y= "+rv.get(1)+" km");
        System.out.println("z= "+rv.get(2)+" km");


    }
}

