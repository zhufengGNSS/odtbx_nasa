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

package jat.test.path_test;


public class Jwhich
{
    public static void main (String [] args) throws Exception
    {
        System.out.println("test");

        String arg="jat.test.path_test.Jwhich";
        
        try
        {
            Class c = Class.forName (arg);
            String classRes = "/"+arg.replace ('.', '/')+".class";
            System.out.println ("Class "+args+" URL: "+c.getResource (classRes));
        }
        catch (Throwable t)
        {
            System.out.println ("Unable to locate class "+arg);
        }
    }
}

