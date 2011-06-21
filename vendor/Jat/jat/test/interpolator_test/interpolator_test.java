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


package jat.test.interpolator_test;


public class interpolator_test
{

    public static void main(String[] args)
    {
        // Test interpolator
        double[] data={1.,2.,3.,5.};
        double[] vari={1.,2.,3.,4.};
        double test_x;
        jat.math.Interpolator in1=new jat.math.Interpolator(vari,data);

        // too low: expect error
        test_x=in1.get_value(0.8);
        System.out.println("Interpolated value: "+test_x);

        // on the lower boundary: expect 1.0
        test_x=in1.get_value(1.0);
        System.out.println("Interpolated value: "+test_x);

        // Expect 2.3
        test_x=in1.get_value(2.3);
        System.out.println("Interpolated value: "+test_x);

        // on the upper boundary: Expect 5.0
        test_x=in1.get_value(4.0);
        System.out.println("Interpolated value: "+test_x);

        // too high: expect error
        test_x=in1.get_value(4.5);
        System.out.println("Interpolated value: "+test_x);
    }
}
