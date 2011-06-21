/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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
package jat.test.propagator;

import jat.alg.integrators.LinePrinter;
import jat.alg.integrators.RungeKutta8;
import jat.cm.ThreeBody;

/**
 * @author rcpage
 *
 */
public class ThreeBodyTest {

	public static void main(String[] args) {
		ThreeBodyTest test = new ThreeBodyTest();
		test.run();
	}
	public void run(){
		double G = 0.0000000000667259;
		double m1 = 1;
		double m2 = 1;
		double m3 = 1;
		ThreeBody tb = new ThreeBody(G,m1,m2,m3);
		RungeKutta8 rk8 = new RungeKutta8(5.0);
		LinePrinter lp = new LinePrinter("C:\\Code\\Jat\\jat\\test\\propagator\\25544.txt");
		lp.setThinning(600);
        // initialize the variables
        double [] x0 = new double[18];
        					
        x0[0] = 6740.83104800;
        x0[1] = -6.23521800;
        x0[2] = -462.64447700;
        x0[3] = 0.41275600;
        x0[4] = 4.78381400;
        x0[5] = 5.99847500;
        x0[6] = 0.0;
        x0[7] = 0.0;
        x0[8] = 0.0;
        x0[9] = 0.0;
        x0[10] = 0.0;
        x0[11] = 0.0;
        x0[12] = 100.0;
        x0[13] = 1.0;
        x0[14] = 100.0;
        x0[15] = 1.0;
        x0[16] = 100.0;
        x0[17] = 0.0;

        // set the final time
        double tf = 345600.0;

        // set the initial time to zero
        double t0 = 0.0;

        // integrate the equations
        rk8.integrate(t0, x0, tf, tb, lp, true);
        System.out.println("Finished");
	}
}
