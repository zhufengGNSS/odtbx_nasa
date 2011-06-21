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
 * Created on May 10, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package jat.test.propagator;
import jat.cm.*;
import jat.alg.integrators.*;

/**
 * @author David Gaylor
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class TwoBodyTest {

	public static void main(String[] args) {
		TwoBody tb = new TwoBody(398600.4415);
		RungeKutta8 rk8 = new RungeKutta8(5.0);
		LinePrinter lp = new LinePrinter("C:\\Documents and Settings\\David Gaylor\\My Documents\\workspace\\Jat\\jat\\test\\propagator\\25544.txt");
		lp.setThinning(600);
        // initialize the variables
        double [] x0 = new double[6];
        					
        x0[0] = 6740.83104800;
        x0[1] = -6.23521800;
        x0[2] = -462.64447700;
        x0[3] = 0.41275600;
        x0[4] = 4.78381400;
        x0[5] = 5.99847500;

        // set the final time
        double tf = 345600.0;

        // set the initial time to zero
        double t0 = 0.0;

        // integrate the equations
        rk8.integrate(t0, x0, tf, tb, lp, true);
		
	}
}
