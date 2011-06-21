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

package jat.demo.RK78F_Test;

import jat.alg.integrators.*;
import jat.cm.*;
import java.io.*;

/**
 * <P>
 * The SimpleIntegrator Class provides an example of how the RungeKutta8 class can be used to
 * integrate a set of ODEs and plot the solution.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RK78F_Test
{

    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args) throws IOException
    {

        System.out.println("Create twobody orbit data");
        // create a Two Body orbit (elliptical)
        TwoBody orbit = new TwoBody(10000.0, 0.2, 0.0, 0.0, 0.0, 0.0);

        // create a Line Printer, print to screen
         LinePrinter lp = new LinePrinter();

        // create an RungeKuttaFehlberg78 integrator
        RungeKuttaFehlberg78 rk78f = new RungeKuttaFehlberg78();
        rk78f.setAccuracy(1.e-9);

        // initialize the state variables
        double [] x0 = orbit.randv();

        // set the final time to one orbit period
        double tf = orbit.period();

        // set the initial time to zero
        double t0 = 0.0;

        // integrate the equations
        double [] xf = rk78f.integrate(t0, x0, tf, orbit, lp, true);

        // close the Line Printer
        lp.close();

    }

}
