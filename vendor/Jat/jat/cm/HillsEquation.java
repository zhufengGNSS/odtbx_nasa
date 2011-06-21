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


package jat.cm;
import java.io.*;
import jat.matvec.data.*;

/** Implements Hills Equations
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class HillsEquation {

    /** Creates a new instance of HillsEquation */
    public HillsEquation() {
    }
    public static void main(java.lang.String args[]) throws IOException {

        FileOutputStream outf = new FileOutputStream("hill3.txt");
        PrintWriter pw = new PrintWriter(outf);

        VectorN x = new VectorN(6);
        Constants c = new Constants();

        double x0 = 0.0;
        double y0 = 0.0;
        double z0 = 1.0;
        double vx0 = 0.0;
        double vy0 = 0.0;
        double vz0 = 0.0;

        TwoBody orbit = new TwoBody(6770.0, 0.0, 51.8, 0.0, 0.0, 0.0);
        double period = orbit.period();
        double w = orbit.meanMotion();

        double t = 0.0;

        while (t < 2.0*period){
            double cwt = Math.cos(w*t);
            double swt = Math.sin(w*t);
            x.x[0] = vx0*swt/w - (3.0*x0 + 2.0*vy0/w)*cwt + 4.0*x0 + 2.0*vy0/w;
            x.x[1] = (6.0*x0 + 4.0*vy0/w)*swt + 2.0*vx0*cwt/w - (6.0*w*x0+3.0*vy0)*t + y0 - 2.0*vx0/w;
            x.x[2] = z0*cwt + vz0*swt/w;
            x.x[3] = vx0*cwt + (3.0*w*x0 + 2.0*vy0)*swt;
            x.x[4] = (6.0*w*x0 + 4.0*vy0)*cwt - 2.0*vx0*swt - (6.0*w*x0 + 3.0*vy0);
            x.x[5] = -1.0*z0*w*swt + vz0*cwt;

            System.out.println(t);

            pw.print(t+"\t");
            x.print(pw);

            t = t + 1.0;
        }
        pw.close();
        outf.close();
    }
}
