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
 *
 */

package jat.traj;

import jat.cm.*;
//import jat.plot.*;
import jat.alg.integrators.*;

/**
 * <P>
 * The TrajTest Class provides an example of how the Trajectory class can be used.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class TrajTest implements Printable
{
	public Trajectory traj = new Trajectory();


    /** Implements the Printable interface to get the data out of the propagator and pass it to the trajectory.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    public void print(double t, double [] y)
    {
    	traj.add(t, y);
    }

    /** Runs the example.
     * @param args Arguments.
     */
    public static void main(String[] args)
    {
        TrajTest x = new TrajTest();

        // create a TwoBody orbit using orbit elements
        TwoBody sat = new TwoBody(7000.0, 0.1, 0.0, 0.0, 0.0, 0.0);

        // find out the period of the orbit
        double period = sat.period();

        // set the final time = one orbit period
        double tf = period;

        // set the initial time to zero
        double t0 = 0.0;

        // propagate the orbit
        sat.propagate(t0, tf, x, true);
        
        // set the attributes
        x.traj.setTitle("Test Trajectory 1");
        x.traj.setCentralBody(CentralBody.EARTH);
        x.traj.setCoordinateSystem(CoordinateSystem.INERTIAL);
        x.traj.setDistanceUnits(DistanceUnits.KILOMETERS);
        x.traj.setEpoch(2003, 3, 27, 12, 16, 0.0);
        x.traj.setTimeUnits(TimeUnits.SECONDS);
        String[] labels = {"t","x","y","z","xdot","ydot","zdot"};
        x.traj.setLabels(labels);

        // Print out the trajectory to the screen
//        LinePrinter lp1 = new LinePrinter();        
//        x.traj.printAll(lp1);
//        lp1.close();
        
        // Print the trajectory to an ASCII file
        System.out.println("Printing to an ASCII file");
        String dir = "C:\\Temp\\";
        String file = "traj1.txt";
        LinePrinter lpout = new LinePrinter(dir + file);
        x.traj.sendToLinePrinter(lpout);
        lpout.close();
        System.out.println("traj1.txt written");
        
        // Serialize the trajectory
        x.traj.serialize(dir+"trajobj.jat");
        System.out.println("trajectory serialized");
        
        // Recover the trajectory and print all to screen       
        Trajectory traj2 = Trajectory.recover(dir+"trajobj.jat");
        traj2.setTitle("Test Trajectory 1 Recovered");
        System.out.println("Printing Recovered Trajectory");
//        LinePrinter lp2 = new LinePrinter();
//        traj2.printAll(lp2);
//        lp2.close();

       // Test 3d stuff
       
//       double[] times = traj2.timeArray();
//       double[] j3d = traj2.j3dArray();
//       int j = 0;      
//       for(int i = 0; i < times.length; i++) {
//       	   System.out.println(times[i]+" "+j3d[j]+" "+j3d[j+1]+" "+j3d[j+2]);
//       	   j = j + 3;
//       }
       	

        // Read in data from the ASCII file, set attributes, print to screen
        Trajectory traj3 = new Trajectory();
        traj3.readFromFile(dir+file);
        traj3.setTitle("Test Trajectory 1 Read From ASCII file");
        traj3.setCentralBody(CentralBody.EARTH);
        traj3.setCoordinateSystem(CoordinateSystem.INERTIAL);
        traj3.setDistanceUnits(DistanceUnits.KILOMETERS);
        traj3.setEpoch(2003, 3, 27, 12, 16, 0.0);
        traj3.setTimeUnits(TimeUnits.SECONDS);
        traj3.setLabels(labels);
        LinePrinter lp3 = new LinePrinter();
        traj3.printAll(lp3);
        lp3.close();
        
    }
}
