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

package jat.demo.ZeroVelCurve;

import jat.matvec.data.*;
import jat.plot.*;

/**
 * <P>
 * 
 * @author Tobias Berthold
 * @version 1.0
 */ 

public class ThreeBodyapp
{
    static SinglePlot traj_plot = new SinglePlot();
    //static double mu=0.01215;
    static double mu=0.2;
    static double L1, L2, L3, L45;

    /** Runs the example.
     * @param args Arguments.
     */
    
    ThreeBodyapp()
    {
    }

    
    public static void main(String[] args)
    {
        // Find libration points:
        // one-dimensional root finding problem
        //L1
        VectorN xinit = new VectorN(1);
        xinit.x[0]=1.-mu/2.;
        LibrationPoint lp=new LibrationPoint(mu,xinit);
        VectorN xout=lp.solveIt();
        L1=xout.x[0];
        xout.print("L1");
        
        //L2
        xinit.x[0]=.5-mu/2.;
        lp.set_initial_guess(xinit);
        xout=lp.solveIt();
        L2=xout.x[0];
        xout.print("L2");

        //L3
        xinit.x[0]=-2.*mu;
        lp.set_initial_guess(xinit);
        xout=lp.solveIt();
        L3=xout.x[0];
        xout.print("L3");

        //L4, L5
        L45=Math.sqrt(3)/2.;

        /*
        // Plot Ux along x-axis for libration points
        double x;
        for(x=-2.;x<2.;x+=0.01)
        {
            //System.out.println(x+"\t  "+Ux(mu,x));
            traj_plot.plot.addPoint(0, x,lp.Ux(mu,x),true);
        }
        */

        
        //System.out.println(ZV_func(1.,1.));
        //System.out.println(ZV_func(0.1,2.));

        // Find zero velocity curve points
        // two-dimensional root finding problem
        /*
        VectorN xinit = new VectorN(2);
        xinit.x[0]=1.0;
        xinit.x[1]=1.0;
        ZeroVelCurve zv=new ZeroVelCurve(xinit);
        VectorN xout=zv.solveIt();
        xout.print("done");
        */                
        
        // Plot the two bodies
        traj_plot.plot.addPoint(1, -mu, 0., false);
        traj_plot.plot.addPoint(1, 1.-mu, 0., false);

        // Plot libration points
        traj_plot.plot.addPoint(0, L1, 0., false);
        traj_plot.plot.addPoint(0, L2, 0., false);
        traj_plot.plot.addPoint(0, L3, 0., false);
        traj_plot.plot.addPoint(0, 0.5-mu, L45, false);
        traj_plot.plot.addPoint(0, 0.5-mu, -L45, false);


        // Plot  along y-axis
        double x=1.6;
        double y;
        double C1=3.8;
        double C2=3.55239;
        double C3=3.1973772;
        ZeroVelCurve zv=new ZeroVelCurve(mu, C1, xinit);
        for(y=-2.;y<2.;y+=0.01)
        {
            traj_plot.plot.addPoint(2, y,zv.ZV_func(x,y),true);
        }
                                                
        traj_plot.plot.setMarksStyle("dots");
        
        traj_plot.setVisible(true);
        
        traj_plot.plot.setXRange(-2.,2.);
        traj_plot.plot.setYRange(-2.,2.);
        
    }
}
